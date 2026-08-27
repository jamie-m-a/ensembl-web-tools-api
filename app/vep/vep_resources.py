"""
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.


Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

import asyncio
from enum import Enum
import json
import logging
import re

from fastapi import Request, status, APIRouter, Query
from pydantic import FilePath
from requests import HTTPError
from starlette.responses import (
    JSONResponse,
    FileResponse,
    Response,
    StreamingResponse,
)
from starlette.concurrency import run_in_threadpool

from core.config import DUMP_INI, DUMP_INI_DIR, LOCAL_RESULTS_VCF
from core.error_response import response_error_handler
from core.logging import InterceptHandler
from vep.models.pipeline_model import (
    ConfigIniParams,
    VEPConfigParams,
    LaunchParams,
    PipelineParams,
    PipelineStatus,
    output_prefix_for,
)
from vep.models.submission_form import Dropdown, FormConfig
from vep.models.upload_vcf_files import (
    Streamer,
    MaxBodySizeException,
    UnsafeFileNameException,
)
from vep.utils.nextflow import launch_workflow, get_workflow_status
from vep.utils.vcf_results import get_results_from_path, stream_filtered_vcf_text
from vep.utils.tsv_export import stream_vep_tsv, flatten_vcf_lines, gzip_text_stream
from vep.utils.results_filters import parse_filters, FilterError, ResultsFilter
from vep.utils.web_metadata import get_genome_explain, get_genome_genebuild
from vep.utils.species_presets import get_species_presets
from vep.utils.spec_loader import (
    resolve_merged_spec,
    write_display_panels_sidecar,
    write_expected_columns_sidecar,
    write_spec_sidecar,
)
from vep.models.display_panels_model import to_display_panels
from vep.form_panels import get_visible_panels
# LOCAL DEV HARNESS (fork only) — see core/config.py
from vep.utils.dump_ini import dump_config_ini

logging.getLogger().handlers = [InterceptHandler()]

router = APIRouter()


class VepStatus(str, Enum):
    submitted = "SUBMITTED"
    running = "RUNNING"
    succeeded = "SUCCEEDED"
    failed = "FAILED"
    cancelled = "CANCELLED"


@router.post("/submissions", name="submit_vep")
async def submit_vep(request: Request):
    try:
        request_streamer = Streamer(request=request)
        stream_result = await request_streamer.stream()
        if not stream_result:
            raise Exception("Failed to upload VEP input files")
        vep_job_parameters = request_streamer.parameters.value.decode()
        genome_id = request_streamer.genome_id.value.decode()
        vep_job_parameters_dict = json.loads(vep_job_parameters)

        job_fields = {
            key: value
            for key, value in vep_job_parameters_dict.items()
            if key in ConfigIniParams.model_fields
        }
        # Extract the selected options from the parameters payload for validation against the spec
        options = {
            key: value
            for key, value in vep_job_parameters_dict.items()
            if key not in ConfigIniParams.model_fields
        }
        ini_parameters = ConfigIniParams(
            **job_fields, genome_id=genome_id, options=options
        )

        # Resolve the merged spec (config + parsing) for this job's assembly
        merged_spec = resolve_merged_spec(ini_parameters.assembly_name)
        # Extract the expected CSQ columns for validating the output VCF
        expected_columns = merged_spec.expected_csq_columns(ini_parameters.options)
        # Get the visible option panels for rendering the results
        display_panels = to_display_panels(
            get_visible_panels(
                species_taxonomy_id=ini_parameters.species_taxonomy_id,
                assembly_name=ini_parameters.assembly_name,
            )
        )
        # LOCAL DEV HARNESS (fork only): dump the generated config.ini and
        # return a fake id, without building launch params or contacting the
        # runner. DUMP_INI_DIR has no per-job subdirectory (unlike the real
        # outdir below), so the sidecars written here are overwritten by the
        # next submission — matching how this harness works: one manually-run
        # job at a time.
        if DUMP_INI:
            write_spec_sidecar(DUMP_INI_DIR, merged_spec)
            write_expected_columns_sidecar(DUMP_INI_DIR, expected_columns)
            write_display_panels_sidecar(DUMP_INI_DIR, display_panels)
            return {
                "submission_id": dump_config_ini(ini_parameters, merged_spec.config)
            }

        ini_file = await run_in_threadpool(
            ini_parameters.create_config_ini_file,
            request_streamer.temp_dir,
            merged_spec.config,
        )
        write_spec_sidecar(request_streamer.temp_dir, merged_spec)
        write_expected_columns_sidecar(request_streamer.temp_dir, expected_columns)
        write_display_panels_sidecar(request_streamer.temp_dir, display_panels)

        vep_job_config_parameters = VEPConfigParams(
            vcf=request_streamer.filepath,
            vep_config=ini_file.name,
            outdir=request_streamer.temp_dir,
            output_prefix=output_prefix_for(request_streamer.filename),
        )
        launch_params = LaunchParams(
            paramsText=vep_job_config_parameters, workDir=request_streamer.temp_dir
        )
        pipeline_params = PipelineParams(launch=launch_params)
        workflow_id = await run_in_threadpool(launch_workflow, pipeline_params)
        return {"submission_id": workflow_id}
    except HTTPError as e:
        try:
            msg = e.response.json()["message"]
        except Exception:
            msg = e.response.text
        logging.error(f"Upstream service error: {msg}: {e}")
        return response_error_handler(result={"status": e.response.status_code})
    except MaxBodySizeException:
        return response_error_handler(result={"status": 413})
    except UnsafeFileNameException as e:
        logging.warning(f"rejected upload file name: {e}")
        return response_error_handler(result={"status": 400})
    except ValueError as e:
        logging.warning("invalid VEP submission: %s", e)
        return response_error_handler(result={"status": 400})
    except Exception as e:
        logging.exception(f"{e.__class__.__name__}: {e}")
        return response_error_handler(result={"status": 500})


@router.get("/submissions/{submission_id}/status", name="submission_status")
async def vep_status(request: Request, submission_id: str):
    try:
        # LOCAL DEV HARNESS (fork only): there is no real pipeline run to poll
        # (submit returned a fake id), so report SUCCEEDED straight away and let
        # the results endpoint serve the local VCF.
        if DUMP_INI or LOCAL_RESULTS_VCF:
            return JSONResponse(
                content={
                    "submission_id": submission_id,
                    "status": VepStatus.succeeded.value,
                },
                headers={"Cache-Control": "no-store"},
            )

        workflow_status = await get_workflow_status(submission_id)
        submission_status = PipelineStatus(
            submission_id=submission_id, status=workflow_status
        )
        if submission_status.status == VepStatus.failed:
            logging.error(
                f"VEP submission f{submission_id} failed: f{workflow_status['workflow']['errorMessage'] or workflow_status['workflow']['errorReport']}")
        return JSONResponse(
            content=submission_status.model_dump(),
            headers={"Cache-Control": "no-store"},
        )

    except HTTPError as e:
        try:
            msg = e.response.json()["message"]
        except Exception:
            msg = e.response.text
        logging.error(f"Upstream service error: {msg}: {e}")
        return response_error_handler(result={"status": e.response.status_code})
    except Exception as e:
        logging.error(f"{e.__class__.__name__}: {e}")
        return response_error_handler(result={"status": 500})


def get_vep_results_file_path(
    input_vcf_file: str, output_prefix: str | None = None
) -> FilePath:
    input_vcf_path = FilePath(input_vcf_file)
    result_name = (
        f"{output_prefix}_VEP.vcf.gz"
        if output_prefix
        else input_vcf_path.stem + "_VEP.vcf.gz"
    )
    return input_vcf_path.with_name(result_name)


def _gzip_download_response(text_stream, filename: str) -> StreamingResponse:
    return StreamingResponse(
        gzip_text_stream(text_stream),
        media_type="application/gzip",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


def _results_download_response(
    results_path: FilePath,
    output_format: str,
    active_filters: list[ResultsFilter] | None = None,
) -> FileResponse | StreamingResponse:
    is_table = output_format in ("tsv", "txt", "table")
    table_extension = "tsv" if output_format == "tsv" else "txt"
    base = re.sub(r"\.vcf(\.gz)?$", "", results_path.name) or "vep_results"
    if active_filters:
        # Filtered download (include only CSQ entries/records passing the filters)
        vcf_text = stream_filtered_vcf_text(results_path, active_filters)
        if is_table:
            # Filtered download in tabular format (send over gzipped stream)
            return _gzip_download_response(
                flatten_vcf_lines(vcf_text),
                f"{base}_filtered.{table_extension}.gz",
            )
        return _gzip_download_response(vcf_text, f"{base}_filtered.vcf.gz")
    if is_table:
        # Full results in tabular format
        return _gzip_download_response(
            stream_vep_tsv(results_path), f"{base}.{table_extension}.gz"
        )
    return FileResponse(
        # Full results in VCF format
        results_path,
        media_type="application/gzip",
        filename=results_path.name,
    )


@router.get("/submissions/{submission_id}/download", name="download_results")
async def download_results(
    request: Request,
    submission_id: str,
    format: str = "vcf",
    filters: str | None = None,
):
    # Optional server-side filtering.
    # `filters` is a JSON array of query-builder conditions.
    try:
        active_filters = parse_filters(filters)
    except FilterError as exc:
        return JSONResponse(
            content={"details": f"Invalid filters: {exc}"},
            status_code=status.HTTP_400_BAD_REQUEST,
        )
    try:
        # LOCAL DEV HARNESS (fork only): serve the VEP output VCF on disk
        # directly, bypassing the Seqera status lookup.
        if LOCAL_RESULTS_VCF:
            return _results_download_response(
                FilePath(LOCAL_RESULTS_VCF), format, active_filters
            )

        workflow_status = await get_workflow_status(submission_id)
        submission_status = PipelineStatus(
            submission_id=submission_id, status=workflow_status
        )
        if submission_status.status == VepStatus.succeeded:
            input_vcf_file = workflow_status["workflow"]["params"]["input"]
            output_prefix = workflow_status["workflow"]["params"].get(
                "output_prefix"
            )
            results_file_path = get_vep_results_file_path(
                input_vcf_file, output_prefix
            )
            if results_file_path.exists():
                return _results_download_response(
                    results_file_path, format, active_filters
                )
            else:
                response_msg = {
                    "details": f"A submission with id {submission_id} succeeded but could not find output file",
                }
                return JSONResponse(
                    content=response_msg, status_code=status.HTTP_404_NOT_FOUND
                )
        else:
            response_msg = {
                "details": f"A submission with id {submission_id} is not yet finished",
            }
            return JSONResponse(
                content=response_msg, status_code=status.HTTP_404_NOT_FOUND
            )

    except FilterError as exc:
        # Compiling the filters against the file's CSQ header failed (e.g. a filter
        # references a column this output doesn't carry) — a client error.
        return JSONResponse(
            content={"details": f"Invalid filters: {exc}"},
            status_code=status.HTTP_400_BAD_REQUEST,
        )
    except HTTPError as e:
        if e.response.status_code in [403, 400]:

            response_msg = {
                "status_code": status.HTTP_404_NOT_FOUND,
                "details": f"A submission with id {submission_id} was not found",
            }
            return JSONResponse(
                content=response_msg, status_code=status.HTTP_404_NOT_FOUND
            )
        else:
            logging.error(f"Upstream service error: {e}")
        return response_error_handler(result={"status": e.response.status_code})
    except Exception as e:
        logging.error(f"{e.__class__.__name__}: {e}")
        return response_error_handler(result={"status": 500})


def _results_response(**kwargs) -> Response:
    """Return a serialized JSONResponse for the results page.
    Pre-serialized to bytes with `model_dump_json` for performance.
    `by_alias=True` is needed to map `source` field in the spec model to `from` in the payload.
    """
    payload = get_results_from_path(**kwargs)
    return Response(
        content=payload.model_dump_json(by_alias=True),
        media_type="application/json",
    )


@router.get("/submissions/{submission_id}/results", name="view_results")
async def fetch_results(
    request: Request,
    submission_id: str,
    page: int = Query(..., ge=1),
    per_page: int = Query(..., ge=1, le=500),
    filters: str | None = None,
):
    results_file_path = None
    try:
        # Optional server-side filtering: `filters` is a JSON array of query-builder
        # conditions. Malformed input is a client error (400), not a 500.
        try:
            active_filters = parse_filters(filters)
        except FilterError as exc:
            return JSONResponse(
                content={"details": f"Invalid filters: {exc}"},
                status_code=status.HTTP_400_BAD_REQUEST,
            )

        def _results(**kwargs):
            """
            `parse_filters` only checks the filter payload shape.
            `_results_response` compiles the filters to check that the requested 
            columns and operators exist in the result VCF.
            """
            try:
                return _results_response(**kwargs)
            except FilterError as exc:
                return JSONResponse(
                    content={"details": f"Invalid filters: {exc}"},
                    status_code=status.HTTP_400_BAD_REQUEST,
                )

        # LOCAL DEV HARNESS (fork only): parse the VEP output VCF on disk
        # directly, bypassing the Seqera status lookup. Reading and (for a
        # filtered request) scanning it is blocking and CPU-bound, so run it in
        # a worker thread rather than stalling the event loop.
        if LOCAL_RESULTS_VCF:
            return await run_in_threadpool(
                _results,
                vcf_path=FilePath(LOCAL_RESULTS_VCF),
                page=page,
                page_size=per_page,
                filters=active_filters,
            )

        workflow_status = await get_workflow_status(submission_id)
        submission_status = PipelineStatus(
            submission_id=submission_id, status=workflow_status
        )
        if submission_status.status == VepStatus.succeeded:
            input_vcf_file = workflow_status["workflow"]["params"]["input"]
            output_prefix = workflow_status["workflow"]["params"].get(
                "output_prefix"
            )
            results_file_path = get_vep_results_file_path(
                input_vcf_file, output_prefix
            )
            if results_file_path.exists():
                return await run_in_threadpool(
                    _results,
                    vcf_path=results_file_path,
                    page=page,
                    page_size=per_page,
                    filters=active_filters,
                )
            else:
                response_msg = {
                    "details": f"A submission with id {submission_id} succeeded but could not find output file",
                }
                return JSONResponse(
                    content=response_msg, status_code=status.HTTP_404_NOT_FOUND
                )
        else:
            response_msg = {
                "details": f"A submission with id {submission_id} is not yet finished",
            }
            return JSONResponse(
                content=response_msg, status_code=status.HTTP_404_NOT_FOUND
            )
    except HTTPError as e:
        if e.response.status_code in [403, 400]:
            response_msg = json.dumps(
                {
                    "status_code": status.HTTP_404_NOT_FOUND,
                    "details": f"A submission with id {submission_id} was not found",
                }
            )
            return JSONResponse(
                content=response_msg, status_code=status.HTTP_404_NOT_FOUND
            )
        else:
            logging.error(f"Upstream service error: {e}")
        return response_error_handler(result={"status": e.response.status_code})
    except Exception as e:
        logging.error(f"{e.__class__.__name__}: {e} (VCF: {results_file_path})")
        return response_error_handler(result={"status": 500})


@router.get("/form_config/{genome_id}", name="get_form_config")
async def get_form_config(
    request: Request,
    genome_id: str,
):
    try:
        attributes, genome = await asyncio.gather(
            get_genome_genebuild(genome_id), get_genome_explain(genome_id)
        )
        species_taxonomy_id = genome.get("species_taxonomy_id")
        assembly_name = (genome.get("assembly") or {}).get("name")
        if not species_taxonomy_id or not assembly_name:
            raise ValueError(
                "get_form_config(): unexpected metadata API explain payload "
                f"for {genome_id}: missing species_taxonomy_id or assembly.name"
            )

        annotation_provider_name = attributes.get("genebuild.provider_name", "")
        annotation_version = attributes.get("genebuild.provider_version", "")
        last_updated_date = attributes.get("genebuild.last_geneset_update", "")

        if (annotation_version or last_updated_date):
            label = f"{annotation_provider_name} {annotation_version or last_updated_date}"
            value = f"{annotation_provider_name}_{annotation_version or last_updated_date}"
        else:
            label = f"{annotation_provider_name}"
            value = f"{annotation_provider_name}"

        options = [{
            "label": label,
            "value": value
        }]

        default_option = options[0]
        transcript_set = Dropdown(
            label="Transcript set",
            options=options,
            default_value=default_option["value"],
        )

        form_config = FormConfig(transcript_set=transcript_set)
        # Panels/options to show for this genome's canonical species/assembly.
        return {
            "parameters": form_config,
            "panels": get_visible_panels(
                attributes,
                species_taxonomy_id=species_taxonomy_id,
                assembly_name=assembly_name,
            ),
        }

    except HTTPError as e:
        if e.response.status_code == 404:
            response_msg = json.dumps(
                {
                    "status_code": status.HTTP_404_NOT_FOUND,
                    "details": f"genome id {genome_id} not found",
                }
            )
            return JSONResponse(
                content=response_msg, status_code=status.HTTP_404_NOT_FOUND
            )
        else:
            logging.error(f"Upstream service error: {e}")
        return response_error_handler(result={"status": e.response.status_code})
    except Exception as e:
        logging.error(f"{e.__class__.__name__}: {e}")
        return response_error_handler(result={"status": 500})


@router.get("/species_presets", name="get_species_presets")
async def species_presets(request: Request):
    """Species presets for the input form quick-select buttons."""
    try:
        return await get_species_presets()
    except Exception as e:
        logging.exception(f"{e.__class__.__name__}: {e}")
        return response_error_handler(result={"status": 500})
