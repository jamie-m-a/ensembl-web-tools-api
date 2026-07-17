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

from enum import Enum
import json
import logging
import re

from fastapi import Request, status, APIRouter
from pydantic import FilePath
from requests import HTTPError
from starlette.responses import JSONResponse, FileResponse, StreamingResponse

from core.config import DUMP_INI, DUMP_INI_DIR, LOCAL_RESULTS_VCF
from core.error_response import response_error_handler
from core.logging import InterceptHandler
from vep.models.pipeline_model import (
    ConfigIniParams,
    VEPConfigParams,
    LaunchParams,
    PipelineParams,
    PipelineStatus,
)
from vep.models.submission_form import Dropdown, FormConfig
from vep.models.upload_vcf_files import Streamer, MaxBodySizeException
from vep.utils.nextflow import launch_workflow, get_workflow_status
from vep.utils.dump_ini import dump_config_ini
from vep.utils.vcf_results import get_results_from_path
from vep.utils.tsv_export import stream_vep_tsv, gzip_text_stream
from vep.utils.results_filters import parse_filters, FilterError
from vep.utils.web_metadata import get_genome_metadata
from vep.utils.spec_loader import resolve_spec, write_spec_sidecar
from vep.form_panels import get_visible_panels

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
        ini_parameters = ConfigIniParams(**vep_job_parameters_dict, genome_id=genome_id)

        # Resolve and pin the parsing spec for this job's assembly now, at
        # submission, rather than waiting for results: it must fail here (while
        # the user is present and nothing has run yet) rather than after the
        # pipeline completes, and it must be the spec used to build the options
        # this submission is based on, not whatever the "current" one is by the
        # time results are parsed (see spec_loader.resolve_spec).
        parsing_spec = resolve_spec(ini_parameters.assembly_name)

        if DUMP_INI:
            # Temporary: dump the generated config.ini to disk and return a fake
            # id, without building launch params or contacting the runner.
            # DUMP_INI_DIR has no per-job subdirectory (unlike the real outdir
            # below), so the sidecar written here is overwritten by the next
            # submission — matching how this dev harness already works: one
            # manually-run job at a time (see write_spec_sidecar).
            write_spec_sidecar(DUMP_INI_DIR, parsing_spec)
            return {"submission_id": dump_config_ini(ini_parameters)}
        ini_file = ini_parameters.create_config_ini_file(request_streamer.temp_dir)
        write_spec_sidecar(request_streamer.temp_dir, parsing_spec)

        vep_job_config_parameters = VEPConfigParams(
            vcf=request_streamer.filepath,
            vep_config=ini_file.name,
            outdir=request_streamer.temp_dir,
        )
        launch_params = LaunchParams(
            paramsText=vep_job_config_parameters, workDir=request_streamer.temp_dir
        )
        pipeline_params = PipelineParams(launch=launch_params)
        if stream_result:
            workflow_id = launch_workflow(pipeline_params)
            return {"submission_id": workflow_id}
        else:
            raise Exception("Failed to upload VEP input files")
    except HTTPError as e:
        try:
            msg = e.response.json()["message"]
        except Exception:
            msg = e.response.text
        logging.error(f"Upstream service error: {msg}: {e}")
        return response_error_handler(result={"status": e.response.status_code})
    except MaxBodySizeException:
        return response_error_handler(result={"status": 413})
    except Exception as e:
        logging.exception(f"{e.__class__.__name__}: {e}")
        return response_error_handler(result={"status": 500})


@router.get("/submissions/{submission_id}/status", name="submission_status")
async def vep_status(request: Request, submission_id: str):
    try:
        # Dev short-circuit: in DUMP_INI / LOCAL_RESULTS_VCF mode there is no real
        # pipeline run to poll (the submission returned a fake id), so report
        # SUCCEEDED straight away and let the results endpoint serve the local
        # VCF. TEMPORARY: paired with DUMP_INI / LOCAL_RESULTS_VCF.
        if DUMP_INI or LOCAL_RESULTS_VCF:
            return JSONResponse(
                content={
                    "submission_id": submission_id,
                    "status": VepStatus.succeeded.value,
                }
            )
        workflow_status = await get_workflow_status(submission_id)
        submission_status = PipelineStatus(
            submission_id=submission_id, status=workflow_status
        )
        if submission_status.status == VepStatus.failed:
            logging.error(
                f"VEP submission f{submission_id} failed: f{workflow_status['workflow']['errorMessage'] or workflow_status['workflow']['errorReport']}")
        return JSONResponse(content=submission_status.model_dump())

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


def get_vep_results_file_path(input_vcf_file: str) -> FilePath:
    input_vcf_path = FilePath(input_vcf_file)
    vep_results_file = input_vcf_path.with_name(
        input_vcf_path.stem + "_VEP.vcf.gz"
    )
    return vep_results_file


# Build the download response for a resolved results VCF path. `format=vcf`
# (default) serves the raw VCF; `format=tsv` streams the flattened, fully
# expanded "columnar" table (spreadsheet-friendly), gzip-compressed.
def _results_download_response(results_path: FilePath, output_format: str):
    if output_format in ("tsv", "txt", "table"):
        base = re.sub(r"\.vcf(\.gz)?$", "", results_path.name) or "vep_results"
        # Compress the table before it's sent (plain gzip, for compatibility).
        # It's served as the payload (application/gzip, .txt.gz filename), not via
        # Content-Encoding, so the browser saves the compressed file rather than
        # transparently decompressing it.
        return StreamingResponse(
            gzip_text_stream(stream_vep_tsv(results_path)),
            media_type="application/gzip",
            headers={"Content-Disposition": f'attachment; filename="{base}.txt.gz"'},
        )
    return FileResponse(
        results_path,
        media_type="application/gzip",
        filename=results_path.name,
    )


@router.get("/submissions/{submission_id}/download", name="download_results")
async def download_results(
    request: Request, submission_id: str, format: str = "vcf"
):
    try:
        # Temporary local-results mode: serve the VEP output VCF on disk directly,
        # bypassing the Seqera status lookup. Enabled by setting LOCAL_RESULTS_VCF
        # (the same file the results view parses). Discrete and easily removed.
        if LOCAL_RESULTS_VCF:
            return _results_download_response(FilePath(LOCAL_RESULTS_VCF), format)
        workflow_status = await get_workflow_status(submission_id)
        submission_status = PipelineStatus(
            submission_id=submission_id, status=workflow_status
        )
        if submission_status.status == VepStatus.succeeded:
            input_vcf_file = workflow_status["workflow"]["params"]["input"]
            results_file_path = get_vep_results_file_path(input_vcf_file)
            if results_file_path.exists():
                return _results_download_response(results_file_path, format)
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


@router.get("/submissions/{submission_id}/results", name="view_results")
async def fetch_results(
    request: Request,
    submission_id: str,
    page: int,
    per_page: int,
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
        # Temporary local-results mode: parse a VEP output VCF on disk directly,
        # bypassing the Seqera status lookup. Enabled by setting LOCAL_RESULTS_VCF.
        if LOCAL_RESULTS_VCF:
            return get_results_from_path(
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
            results_file_path = get_vep_results_file_path(input_vcf_file)
            if results_file_path.exists():
                return get_results_from_path(
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
    species_taxonomy_id: str | None = None,
    assembly_name: str | None = None,
):
    try:
        attributes = await get_genome_metadata(genome_id)
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
        # Panels/options to show on the input form for this genome. Currently the
        # always-visible set for every species; the genome `attributes` are the
        # hook for species-conditional visibility later.
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
