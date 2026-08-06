import requests
from starlette.concurrency import run_in_threadpool

from vep.models.submission_form import GenomeAnnotationProvider

from core.config import (
    WEB_METADATA_API,
    VEP_SUPPORT_PATH,
    GENOME_METADATA_API,
    METADATA_TIMEOUT,
)


def get_vep_support_location(genome_id: str) -> dict:
    """Synchronous, and blocking: it is called from `create_config_ini_file`,
    itself sync. Async callers must reach it through a threadpool (see the
    submit route) so the blocking request does not stall the event loop."""
    try:
        response = requests.get(
            WEB_METADATA_API
            + "genome/"
            + genome_id
            + "/vep/file_paths",
            timeout=METADATA_TIMEOUT,
        )
        response.raise_for_status()
        data = response.json()
        return {
            "faa_location": f"{VEP_SUPPORT_PATH}{data['faa_location']}",
            "gff_location": f"{VEP_SUPPORT_PATH}{data['gff_location']}",
        }
    except KeyError as e:
        e.args = (
            f"get_vep_support_location(): unexpected metadata API payload for f{genome_id}:",
            *e.args,
        )
        raise
    except requests.HTTPError as e:
        e.args = (
            f"get_vep_support_location(): error response from metadata API for f{genome_id}:",
            *e.args,
        )
        raise
    except Exception as e:
        e.args = (f"{type(e).__name__} in get_vep_support_location():", *e.args)
        raise


async def get_genome_metadata(genome_id: str) -> GenomeAnnotationProvider:
    try:
        # `requests` is synchronous, so calling it directly from this coroutine
        # would block the event loop for the whole round trip — stalling every
        # other in-flight request, not just this one.
        response = await run_in_threadpool(
            requests.get,
            GENOME_METADATA_API
            + "genome/"
            + genome_id
            + "/dataset/genebuild/attributes?"
            + "attribute_names=genebuild.provider_name&"
            + "attribute_names=genebuild.provider_version&"
            + "attribute_names=genebuild.last_geneset_update",
            timeout=METADATA_TIMEOUT,
        )
        response.raise_for_status()
        attributes = {}
        for attribute in response.json()["attributes"]:
            name = attribute["name"]
            value = attribute["value"]
            attributes[name] = value
        return attributes
    except KeyError as e:
        e.args = (
            f"get_genome_metadata(): unexpected metadata API payload for f{genome_id}:",
            *e.args,
        )
        raise
    except requests.HTTPError as e:
        e.args = (
            f"get_genome_metadata(): error response from metadata API for f{genome_id}:",
            *e.args,
        )
        raise
    except Exception as e:
        e.args = (f"{type(e).__name__} in get_genome_metadata():", *e.args)
        raise
