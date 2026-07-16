import requests

from vep.models.submission_form import GenomeAnnotationProvider

from core.config import WEB_METADATA_API, VEP_SUPPORT_PATH, GENOME_METADATA_API


def get_vep_support_location(genome_id: str) -> dict:
    try:
        response = requests.get(
            WEB_METADATA_API
            + "genome/"
            + genome_id
            + "/vep/file_paths"
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
        response = requests.get(
            GENOME_METADATA_API
            + "genome/"
            + genome_id
            + "/dataset/genebuild/attributes?"
            + "attribute_names=genebuild.provider_name&"
            + "attribute_names=genebuild.provider_version&"
            + "attribute_names=genebuild.last_geneset_update"
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
