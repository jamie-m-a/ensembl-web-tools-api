"""TEMPORARY dev helper.

Dumps the generated VEP ``config.ini`` to disk instead of launching the
Nextflow/Seqera pipeline, so the "form -> ini" stage can be inspected end to
end while the stages are being wired up. This is intentionally kept separate
from the real submission path and is expected to be removed once the pipeline
integration is exercised for real.
"""

import datetime
import logging
import os
import uuid

from core.config import DUMP_INI_DIR
from vep.models.config_spec_model import ConfigSpec
from vep.models.pipeline_model import ConfigIniParams


def dump_config_ini(ini_parameters: ConfigIniParams, config_spec: ConfigSpec) -> str:
    """Write the config.ini built from ``ini_parameters`` (and the config half of
    the job's pinned spec) into ``DUMP_INI_DIR`` under a unique, timestamped
    filename and return a fake submission id.

    Does not contact the pipeline runner.
    """
    os.makedirs(DUMP_INI_DIR, exist_ok=True)
    submission_id = f"dump-{uuid.uuid4().hex[:12]}"

    # create_config_ini_file always writes "config.ini" into the given dir;
    # rename it to a unique file so successive dumps don't clobber each other.
    ini_file = ini_parameters.create_config_ini_file(DUMP_INI_DIR, config_spec)
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    destination = os.path.join(
        DUMP_INI_DIR, f"config-{timestamp}-{submission_id}.ini"
    )
    os.replace(ini_file.name, destination)

    logging.info("DUMP_INI: wrote VEP config to %s", destination)
    return submission_id
