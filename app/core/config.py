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

import logging
import sys

from loguru import logger
from starlette.config import Config
from starlette.datastructures import CommaSeparatedStrings

from .logging import InterceptHandler
import json

VERSION = "0.0.0"
API_PREFIX = "/api/tools"

config = Config(".env")
DEBUG: bool = config("DEBUG", cast=bool, default=False)
TRUST_ENV: bool = config("TRUST_ENV", cast=bool, default=True)
PROJECT_NAME: str = config("PROJECT_NAME", default="Ensembl Web Tools API")
ALLOWED_HOSTS: list[str] = config(
    "ALLOWED_HOSTS",
    cast=CommaSeparatedStrings,
    default="*",
)
# In the container the data dir is mounted at /data; for local dev fall back to
# the repo's own data/ directory so the app boots without that mount.
import os as _os

_blast_config_path = "/data/blast_config.json"
if not _os.path.exists(_blast_config_path):
    _blast_config_path = _os.path.join(
        _os.path.dirname(__file__), "..", "..", "data", "blast_config.json"
    )
with open(_blast_config_path) as f:
    BLAST_CONFIG = json.load(f)


# logging configuration
logging.basicConfig(level=logging.DEBUG)
LOGGING_LEVEL = logging.DEBUG if DEBUG else logging.INFO
LOGGERS = ("uvicorn.asgi", "uvicorn.access")
logging.getLogger().handlers = [InterceptHandler()]
for logger_name in LOGGERS:
    logging_logger = logging.getLogger(logger_name)
    logging_logger.handlers = [InterceptHandler(level=LOGGING_LEVEL)]

logger.configure(handlers=[{"sink": sys.stderr, "level": LOGGING_LEVEL}])

import os

# Dump-ini switch (dev/testing, temporary): when enabled, a submission builds
# the VEP config.ini and writes it to DUMP_INI_DIR instead of launching the
# pipeline, returning a fake submission id. Used to inspect the form -> ini
# stage end to end. DUMP_INI_DIR defaults to the shared repo data/output dir
# (sibling of this repo), overridable via the env var.
DUMP_INI: bool = config("DUMP_INI", cast=bool, default=False)
_default_dump_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "data", "output")
)
DUMP_INI_DIR: str = config("DUMP_INI_DIR", default=_default_dump_dir)

# Local results mode (dev/testing, temporary): when LOCAL_RESULTS_VCF is set to
# a VEP output VCF path, the results endpoint parses that file directly instead
# of resolving the submission via Seqera. Off by default; discrete and easily
# removed. Example: LOCAL_RESULTS_VCF=/path/to/data/output/output.vcf.gz
LOCAL_RESULTS_VCF: str = config("LOCAL_RESULTS_VCF", default="")

# Nextflow Configurations. Defaults are empty so the app can run in mock mode
# without Seqera credentials; the real runner requires these to be set.
NF_TOKEN = config("NF_TOKEN", default="")
NF_COMPUTE_ENV_ID = config("NF_COMPUTE_ENV_ID", default="")
NF_PIPELINE_URL = config("NF_PIPELINE_URL", default="")
NF_WORK_DIR = config("NF_WORK_DIR", default="")
SEQERA_API = config("SEQERA_API", default="")
NF_WORKSPACE_ID = config("NF_WORKSPACE_ID", default="")

WEB_METADATA_API = config(
    "WEB_METADATA_API", default="https://beta.ensembl.org/api/metadata/"
)
VEP_SUPPORT_PATH = config("VEP_SUPPORT_PATH", default="/tmpdir")

# Genome-metadata API base used for form config (get_genome_metadata). Defaults
# to staging so it matches the species-search source; flip GENOME_METADATA_LIVE
# on to use the live API. Either URL can be overridden explicitly.
GENOME_METADATA_API_STAGING = config(
    "GENOME_METADATA_API_STAGING",
    default="https://staging-2020.ensembl.org/api/metadata/",
)
GENOME_METADATA_API_LIVE = config(
    "GENOME_METADATA_API_LIVE",
    default="https://beta.ensembl.org/api/metadata/",
)
GENOME_METADATA_LIVE: bool = config(
    "GENOME_METADATA_LIVE", cast=bool, default=False
)
GENOME_METADATA_API = (
    GENOME_METADATA_API_LIVE if GENOME_METADATA_LIVE else GENOME_METADATA_API_STAGING
)
