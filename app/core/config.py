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
import os
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

_blast_config_path = "/data/blast_config.json"
if not os.path.exists(_blast_config_path):
    _blast_config_path = os.path.join(
        os.path.dirname(__file__), "..", "..", "data", "blast_config.json"
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

# Nextflow configuration. Supplied by envvars (k8s secret)
NF_TOKEN = config("NF_TOKEN", default="")
NF_COMPUTE_ENV_ID = config("NF_COMPUTE_ENV_ID", default="")
NF_PIPELINE_URL = config("NF_PIPELINE_URL", default="")
NF_WORK_DIR = config("NF_WORK_DIR", default="")
SEQERA_API = config("SEQERA_API", default="")
NF_WORKSPACE_ID = config("NF_WORKSPACE_ID", default="")

WEB_METADATA_API = config(
    "WEB_METADATA_API", default="https://www.ensembl.org/api/metadata/"
)
VEP_SUPPORT_PATH_ROOT = config("VEP_SUPPORT_PATH", default="/tmpdir")
VEP_SUPPORT_PATH = os.path.join(VEP_SUPPORT_PATH_ROOT, "organisms")
VEP_PLUGIN_DATA_PATH = os.path.join(VEP_SUPPORT_PATH_ROOT, "vep-plugins-data")

# ---------------------------------------------------------------------------
# LOCAL DEV HARNESS — not for upstream.
#
# Upstream removed these in 68898ea ("Remove dev mode, cleanup"). They are kept
# on this fork only so the frontend can be developed against a local backend
# without Seqera credentials or a pipeline run. Nothing in production reads
# them: both are off unless explicitly set, and the branches they gate are
# guarded on the values below rather than on DEBUG.
#
# Drop this block, `vep/utils/dump_ini.py`, the `dev` service in
# docker-compose.yaml and the branches in `vep_resources.py` if upstream ever
# ships a supported local-development path.
# ---------------------------------------------------------------------------

# When enabled, a submission builds the VEP config.ini, writes it to
# DUMP_INI_DIR and returns a fake submission id, instead of launching the
# pipeline. Lets the form -> ini stage be inspected end to end.
DUMP_INI: bool = config("DUMP_INI", cast=bool, default=False)
_default_dump_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "data", "output")
)
DUMP_INI_DIR: str = config("DUMP_INI_DIR", default=_default_dump_dir)

# When set to a VEP output VCF path, the results and download endpoints parse
# that file directly instead of resolving the submission via Seqera.
LOCAL_RESULTS_VCF: str = config("LOCAL_RESULTS_VCF", default="")
