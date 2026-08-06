import requests
from starlette.concurrency import run_in_threadpool

from vep.models.pipeline_model import PipelineParams
from core.config import (
    NF_TOKEN,
    SEQERA_API,
    NF_WORKSPACE_ID,
    SEQERA_LAUNCH_TIMEOUT,
    SEQERA_STATUS_TIMEOUT,
)


def launch_workflow(pipeline_params: PipelineParams):
    """Synchronous, and blocking. Async callers must reach it through a
    threadpool (see the submit route) so it does not stall the event loop."""
    try:
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {NF_TOKEN}",
        }
        params = {"workspaceId": NF_WORKSPACE_ID}
        SEQERA_WORKFLOW_LAUNCH_URL = SEQERA_API + "/workflow/launch"
        payload = pipeline_params.model_dump()
        response = requests.post(
            SEQERA_WORKFLOW_LAUNCH_URL,
            params=params,
            headers=headers,
            json=payload,
            timeout=SEQERA_LAUNCH_TIMEOUT,
        )
        response.raise_for_status()
        response_json = response.json()
        return response_json["workflowId"]
    except KeyError as e:
        e.args = (
            f"launch_workflow(): unexpected payload from Seqera: f{response.text}",
            *e.args,
        )
        raise
    except requests.HTTPError as e:
        e.args = (
            "launch_workflow(): error response from Seqera:",
            *e.args,
        )
        raise
    except (requests.ConnectionError, requests.Timeout) as e:
        e.args = (
            "launch_workflow(): network error while connecting to Seqera:",
            *e.args,
        )
        raise
    except Exception as e:
        e.args = (f"{type(e).__name__} in launch_workflow():", *e.args)
        raise


async def get_workflow_status(submission_id):
    try:
        _headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {NF_TOKEN}",
        }
        _seqera_workflow_status_url = f"{SEQERA_API}/workflow/{submission_id}"
        params = {"workspaceId": NF_WORKSPACE_ID}
        # Synchronous `requests` inside a coroutine would block the event loop.
        # This is the hot one: it runs every 15s for every active submission.
        response = await run_in_threadpool(
            requests.get,
            _seqera_workflow_status_url,
            params=params,
            headers=_headers,
            timeout=SEQERA_STATUS_TIMEOUT,
        )

        response.raise_for_status()
        response_json = response.json()
        return response_json
    except KeyError as e:
        e.args = (
            f"get_workflow_status(): unexpected payload from Seqera: {response.text}",
            *e.args,
        )
        raise
    except requests.HTTPError as e:
        e.args = (
            "get_workflow_status(): error response from Seqera:",
            *e.args,
        )
        raise
    except (requests.ConnectionError, requests.Timeout) as e:
        e.args = (
            "get_workflow_status(): network error while connecting to Seqera:",
            *e.args,
        )
        raise
    except Exception as e:
        e.args = (f"{type(e).__name__} in get_workflow_status():", *e.args)
        raise
