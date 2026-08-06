"""Tests for the outbound HTTP calls to Seqera and the metadata APIs.

Two properties, neither of which the code enforced before:

1. Every request carries an explicit timeout. Without one, `requests` waits
   forever, and a hung upstream holds the worker indefinitely — worst for the
   status poll, which runs every 15s for every active submission.
2. The blocking calls do not run on the event loop. `requests` is synchronous,
   so awaiting it directly inside a coroutine stalls *every* in-flight request
   for the duration of the round trip, not just the one that made the call.
"""

import asyncio
import threading

import pytest

from vep.utils import nextflow, web_metadata


class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload
        self.text = "fake response"

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


class _Recorder:
    """Stands in for requests.get/post, capturing how it was called and which
    thread it ran on."""

    def __init__(self, payload):
        self._payload = payload
        self.calls = []

    def __call__(self, *args, **kwargs):
        self.calls.append(
            {"args": args, "kwargs": kwargs, "thread_id": threading.get_ident()}
        )
        return _FakeResponse(self._payload)

    @property
    def timeout(self):
        return self.calls[0]["kwargs"].get("timeout")


SUPPORT_PATHS = {"faa_location": "/f.faa", "gff_location": "/g.gff"}
GENOME_ATTRS = {"attributes": [{"name": "genebuild.provider_name", "value": "Ensembl"}]}
WORKFLOW_LAUNCH = {"workflowId": "wf-123"}
WORKFLOW_STATUS = {"status": "SUCCEEDED"}


class _StubPipelineParams:
    """launch_workflow only ever calls .model_dump() on its argument, so the
    full PipelineParams tree (which validates real paths) is not needed here."""

    def model_dump(self):
        return {"launch": {}}


# --- 1. every outbound call sets a timeout ---------------------------------


def _assert_connect_read_pair(timeout):
    assert timeout is not None, "request was made with no timeout"
    assert isinstance(timeout, tuple) and len(timeout) == 2, (
        f"expected a (connect, read) pair, got {timeout!r}"
    )
    assert all(isinstance(v, float) and v > 0 for v in timeout)


def test_get_vep_support_location_sets_a_timeout(monkeypatch):
    recorder = _Recorder(SUPPORT_PATHS)
    monkeypatch.setattr(web_metadata.requests, "get", recorder)

    web_metadata.get_vep_support_location("genome-1")

    _assert_connect_read_pair(recorder.timeout)


def test_get_genome_metadata_sets_a_timeout(monkeypatch):
    recorder = _Recorder(GENOME_ATTRS)
    monkeypatch.setattr(web_metadata.requests, "get", recorder)

    asyncio.run(web_metadata.get_genome_metadata("genome-1"))

    _assert_connect_read_pair(recorder.timeout)


def test_launch_workflow_sets_a_timeout(monkeypatch):
    recorder = _Recorder(WORKFLOW_LAUNCH)
    monkeypatch.setattr(nextflow.requests, "post", recorder)

    assert nextflow.launch_workflow(_StubPipelineParams()) == "wf-123"

    _assert_connect_read_pair(recorder.timeout)


def test_get_workflow_status_sets_a_timeout(monkeypatch):
    recorder = _Recorder(WORKFLOW_STATUS)
    monkeypatch.setattr(nextflow.requests, "get", recorder)

    asyncio.run(nextflow.get_workflow_status("wf-123"))

    _assert_connect_read_pair(recorder.timeout)


# --- 2. the async helpers keep blocking work off the event loop ------------


def _thread_the_request_ran_on(recorder, coroutine_factory):
    """Run the coroutine and report the event-loop thread alongside the thread
    the (fake) request actually executed on."""
    loop_thread_id = {}

    async def run():
        loop_thread_id["value"] = threading.get_ident()
        return await coroutine_factory()

    asyncio.run(run())
    return loop_thread_id["value"], recorder.calls[0]["thread_id"]


def test_get_genome_metadata_does_not_block_the_event_loop(monkeypatch):
    recorder = _Recorder(GENOME_ATTRS)
    monkeypatch.setattr(web_metadata.requests, "get", recorder)

    loop_thread, request_thread = _thread_the_request_ran_on(
        recorder, lambda: web_metadata.get_genome_metadata("genome-1")
    )

    assert request_thread != loop_thread, (
        "requests.get ran on the event-loop thread; a slow metadata API would "
        "stall every other in-flight request"
    )


def test_get_workflow_status_does_not_block_the_event_loop(monkeypatch):
    recorder = _Recorder(WORKFLOW_STATUS)
    monkeypatch.setattr(nextflow.requests, "get", recorder)

    loop_thread, request_thread = _thread_the_request_ran_on(
        recorder, lambda: nextflow.get_workflow_status("wf-123")
    )

    assert request_thread != loop_thread, (
        "requests.get ran on the event-loop thread; the 15s status poll would "
        "stall every other in-flight request"
    )


# --- 3. errors still surface with the right function named ------------------


def test_status_errors_name_their_own_function(monkeypatch):
    """The status helper's error messages used to say launch_workflow()."""

    def _boom(*_args, **_kwargs):
        raise nextflow.requests.ConnectionError("refused")

    monkeypatch.setattr(nextflow.requests, "get", _boom)

    with pytest.raises(nextflow.requests.ConnectionError) as excinfo:
        asyncio.run(nextflow.get_workflow_status("wf-123"))

    assert "get_workflow_status()" in excinfo.value.args[0]
