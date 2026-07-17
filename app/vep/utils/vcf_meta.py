"""VCF metadata (variant / header counts) with a local file cache.

Extracted from vcf_results. Used by the bcftools fallback path in
get_results_from_path when no page-index sidecar is present.
"""

import subprocess
from pydantic import FilePath
from vep.models import vcf_results_model as model

META_FILE = "results_meta.json"


# ---------------------------------------------------------------------------
# DEV/LOCAL-ONLY: stale metadata-cache guard
#
# In the full pipeline, every run writes its output VCF (and this metadata
# cache) into its own directory, so the cache can never go stale. This guard
# only matters when a single local output file is parsed repeatedly and then
# regenerated (e.g. re-running the pipeline against the fixed LOCAL_RESULTS_VCF
# path during development): without it the old header/record counts would be
# reused and mis-slice the new file. Safe to remove once outputs always land in
# per-run directories.
# ---------------------------------------------------------------------------
def _is_meta_cache_stale(meta_path: FilePath, vcf_path: FilePath) -> bool:
    """True if the metadata cache exists but predates its VCF."""
    return (
        meta_path.exists()
        and meta_path.stat().st_mtime < vcf_path.stat().st_mtime
    )


def _get_vcf_meta(vcf_path: FilePath) -> model.VcfMetadata:
    """Helper method to manage metainfo for a VCF file"""

    meta_path = vcf_path.with_name(META_FILE)
    if not meta_path.exists() or _is_meta_cache_stale(meta_path, vcf_path):
        variant_count_str = subprocess.check_output(
            f"bcftools stats {vcf_path} | grep 'number of records:'",
            shell=True, text=True
        )
        header_count_str = subprocess.check_output(
            f"bcftools view -h {vcf_path} | wc -l",
            shell=True, text=True
        )
        try:
            vcf_info = model.VcfMetadata(
                variant_count=int(variant_count_str.split(":")[-1]),
                header_count=int(header_count_str)
            )
        except ValueError as e:
            e.args = (
                f"_get_vcf_meta: unexpected bcftools output: variant_count: {variant_count_str} | header_count: {header_count_str}",
                *e.args,
            )
            raise

        with open(meta_path, "w") as meta_file:
            meta_file.write(vcf_info.model_dump_json())
    else:
        with open(meta_path, "r") as meta_file:
            vcf_info = model.VcfMetadata.model_validate_json(meta_file.read())
    return vcf_info
