"""Recreate a previously curated count carrier without running annotation."""
from __future__ import annotations

import importlib
import copy
import hashlib
import json
import shutil
import tarfile
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from http.client import IncompleteRead
from pathlib import Path

from scripts.download_pbmc_asset import ROOT, LIMIT, download, footprint


MODULES = {
    "pinned_h5_assignments": "pinned_h5_assignments",
    "pinned_cell_assignments": "pinned_cell_assignments",
    "pinned_author_counts": "pinned_author_counts",
    "cmo_assignments": "cmo_assignments",
    "atlas_geo_assignments": "atlas_geo_assignments",
    "paired_pool_metadata": "paired_pool_metadata",
    "pooled_pbmc_181897": "pooled_pbmc_181897",
    "gzip_h5ad": "gzip_h5ad",
    "paired_mtx_tables": "paired_mtx_tables",
    "rds_seurat": "rds_seurat",
    "ivanova_original": "ivanova_original",
    "native_mtx_bundle": "native_mtx_bundle",
    "ircolitis_original": "ircolitis_original",
    "rag_pbmc_original": "rag_pbmc_original",
    "lymphoma_original": "lymphoma_original",
    "cillo_original": "cillo_original",
    "reviewed_native_bundle": "reviewed_native_bundle",
    "reviewed_public_native_bundle": "public_native_bundle",
    "memory_b_original": "memory_b_original",
    "sorted_memory_original": "sorted_memory_original",
    "monocyte_plate_original": "monocyte_plate_original",
    "hsct_pbmc_original": "hsct_pbmc_original",
    "diabetes_initial": "diabetes_initial",
    "lyme_vireo": "lyme_vireo",
    "critical_illness_raw": "critical_illness_raw",
    "partitioned_h5ad_counts": "partitioned_h5ad_counts",
    "count_layer_h5ad": "count_layer_h5ad",
    "typed_h5ad": "typed_h5ad",
    "rdata_counts": "rdata_counts",
    "hca_eqtl45": "hca_eqtl45",
    "influenza_hashing": "influenza_hashing",
    "cll_hashing": "cll_hashing",
    "bcell_hashing": "bcell_hashing",
    "molecule_info": "molecule_info",
    "reviewed_panel_h5": "reviewed_panel_h5",
    "prefix_mtx": "prefix_mtx",
    "benchmark_archive_RNA": "benchmark_archive_RNA",
    "parse_geo": "parse_geo",
    "reviewed_series_mtx": "reviewed_series_mtx",
    "geo_mtx_table": "geo_mtx_table",
}


def download_ranges(url: str, destination: Path, n_bytes: int, expected_sha256: str) -> None:
    """Fetch a BioStudies file in bounded ranges when full-body GET stalls."""
    if destination.exists():
        with destination.open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() == expected_sha256:
                return
        raise ValueError("Existing asset differs from original curated checksum")
    part = destination.with_suffix(destination.suffix + ".part")
    offset = part.stat().st_size if part.exists() else 0
    if offset > n_bytes:
        raise ValueError('Partial range download exceeds pinned source size')
    remaining = n_bytes - offset
    if footprint() + remaining + 2_000_000_000 >= LIMIT or shutil.disk_usage(ROOT).free < remaining + 2_000_000_000:
        raise RuntimeError("Range download exceeds disk budget")
    destination.parent.mkdir(parents=True, exist_ok=True)
    chunk = 8 * 1024 * 1024
    windows = [(start, min(start + chunk, n_bytes) - 1) for start in range(offset, n_bytes, chunk)]

    def fetch(window):
        start, end = window
        request = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
        with urllib.request.urlopen(request, timeout=120) as response:
            if response.status != 206 or response.headers.get("Content-Range", "").split("/")[-1] != str(n_bytes):
                raise ValueError("Server did not honor the requested byte range")
            data = response.read()
        if len(data) != end - start + 1:
            raise ValueError("Incomplete ranged response")
        return data

    digest = hashlib.sha256()
    if offset:
        # Verify the existing prefix independently before appending. The
        # final pinned SHA checks every byte, including the saved partial.
        with part.open('rb') as saved:
            head=saved.read(min(offset,chunk))
        if fetch((0,len(head)-1))!=head:
            raise ValueError('Partial source prefix differs from pinned URL')
        with part.open('rb') as saved:
            for block in iter(lambda:saved.read(chunk),b''):
                digest.update(block)
    with ThreadPoolExecutor(max_workers=4) as pool, part.open("ab") as handle:
        for data in pool.map(fetch, windows):
            handle.write(data)
            digest.update(data)
    if digest.hexdigest() != expected_sha256:
        raise ValueError("Ranged source bytes differ from original curation")
    part.replace(destination)
    manifest = destination.with_suffix(destination.suffix + ".manifest.json")
    manifest.write_text(json.dumps({"url": url, "retrieved_at": datetime.now(timezone.utc).isoformat(),
                                    "bytes": n_bytes, "sha256": digest.hexdigest(),
                                    "response_headers": {"retrieval": "verified HTTP byte ranges"}}, indent=2),
                        encoding="utf-8")


def current_ebi_asset_urls(cfg: dict) -> dict:
    """Resolve legacy BioStudies file URLs while preserving pinned source bytes."""
    if not str(cfg.get("study_id", "")).startswith("E-MTAB-") or not cfg.get("geo_samples"):
        return cfg
    accession = cfg["study_id"]
    info_url = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{accession}/info"
    with urllib.request.urlopen(info_url, timeout=90) as response:
        info = json.load(response)
    replay_cfg = copy.deepcopy(cfg)
    asset_dir = ROOT / "data/downloads" / (accession + "_assets")
    saved_path = ROOT / "results/star_maps" / cfg["dataset_id"] / "source_manifest.json"
    saved_assets = json.loads(saved_path.read_text(encoding="utf-8")).get("assets", [])
    saved_by_name = {item["url"].rsplit("/", 1)[-1]: item for item in saved_assets}
    for sample in replay_cfg["geo_samples"]:
        for asset in sample["assets"]:
            existing_path = asset_dir / asset["filename"]
            if existing_path.exists():
                existing_manifest = existing_path.with_suffix(existing_path.suffix + ".manifest.json")
                if existing_manifest.exists():
                    asset["url"] = json.loads(existing_manifest.read_text(encoding="utf-8"))["url"]
                continue
            if asset["url"].startswith("https://www.ebi.ac.uk/biostudies/files/"):
                asset["url"] = info["httpLink"] + "/Files/" + asset["filename"]
                saved_item = saved_by_name.get(asset["filename"])
                if saved_item and saved_item["bytes"] >= 20_000_000:
                    download_ranges(asset["url"], existing_path,
                                    int(saved_item["bytes"]), saved_item["sha256"])
    return replay_cfg


def resume_pinned_count_layer_asset(cfg: dict) -> None:
    """Resume an interrupted large audited H5AD without changing its bytes."""
    asset = cfg["source_asset"]
    destination = ROOT / "data/downloads" / asset["filename"]
    part = destination.with_suffix(destination.suffix + ".part")
    if destination.exists() or not part.exists() or part.stat().st_size == 0:
        return
    total = int(asset["expected_bytes"])
    size = part.stat().st_size
    if size >= total:
        raise ValueError("Partial audited source has an unexpected length")
    saved = json.loads((ROOT / "results/star_maps" / cfg["dataset_id"] /
                        "source_manifest.json").read_text(encoding="utf-8"))["source"]
    if saved["url"] != asset["url"] or saved["bytes"] != total:
        raise ValueError("Audited source identity differs from current configuration")
    etag = saved.get("response_headers", {}).get("ETag")
    headers = {"Range": f"bytes={size}-{total - 1}"}
    if etag:
        headers["If-Range"] = etag
    if footprint() + (total - size) + 2_000_000_000 >= LIMIT or shutil.disk_usage(ROOT).free < total - size + 2_000_000_000:
        raise RuntimeError("Audited source resume exceeds disk budget")
    chunk = 8 * 1024 * 1024
    windows = [(start, min(start + chunk, total) - 1) for start in range(size, total, chunk)]

    def fetch(window):
        start, end = window
        for attempt in range(4):
            try:
                request = urllib.request.Request(asset["url"], headers={**headers,
                                                                         "Range": f"bytes={start}-{end}"})
                with urllib.request.urlopen(request, timeout=180) as response:
                    if response.status != 206 or response.headers.get("Content-Range") != f"bytes {start}-{end}/{total}":
                        raise ValueError("Server did not honor the pinned byte range")
                    if etag and response.headers.get("ETag") not in (None, etag):
                        raise ValueError("Audited source ETag changed")
                    data = response.read()
                if len(data) != end - start + 1:
                    raise IncompleteRead(data, end - start + 1 - len(data))
                return data
            except (OSError, IncompleteRead):
                if attempt == 3:
                    raise
                time.sleep(2 * (attempt + 1))

    with ThreadPoolExecutor(max_workers=8) as pool, part.open("ab") as handle:
        for index, data in enumerate(pool.map(fetch, windows), 1):
            handle.write(data)
            if index % 16 == 0:
                print(f"RESUME {cfg['dataset_id']} bytes={handle.tell()}/{total}", flush=True)
    if part.stat().st_size != total:
        raise ValueError("Resumed audited source is incomplete")
    with part.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    if digest != saved["sha256"]:
        raise ValueError("Resumed source differs from original audited checksum")
    part.replace(destination)
    manifest = destination.with_suffix(destination.suffix + ".manifest.json")
    manifest.write_text(json.dumps({"url": asset["url"],
                                    "retrieved_at": datetime.now(timezone.utc).isoformat(),
                                    "bytes": total, "sha256": digest,
                                    "response_headers": {"retrieval": "verified resume range",
                                                         "ETag": etag}}, indent=2) + "\n",
                        encoding="utf-8")


def prepare(cfg: dict):
    """Return source H5AD, original source manifest, and cleanup paths.

    It follows the original curated input-format branch and verifies preserved
    remote assets through the existing downloader and preparer manifests.
    """
    dataset_id = cfg["dataset_id"]
    source = ROOT / "data/downloads" / (dataset_id + ".h5ad")
    fmt = cfg.get("input_format", "h5ad")
    cleanup = [source]
    if fmt == "h5ad_tar":
        archive = ROOT / "data/downloads" / (dataset_id + ".h5ad.tar.gz")
        download(cfg["url"], archive, cfg["expected_bytes"])
        manifest = json.loads(archive.with_suffix(archive.suffix + ".manifest.json").read_text())
        with tarfile.open(archive) as handle:
            members = [m for m in handle if m.isfile() and m.name.endswith(".h5ad")]
            if len(members) != 1:
                raise ValueError("Expected exactly one H5AD member")
            member = members[0]
            if footprint() + member.size >= LIMIT or shutil.disk_usage(ROOT).free < member.size:
                raise RuntimeError("Archive expansion exceeds budget")
            with source.open("wb") as target:
                shutil.copyfileobj(handle.extractfile(member), target)
            manifest["selected_h5ad_member"] = member.name
        cleanup.append(archive)
    elif fmt == "geo_10x_tar":
        from src.geo_inputs import convert_geo_tar
        archive = ROOT / "data/downloads" / (dataset_id + "_RAW.tar")
        download(cfg["url"], archive, cfg["expected_bytes"])
        manifest = json.loads(archive.with_suffix(".tar.manifest.json").read_text())
        convert_geo_tar(archive, source, cfg)
        cleanup.append(archive)
    elif fmt == "geo_10x_assets":
        from src.geo_assets import prepare_assets
        from src.geo_inputs import convert_geo_tar
        replay_cfg = current_ebi_asset_urls(cfg)
        archive, manifest, assets = prepare_assets(replay_cfg)
        saved_path = ROOT / "results/star_maps" / dataset_id / "source_manifest.json"
        saved = json.loads(saved_path.read_text(encoding="utf-8"))
        if saved.get("kind") == "selected_source_assets":
            saved_checksums = sorted(item["sha256"] for item in saved["assets"])
            current_checksums = sorted(item["sha256"] for item in manifest["assets"])
            if saved_checksums != current_checksums:
                raise ValueError("Current source assets differ from originally curated bytes")
        convert_geo_tar(archive, source, cfg)
        cleanup.extend([archive, *assets])
    elif fmt in MODULES:
        if fmt == "count_layer_h5ad":
            resume_pinned_count_layer_asset(cfg)
        module = importlib.import_module("src." + MODULES[fmt])
        manifest, assets = module.prepare(current_ebi_asset_urls(cfg), source)
        cleanup.extend(assets)
    elif fmt in (None, "h5ad"):
        download(cfg["url"], source, cfg["expected_bytes"])
        manifest = json.loads(source.with_suffix(".h5ad.manifest.json").read_text())
    else:
        raise ValueError("Unsupported source format: " + str(fmt))
    if not source.exists():
        raise ValueError("Preparer did not create a count carrier")
    return source, manifest, [Path(item) for item in cleanup]
