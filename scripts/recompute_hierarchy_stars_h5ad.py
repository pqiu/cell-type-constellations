"""Recompute exact parent stars using cached labels and a verified source H5AD.

This source adapter handles directly downloaded and verified prepared H5ADs.
The full source gene axis is used for each cell's library size. The unfiltered
mode materializes all present parent profiles on the fixed 2,000-gene axis and
saves direct-label sums needed to reconstruct their native-gene profiles.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from inspect_pooled_star_axes import ROOT
from pooled_star_umap import axis
sys.path.insert(0, str(ROOT))
from src.hierarchy_stars import log_normalized_sum, pooled_parent_means


SOURCE = ROOT / "results/pooled_star_translation"
DESTINATION = ROOT / "results/hierarchy_star_maps"


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def parent_design(extra_five: bool = False, all_ancestors: bool = False):
    if extra_five and all_ancestors:
        raise ValueError('Select one parent inventory')
    selected = (SOURCE / "filtered_min20/retained_labels.txt").read_text(encoding="utf-8").splitlines()
    mapping = pd.read_csv(ROOT / "references/vocabulary/label_mapping.tsv", sep="\t")
    membership = pd.read_csv(ROOT / "references/vocabulary/aggregation_membership.tsv", sep="\t")
    source_labels = sorted(mapping.original_label)
    parent_ids = dict(zip(mapping.original_label, mapping.canonical_id))
    if all_ancestors:
        parent_labels = sorted(membership.aggregate_id.unique())
        parent_ids.update({identity: identity for identity in parent_labels})
    elif extra_five:
        coverage = pd.read_csv(ROOT / 'results/hierarchy_global_reanalysis/hierarchy_label_coverage.csv')
        old_ids = set(mapping.loc[mapping.original_label.isin(selected), 'canonical_id'])
        new = coverage.loc[coverage.retained & ~coverage.canonical_id.isin(old_ids)]
        if len(new) != 5:
            raise AssertionError('New parent inventory changed')
        parent_labels = new.canonical_name.tolist()
        parent_ids.update(zip(new.canonical_name, new.canonical_id))
    else:
        parent_labels = selected
    members = membership.groupby("aggregate_id").original_label.apply(set)
    incidence = np.zeros((len(source_labels), len(parent_labels)), dtype=np.int8)
    source_index = {label: i for i, label in enumerate(source_labels)}
    for j, label in enumerate(parent_labels):
        for child in members[parent_ids[label]]:
            incidence[source_index[child], j] = 1
    expected = 379 if all_ancestors else (5 if extra_five else 79)
    if len(parent_labels) != expected or len(source_labels) != 203:
        raise AssertionError("Frozen label inventory changed")
    return source_labels, parent_labels, parent_ids, sparse.csr_matrix(incidence)


def process(dataset_id: str, source: Path, verified_manifest: dict | None = None,
            extra_five: bool = False, all_ancestors: bool = False) -> dict:
    original = ROOT / "results/star_maps" / dataset_id
    run = json.loads((original / "run.json").read_text(encoding="utf-8"))
    cfg = run["configuration"]
    saved_manifest = json.loads((original / "source_manifest.json").read_text(encoding="utf-8"))
    source_digest = sha256(source)
    if verified_manifest is None:
        if cfg.get("input_format", "h5ad") not in ("h5ad", None):
            raise ValueError("A prepared source needs its verified manifest")
        if saved_manifest.get("sha256") != source_digest:
            raise ValueError("Redownloaded count file differs from curated source")
    elif saved_manifest.get("sha256") and saved_manifest["sha256"] != verified_manifest.get("sha256"):
        raise ValueError("Redownloaded source asset differs from curated source")
    out_root = (ROOT/'results/hierarchy_all_star_maps' if all_ancestors else
                ROOT/'results/hierarchy_star_maps_extra5' if extra_five else DESTINATION)
    out = out_root / dataset_id
    out.mkdir(parents=True, exist_ok=True)
    if (out / "report.json").exists():
        prior = json.loads((out / "report.json").read_text(encoding="utf-8"))
        if prior.get('parent_labels') != (379 if all_ancestors else 5 if extra_five else 79):
            raise ValueError('Existing output uses a different parent inventory')
        if prior.get("source_sha256") != source_digest:
            raise ValueError("Existing hierarchy output has different source provenance")
        if not all((out / name).exists() and sha256(out / name) == digest
                   for name, digest in prior.get("output_sha256", {}).items()):
            raise ValueError("Existing hierarchy output failed checksum verification")
        expected_native_file = 'under20_direct_sums_native.h5' if all_ancestors else 'expression_native.h5'
        layout_ok = (prior.get('storage_layout')=='present_only_with_expression_row'
                     if all_ancestors else True)
        if layout_ok and expected_native_file in prior.get("output_sha256", {}) and all(
            (out / name).exists() and sha256(out / name) == digest
            for name, digest in prior.get("output_sha256", {}).items()
        ):
            return prior
    source_labels, parent_labels, parent_ids, incidence = parent_design(
        extra_five=extra_five, all_ancestors=all_ancestors)
    label_index = {label: i for i, label in enumerate(source_labels)}
    genes = (SOURCE / "selected_genes.txt").read_text(encoding="utf-8").splitlines()
    gene_positions = axis(original / "genes.csv")
    columns = np.asarray([gene_positions[gene] for gene in genes], dtype=np.int64)
    metadata_path = original / "cell_metadata.csv.gz"
    index_name = pd.read_csv(metadata_path, nrows=0).columns[0]
    metadata_columns = [index_name, "retained", "sample_id", "scimilarity_label"]
    if cfg.get("use_source_gene_availability", False):
        metadata_columns.append("source_library_index")
    cached = pd.read_csv(metadata_path, index_col=0,
                         usecols=metadata_columns, low_memory=False)
    if not cached.retained.dtype == bool:
        raise ValueError("Cached retention flags are not Boolean")
    stars = pd.read_csv(original / "stars.csv", usecols=["sample_id", "label", "n_cells", "status"])
    samples = pd.read_csv(original / "samples.csv", usecols=["sample_id", "donor_uid"])
    expected_parent_counts = None
    expected_supported_rows = None
    expected_present_rows = None
    parent_names = None
    if all_ancestors:
        direct_count_table = stars.pivot(index='sample_id', columns='label', values='n_cells')
        direct_count_table = direct_count_table.reindex(index=samples.sample_id, columns=source_labels)
        if direct_count_table.isna().any().any():
            raise ValueError('Original direct-label count inventory is incomplete')
        expected_parent_counts = np.asarray(
            direct_count_table.to_numpy(dtype=np.int64) @ incidence)
        expected_supported_rows = int((expected_parent_counts >= cfg['min_cells_per_star']).sum())
        expected_present_rows = int((expected_parent_counts > 0).sum())
        expected_under20_rows = int(stars.status.eq('insufficient_support').sum())
        membership = pd.read_csv(ROOT/'references/vocabulary/aggregation_membership.tsv',sep='\t')
        parent_names = (membership[['aggregate_id','aggregate_name']]
                        .drop_duplicates('aggregate_id').set_index('aggregate_id').aggregate_name.to_dict())
    a = ad.read_h5ad(source, backed="r")
    try:
        if len(cached) != a.n_obs or not np.array_equal(cached.index.astype(str), a.obs.index.astype(str)):
            raise ValueError("Cached cell annotations do not match source cell order")
        if cfg["raw_location"] == "raw.X" and a.raw is not None:
            count_matrix, var = a.raw.X, a.raw.var
        elif cfg["raw_location"] == "X" or (
            cfg["raw_location"] == "raw.X" and a.raw is None and
            (cfg["dataset_id"] in {"GSE309595_genes001", "GSE318371"} or
             bool(a.uns.get("count_carrier_X_original", False)))
        ):
            count_matrix, var = a.X, a.var
        else:
            raise ValueError("Source raw count location changed")
        genes_table = pd.read_csv(original / "genes.csv", usecols=["source_gene_id"],
                                  keep_default_na=False)
        if not np.array_equal(var.index.astype(str), genes_table.source_gene_id.astype(str)):
            raise ValueError("Count gene order does not match original export")
        availability = None
        if cfg.get("use_source_gene_availability", False):
            availability = np.asarray(a.uns["source_gene_availability"], dtype=bool)
            if not np.array_equal(np.asarray(a.uns["source_gene_ids"]).astype(str),
                                  var.index.astype(str)):
                raise ValueError("Source gene availability axis mismatch")
        if len(genes) != 2000:
            raise AssertionError("Expected the frozen 2,000-gene axis")
        native_genes = len(var)
        native_columns = np.arange(native_genes, dtype=np.int64)
        if all_ancestors:
            from scripts.download_pbmc_asset import LIMIT, footprint
            # A conservative uncompressed bound; actual gzip output is smaller.
            upper_bound = 4*(expected_present_rows*len(genes) +
                             expected_under20_rows*native_genes) + 1_000_000_000
            used=footprint()
            free=shutil.disk_usage(ROOT).free
            if used+upper_bound+1_000_000_000>=LIMIT or free<upper_bound+1_000_000_000:
                raise RuntimeError('Unfiltered output cannot fit inside the 1 TB project budget')
        by_sample = cached.groupby("sample_id", sort=False).indices
        direct_expression = original / "expression.h5"
        if not direct_expression.exists() and all_ancestors:
            direct_expression = out / "native_direct_mean_basis.h5"
            receipt = json.loads((out / "retired_direct_label_expression.json").read_text())
            if (receipt['original_sha256'] != run['output_sha256']['expression.h5'] or
                    sha256(direct_expression) != receipt['basis_sha256']):
                raise ValueError('Retired direct-label mean basis fails provenance checks')
        old = h5py.File(direct_expression, "r")
        out_h5 = h5py.File(out / "expression.h5", "w")
        native_filename = 'under20_direct_sums_native.h5' if all_ancestors else 'expression_native.h5'
        native_h5 = h5py.File(out / native_filename, "w")
        total_rows = expected_present_rows if all_ancestors else len(samples) * len(parent_labels)
        compression_options = dict(compression='gzip', compression_opts=1, shuffle=True) if all_ancestors else dict(compression='gzip')
        result = out_h5.create_dataset("norm10k_log1p_mean", shape=(total_rows, len(genes)),
                                       dtype="float32", chunks=(1, len(genes)),
                                       fillvalue=np.nan,
                                       **({'maxshape':(None,len(genes))} if all_ancestors else {}),
                                       **compression_options)
        if all_ancestors:
            basis_result = native_h5.create_dataset('normalized_log_sum',
                                                    shape=(expected_under20_rows,native_genes),
                                                    maxshape=(None,native_genes),
                                                    dtype='float32',chunks=(1,native_genes),
                                                    fillvalue=np.nan,**compression_options)
            native_result = None
        else:
            native_result = native_h5.create_dataset("norm10k_log1p_mean",
                                                     shape=(total_rows, native_genes),
                                                     dtype="float32", chunks=(1, native_genes),
                                                     fillvalue=np.nan, **compression_options)
            basis_result = None
        records = []
        basis_records = []
        basis_cursor = 0
        expression_cursor = 0
        under20_cells = 0
        direct_cells = 0
        checked_supported_means = 0
        try:
            star_groups = {sid: group for sid, group in stars.groupby("sample_id", sort=False)}
            for sample_no, sample in enumerate(samples.itertuples(index=False)):
                sid = sample.sample_id
                sample_rows = star_groups[sid]
                source_row_by_label = dict(zip(sample_rows.label, sample_rows.index))
                sample_rows = sample_rows.set_index("label")
                if set(sample_rows.index) != set(source_labels):
                    raise ValueError("Original star label inventory is incomplete")
                row_order = sample_rows.loc[source_labels]
                counts = row_order.n_cells.to_numpy(dtype=np.int64)
                direct_cells += int(counts.sum())
                positions = np.asarray(by_sample[sid], dtype=np.int64)
                positions = positions[cached.retained.to_numpy()[positions]]
                observed_labels = cached.scimilarity_label.to_numpy()[positions]
                observed_counts = pd.Series(observed_labels).value_counts()
                for label, count in zip(source_labels, counts):
                    if int(observed_counts.get(label, 0)) != count:
                        raise ValueError(f"Saved labels and star counts differ for {sid}, {label}")
                direct_sums = np.zeros((len(source_labels), native_genes), dtype=np.float64)
                # Supported direct means already summarize all their cells exactly.
                for label_i, label in enumerate(source_labels):
                    n = counts[label_i]
                    if n == 0:
                        continue
                    if row_order.loc[label, "status"] == "supported":
                        saved_mean = old["norm10k_log1p_mean"][source_row_by_label[label], :]
                        direct_sums[label_i] = np.asarray(saved_mean, dtype=np.float64) * n
                        if sample_no == 0 and checked_supported_means < 3 and np.isfinite(saved_mean).all():
                            subset = positions[observed_labels == label]
                            reconstructed = log_normalized_sum(count_matrix[subset], columns) / n
                            if np.max(np.abs(reconstructed - saved_mean[columns])) > 2e-6:
                                raise ValueError("Reopened counts do not reproduce saved direct star")
                            checked_supported_means += 1
                    else:
                        subset = positions[observed_labels == label]
                        if len(subset) != n or n >= cfg["min_cells_per_star"]:
                            raise AssertionError("Unsupported label support mismatch")
                        direct_sums[label_i] = log_normalized_sum(count_matrix[subset], native_columns)
                        if availability is not None:
                            libraries = cached.source_library_index.to_numpy(dtype=np.int32)[subset]
                            absent_genes = ~availability[libraries].all(axis=0)
                            direct_sums[label_i, absent_genes] = np.nan
                        if all_ancestors:
                            basis_result[basis_cursor] = direct_sums[label_i].astype(np.float32)
                            basis_records.append({'sample_id':sid,'label':label,'n_cells':int(n),
                                                  'basis_row':basis_cursor})
                            basis_cursor += 1
                        under20_cells += int(n)
                means, parent_counts = pooled_parent_means(direct_sums, counts, incidence,
                                                            min_cells=1 if all_ancestors else cfg["min_cells_per_star"])
                if not np.array_equal(parent_counts, np.asarray(counts @ incidence).ravel()):
                    raise AssertionError("Parent cell counts changed unexpectedly")
                if all_ancestors and not np.array_equal(parent_counts, expected_parent_counts[sample_no]):
                    raise AssertionError('Precomputed parent support changed')
                for j, label in enumerate(parent_labels):
                    n = int(parent_counts[j])
                    row_i = expression_cursor if all_ancestors else sample_no * len(parent_labels) + j
                    if (n > 0 if all_ancestors else n >= cfg["min_cells_per_star"]):
                        if not all_ancestors:
                            native_result[row_i] = means[j].astype(np.float32)
                        result[row_i] = means[j, columns].astype(np.float32)
                        if all_ancestors:
                            expression_cursor += 1
                    records.append({"sample_id": sid, "donor_uid": sample.donor_uid,
                                    "label": parent_names[label] if all_ancestors else label,
                                    "parent_id": parent_ids[label],
                                    "n_cells": n, "status": "supported" if n >= cfg["min_cells_per_star"] else
                                    ("absent" if n == 0 else "insufficient_support"),
                                    **({'expression_row': row_i if n > 0 else None}
                                       if all_ancestors else {})})
            if all_ancestors and expression_cursor != expected_present_rows:
                raise AssertionError('Stored expression-row count changed')
            if all_ancestors and basis_cursor != expected_under20_rows:
                raise AssertionError('Stored under-20 direct sum count changed')
            out_h5.flush()
            native_h5.flush()
        finally:
            old.close()
            out_h5.close()
            native_h5.close()
    finally:
        a.file.close()
    pd.DataFrame(records).to_csv(out / "stars.csv", index=False)
    if all_ancestors:
        pd.DataFrame(basis_records,columns=['sample_id','label','n_cells','basis_row']).to_csv(
            out/'under20_direct_rows.csv',index=False)
    report = {
        "dataset_id": dataset_id,
        "source_sha256": source_digest,
        "curated_source_manifest_sha256": saved_manifest.get("sha256"),
        "cached_annotation": str(original / "cell_metadata.csv.gz"),
        "annotation_model_rerun": False,
        "gene_count": len(genes),
        "native_gene_count": native_genes,
        "parent_labels": len(parent_labels),
        "candidate_parent_rows": len(samples)*len(parent_labels),
        "present_parent_rows": expression_cursor if all_ancestors else int(sum(x['status']!='absent' for x in records)),
        "supported_parent_rows": expected_supported_rows if all_ancestors else int(sum(x['status']=='supported' for x in records)),
        "low_support_parent_rows": (expression_cursor-expected_supported_rows) if all_ancestors else None,
        "storage_layout": 'present_only_with_expression_row' if all_ancestors else 'all_candidate_rows',
        "native_profile_storage": 'factorized_direct_sums' if all_ancestors else 'materialized_parent_means',
        "under20_direct_basis_rows": basis_cursor if all_ancestors else None,
        "original_direct_expression_sha256": run['output_sha256']['expression.h5'] if all_ancestors else None,
        "samples": len(samples),
        "direct_labeled_cells": direct_cells,
        "cells_recomputed_from_source_counts_because_no_direct_mean": under20_cells,
        "supported_direct_means_independently_checked": checked_supported_means,
        "parent_centroid": "Mean of each contributing cell's log1p(TP10K), with union of deduplicated ontology descendant labels",
        "parent_support_threshold_cells": int(cfg["min_cells_per_star"]),
        "output_sha256": {name: sha256(out / name) for name in
                          (("stars.csv", "expression.h5", "under20_direct_sums_native.h5",
                            "under20_direct_rows.csv") if all_ancestors else
                           ("stars.csv", "expression.h5", "expression_native.h5"))},
    }
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_id")
    parser.add_argument("--source", type=Path)
    args = parser.parse_args()
    source = args.source or ROOT / "data/downloads" / (args.dataset_id + ".h5ad")
    print(json.dumps(process(args.dataset_id, source), indent=2))


if __name__ == "__main__":
    main()
