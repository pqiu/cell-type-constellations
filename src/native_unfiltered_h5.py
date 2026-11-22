"""Audit original human RNA droplet H5 without treating empty barcodes as cells."""
import re
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from src.lyme_vireo import audit_h5


def audit(path, seen, allowed_non_rna_types=()):
    with h5py.File(path, 'r') as f:
        m = f['matrix']; ft = m['features']
        ids = ft['id'].asstr()[()]; symbols = ft['name'].asstr()[()]
        types = ft['feature_type'].asstr()[()]; bars = m['barcodes'].asstr()[()]
        rna = types == 'Gene Expression'
        other_types = set(types[~rna])
        if (not rna.any() or not other_types.issubset(set(allowed_non_rna_types))
                or not all(re.fullmatch(r'ENSG\d+(?:\.\d+)?', x) for x in ids[rna])):
            raise ValueError('Only reviewed native human RNA is eligible')
        if sum(x.startswith('MT-') for x in symbols[rna]) != 13:
            raise ValueError('Native RNA mitochondrial gene coverage changed')
        x = sparse.csc_matrix((m['data'][()], m['indices'][()], m['indptr'][()]), shape=tuple(m['shape'][()])).T.tocsr()
        x = x[:, rna].tocsr()
        totals = np.asarray(x.sum(axis=1, dtype=np.int64)).ravel()
        selected = pd.DataFrame(index=bars[totals >= 500])
        source_cells = len(bars); positive = int((totals > 0).sum())
        del x
    # Complete native count validation includes subthreshold/empty droplets.
    # Early selection only limits duplicate fingerprints and AnnData size.
    review = audit_h5(path, selected, seen)
    review.update(source_barcode_rows=source_cells, RNA_positive_barcode_rows=positive,
        duplicate_fingerprint_scope='All barcodes meeting unchanged minRNAcounts500; complete native matrix integer/coordinate checks precede RNA QC',
        early_RNA_QC_eligible_barcodes=len(selected), mitochondrial_features=13,
        excluded_non_RNA_features=int((~rna).sum()), reviewed_non_RNA_types=sorted(other_types))
    return review
