"""Independently validate complete labeled dense count rows with bounded workers."""
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from io import BytesIO
import hashlib
import multiprocessing
import numpy as np


def _row_counts(task):
    index, expected_label, cells, line = task
    label, separator, values = line.partition(b'\t')
    if not separator or label.decode('utf-8') != expected_label:
        raise ValueError('Original feature-row identity differs')
    # loadtxt supplies a second parser, distinct from the staging fromstring path.
    counts = np.loadtxt(BytesIO(values), delimiter='\t', dtype=np.float64, ndmin=1)
    if counts.shape != (cells,) or not np.isfinite(counts).all():
        raise ValueError('Wrong count-row width or nonfinite values')
    if (counts < 0).any() or not np.equal(counts,np.floor(counts)).all() or counts.max(initial=0) > np.iinfo(np.int32).max:
        raise ValueError('Counts must be finite nonnegative int32-representable integers')
    return counts.astype('<i8')


def _signature(index,label,counts):
    return dict(original_feature_row=index,original_feature_label=label,
                nonzero_cells=int(np.count_nonzero(counts)),total_raw_counts=int(counts.sum()),
                maximum_raw_count=int(counts.max(initial=0)),
                complete_original_row_sha256=hashlib.sha256(counts.tobytes()).hexdigest())


def row_signature(task):
    counts=_row_counts(task)
    return _signature(task[0],task[1],counts)


def row_signature_and_sparse(task):
    counts=_row_counts(task);indices=np.flatnonzero(counts).astype(np.int32)
    return _signature(task[0],task[1],counts),counts[indices].astype(np.int32),indices


def verified_sparse_rows(stream,labels,cells,expected,workers=8):
    """Recreate sparse rows from a second parser, checking every original-row hash."""
    if not 1<=workers<=8:raise ValueError('Use one to eight bounded workers')
    pending=deque();seen=verified=0
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        for index,line in enumerate(stream):
            if index>=len(labels):raise ValueError('Extra original feature row')
            pending.append(pool.submit(row_signature_and_sparse,(index,labels[index],cells,line)));seen=index+1
            if len(pending)>=3*workers:
                signature,data,indices=pending.popleft().result();assert signature==expected[verified];yield verified,data,indices;verified+=1
        while pending:
            signature,data,indices=pending.popleft().result();assert signature==expected[verified];yield verified,data,indices;verified+=1
    if seen!=len(labels) or verified!=len(labels):raise ValueError('Missing original feature rows')


def validate_rows(stream, labels, cells, expected, workers=8, progress=None):
    """Check all original rows in order; queue memory is bounded to 3 rows/worker."""
    if not 1 <= workers <= 8:raise ValueError('Use one to eight bounded workers')
    pending=deque();next_index=0;verified=0
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        for index,line in enumerate(stream):
            if index >= len(labels):raise ValueError('Extra original feature row')
            pending.append(pool.submit(row_signature,(index,labels[index],cells,line)))
            next_index=index+1
            if len(pending)>=3*workers:
                result=pending.popleft().result();assert result==expected[verified], 'Original values differ from sparse-block proof';verified+=1
                if progress is not None and verified%512==0:progress(verified)
        while pending:
            result=pending.popleft().result();assert result==expected[verified], 'Original values differ from sparse-block proof';verified+=1
            if progress is not None and verified%512==0:progress(verified)
    if next_index!=len(labels) or verified!=len(labels):raise ValueError('Missing original feature rows')
    return verified
