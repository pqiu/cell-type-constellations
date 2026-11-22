"""Exact raw-count fingerprints for candidate duplicate detection, never donor identification."""
import hashlib
import numpy as np
from scipy import sparse

def count_fingerprints(x,minimum_detected=100,minimum_sum=500):
 x=sparse.csr_matrix(x,dtype=np.float64,copy=True);x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
 if not np.isfinite(x.data).all() or np.any(x.data<0):raise ValueError('Invalid source count values')
 out=[]
 for i in range(x.shape[0]):
  start,stop=x.indptr[i:i+2];values=x.data[start:stop];indices=x.indices[start:stop]
  if len(values)<minimum_detected or values.sum()<minimum_sum:out.append(None);continue
  h=hashlib.sha256();h.update(np.asarray(indices,dtype='<i4').tobytes());h.update(np.asarray(values,dtype='<f8').tobytes());out.append(h.digest())
 return out
