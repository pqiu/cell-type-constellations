"""Validate source count units without treating normalized expression as raw counts."""
import numpy as np

def validate_count_values(values,cfg):
 semantics=cfg.get('count_semantics','integer_counts')
 if semantics not in ('integer_counts','salmon_alevin_expected_counts','author_corrected_fractional_counts','author_featurecounts_fractional_assignments'):
  raise ValueError('Unsupported source count semantics')
 if semantics=='salmon_alevin_expected_counts':
  evidence=cfg.get('count_semantics_evidence')
  if not isinstance(evidence,dict) or not all(evidence.get(k) for k in ('source_pipeline_url','source_pipeline_sha256','original_product','selected_layer','interpretation')):
   raise ValueError('Expected counts require explicit original product, layer and pipeline evidence')
 if semantics=='author_corrected_fractional_counts':
  evidence=cfg.get('count_semantics_evidence')
  required=('source_url','source_asset_sha256','original_product','selected_layer','interpretation')
  if not isinstance(evidence,dict) or not all(evidence.get(k) for k in required):
   raise ValueError('Author-corrected fractional counts require explicit source asset, product, layer and interpretation evidence')
 if semantics=='author_featurecounts_fractional_assignments':
  evidence=cfg.get('count_semantics_evidence')
  required=('source_url','source_asset_sha256','original_product','selected_layer','interpretation')
  if not isinstance(evidence,dict) or not all(evidence.get(k) for k in required):
   raise ValueError('Fractional featureCounts assignments require explicit source asset, product, layer and interpretation evidence')
 if not np.isfinite(values).all() or np.any(values<0):raise ValueError('Count values must be finite and nonnegative')
 if semantics=='integer_counts' and not np.allclose(values,np.round(values),atol=1e-6,rtol=0):
  raise ValueError('Noninteger values: cannot assume integer counts')
 # Values are validated in place, never rounded or transformed here.
