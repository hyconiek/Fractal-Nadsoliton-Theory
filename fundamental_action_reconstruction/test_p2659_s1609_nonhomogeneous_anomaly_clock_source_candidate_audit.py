from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SCRIPT=ROOT/'p2659_s1609_nonhomogeneous_anomaly_clock_source_candidate_audit.py'
OUT=ROOT/'generated'/'p2659_s1609_nonhomogeneous_anomaly_clock_source_candidate_audit.json'
MD=ROOT/'generated'/'p2659_s1609_nonhomogeneous_anomaly_clock_source_candidate_audit.md'
class P2659NonhomogeneousAnomalyClockSourceCandidateAuditTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  subprocess.run([sys.executable,str(SCRIPT)],cwd=ROOT,check=True)
  cls.payload=json.loads(OUT.read_text())
 def test_content_first_and_upstream(self):
  audit=self.payload['semantic_rg_antiduplication_audit']; self.assertIn('content-first',audit['mode'])
  for key in ('nonhomogeneous_anomaly_content','scale_breaking_content','homogeneous_no_go_content','nonclosure_guard_content'):
   self.assertIn(key,audit['patterns']); self.assertGreater(audit['patterns'][key]['count'],0)
  up=self.payload['upstream_consistency']; self.assertTrue(up['p2657_unique_scale_not_selected']); self.assertTrue(up['p2658_homogeneous_no_go_verified']); self.assertTrue(up['p2658_anomaly_target_left_open'])
 def test_anomaly_witness_is_not_source(self):
  a=self.payload['anomaly_candidate_audit']; self.assertTrue(a['all_integer_phase_conditions_satisfied']); self.assertLess(a['max_phase_error'],1e-11)
  self.assertTrue(a['nonzero_lambda_breaks_pure_homogeneous_covariance']); self.assertTrue(a['declared_absolute_action_selectors_are_external']); self.assertFalse(a['uv_unit_selected_by_intrinsic_anomaly_now'])
 def test_source_matrix_blocks_false_passes(self):
  m={r['candidate']:r for r in self.payload['source_candidate_matrix']}
  self.assertFalse(m['affine_nonhomogeneous_trace_plus_declared_lambda']['passes_as_uv_unit_source_now']); self.assertTrue(m['affine_nonhomogeneous_trace_plus_declared_lambda']['uses_external_clock_or_scale_anchor'])
  self.assertFalse(m['integer_phase_with_nonhomogeneous_action']['passes_as_uv_unit_source_now'])
  self.assertFalse(m['declared_absolute_action_quantum_for_anomaly']['passes_as_uv_unit_source_now'])
  self.assertFalse(m['derived_boundary_anomaly_with_fixed_coefficient']['covered_by_audit'])
 def test_no_closure_and_docs(self):
  d=self.payload['closure_decision']; self.assertEqual(d['passing_uv_unit_source_candidates'],[]); self.assertFalse(d['uv_unit_selected_now']); self.assertFalse(d['beta_source_exported_now']); self.assertFalse(d['role_bearing_ltotal_now']); self.assertFalse(d['toe_closure_now'])
  self.assertTrue(all(v is False for v in self.payload['negative_export_flags'].values()))
  self.assertIn('P2659/S1609',MD.read_text()); self.assertIn('P2659/S1609',(ROOT/'STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md').read_text()); self.assertIn('P2659/S1609',(ROOT/'STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md').read_text())
if __name__=='__main__': unittest.main()
