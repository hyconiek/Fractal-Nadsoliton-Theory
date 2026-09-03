#!/usr/bin/env python3
"""Independent checks for FIN ST522--ST536."""
import hashlib,json,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
class ST522ST536Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST522_ST536_Results.json").read_text())
 def test_programs_hashes(self):
  self.assertEqual(set(self.r),{f"ST{k}" for k in range(522,537)})
  for k in range(522,537):
   x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((ROOT/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
 def test_conditioned_exact_bases(self):
  self.assertEqual(self.r["ST522"]["rank"],1791);self.assertEqual(self.r["ST522"]["free_count"],237)
  x=self.r["ST523"];self.assertTrue(x["all_ranks_1791"]);self.assertTrue(x["all_conditioned_pivots_admitted"]);self.assertTrue(x["all_unseen_seed_replays_zero"])
  p=ROOT/x["basis_file"];self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),x["basis_sha256"]);z=np.load(p);self.assertEqual(len(z.files),5)
 def test_conditioning_did_not_fake_Q_closure(self):
  self.assertEqual(self.r["ST524"]["reconstructed_count"],0);self.assertFalse(self.r["ST524"]["conditioning_improved_over_ST509"])
  self.assertFalse(self.r["ST526"]["rank_Q_equals_1791"]);self.assertEqual(self.r["ST526"]["accepted_rank_interval"],[1791,1892])
 def test_no_go_dichotomies(self):
  self.assertFalse(self.r["ST528"]["gain_from_zero_without_source"]);self.assertEqual(self.r["ST529"]["stochastic_invariant_selector"],"uniform 1/12");self.assertFalse(self.r["ST530"]["strict_scale_source"])
  self.assertEqual(set(self.r["ST527"]["hypothesis_removal_witnesses"]),{"equivariance","exact_uniform_initial_state","determinism_or_uniqueness","transitivity","tangency"})
 def test_three_channels_and_tester(self):
  self.assertTrue(all(v>0 for v in self.r["ST531"]["pairwise_distances"].values()));self.assertFalse(self.r["ST531"]["single_generator_selects_channel"])
  self.assertTrue(self.r["ST532"]["positive_after_payment"]);self.assertFalse(self.r["ST532"]["finite_shot_likelihood"])
 def test_ir_and_transition_status(self):
  self.assertGreater(self.r["ST533"]["certified_stop"],.24);self.assertGreater(self.r["ST533"]["minimum_margin"],0)
  x=self.r["ST534"];self.assertEqual(x["candidate_orbit_size"],12);self.assertGreater(x["candidate_minimum_paid_local_curvature"],5.6);self.assertFalse(x["candidate_globally_first"])
 def test_gates(self):
  x=self.r["ST535"];self.assertFalse(x["new_strict_gain_source"]);self.assertFalse(x["new_nonpremise_selector"]);self.assertFalse(x["new_scale_charged_source"]);self.assertFalse(x["characteristic_zero_degree8_rank_exact"])
  self.assertEqual(self.r["ST536"]["independent_laboratory_record"],"absent")
if __name__=="__main__":unittest.main(verbosity=2)
