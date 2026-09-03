#!/usr/bin/env python3
"""Independent checks for FIN ST537--ST551."""
import hashlib,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
class ST537ST551Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST537_ST551_Results.json").read_text())
 def test_programs_hashes(self):
  self.assertEqual(set(self.r),{f"ST{k}" for k in range(537,552)})
  for k in range(537,552):
   x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((ROOT/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
 def test_blockwise_exact_ranks(self):
  rows=self.r["ST537"]["prime_runs"]["1000003"]
  self.assertEqual([(x["rank"],x["shape"][1]) for x in rows],[(126,126),(672,672),(528,528),(702,702)])
  self.assertTrue(self.r["ST537"]["cross_prime_rank_agreement"])
  self.assertTrue(self.r["ST539"]["all_families_individually_Q_independent"])
 def test_cumulative_dependency_localization(self):
  rows=self.r["ST538"]["prime_runs"]["1000003"]
  self.assertEqual([x["rank"] for x in rows],[126,798,1288,1791])
  self.assertEqual([x["nullity"] for x in rows],[0,0,38,237])
 def test_transition_subfamilies(self):
  self.assertEqual(self.r["ST540"]["canonical_subsets"],2047);self.assertFalse(self.r["ST540"]["candidate_beaten"])
  self.assertEqual(self.r["ST541"]["starts"],240);self.assertFalse(self.r["ST541"]["candidate_beaten"])
 def test_noise_and_nonunique_flow(self):
  self.assertLess(self.r["ST542"]["maximum_frequency_deviation_from_1_12"],.001)
  self.assertFalse(self.r["ST543"]["strict_selector"]);self.assertIn("solution-selection",self.r["ST543"]["hidden_resource"])
 def test_finite_shot_and_CP_conjugacy(self):
  self.assertGreaterEqual(self.r["ST544"]["sufficient_sample_count"],1);self.assertFalse(self.r["ST544"]["event_level_physical_model"])
  self.assertLess(self.r["ST545"]["record_residual"],1e-12);self.assertFalse(self.r["ST545"]["physical_process_tensor"])
 def test_calibration_and_IR(self):
  self.assertEqual(self.r["ST546"]["residual_torsor_dimension"],1);self.assertFalse(self.r["ST546"]["absolute_seconds"])
  self.assertGreater(self.r["ST547"]["certified_stop"],.25);self.assertGreater(self.r["ST547"]["minimum_margin"],0)
  self.assertEqual(self.r["ST548"]["maximum_positive_clusters"],1);self.assertFalse(self.r["ST548"]["global_complement_theorem"])
 def test_gates_and_package_roles(self):
  self.assertFalse(self.r["ST549"]["new_strict_gain_source"]);self.assertFalse(self.r["ST549"]["new_nonpremise_selector"])
  self.assertFalse(self.r["ST550"]["strict_derivation"]);self.assertEqual(set(self.r["ST550"]["removal_matrix"]),{"remove_CA","remove_SA","remove_OA"})
  self.assertEqual(self.r["ST551"]["independent_laboratory_record"],"absent")
if __name__=="__main__":unittest.main(verbosity=2)
