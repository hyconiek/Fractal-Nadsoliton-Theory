import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_exclusion_exact_cover_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_exclusion_exact_cover_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_exclusion_exact_cover_certificate_report.md"


class StrictAlphaCyclicExclusionExactCoverCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_CYCLIC_EXCLUSION_EXACT_COVER_CERTIFICATE_PROBE__CONDITIONAL_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "exact-cover-nearest-plus-antipodal-constraints-select-A5-d5-orbit-conditionally-source-not-derived",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_constraint_system(self):
        constraints = self.payload["constraint_system"]
        self.assertEqual(len(constraints["variables"]), 12)
        self.assertEqual(constraints["cardinality"], "sum x_i = 5")
        self.assertEqual(len(constraints["nearest_forbidden_pairs"]), 12)
        self.assertEqual(len(constraints["antipodal_forbidden_pairs"]), 6)
        self.assertIn([0, 1], constraints["nearest_forbidden_pairs"])
        self.assertIn([11, 0], constraints["nearest_forbidden_pairs"])
        self.assertIn([0, 6], constraints["antipodal_forbidden_pairs"])

    def test_exact_cover_solution_audit(self):
        audit = self.payload["exact_cover_solution_audit"]
        self.assertEqual(audit["variable_count"], 12)
        self.assertEqual(audit["nearest_edge_constraint_count"], 12)
        self.assertEqual(audit["antipodal_edge_constraint_count"], 6)
        self.assertEqual(audit["forbidden_pair_count"], 18)
        self.assertEqual(audit["solution_count"], 12)
        self.assertEqual(audit["solution_histogram_rows"], [{"distance_histogram_d1_to_d6": [0, 3, 2, 1, 4, 0], "support_count": 12}])
        self.assertEqual(audit["solution_orbit_count"], 1)
        self.assertTrue(audit["selects_A5_d5_orbit"])
        orbit = audit["solution_orbits"][0]
        self.assertEqual(orbit["orbit_size"], 12)
        self.assertTrue(orbit["is_A5_d5_orbit"])
        self.assertFalse(orbit["is_A1_contiguous_orbit"])
        necklace = audit["solution_gap_necklaces"][0]
        self.assertEqual(necklace["canonical_gap_necklace"], [2, 2, 3, 2, 3])
        self.assertFalse(necklace["has_antipodal_contact"])
        self.assertTrue(necklace["is_A5_d5_necklace"])

    def test_d1_zero_gap_boundary(self):
        gap = self.payload["d1_zero_gap_classification"]
        self.assertEqual(gap["d1_zero_support_count"], 36)
        self.assertEqual(gap["d1_zero_gap_necklace_count"], 3)
        necklaces = {tuple(row["canonical_gap_necklace"]): row for row in gap["d1_zero_gap_necklaces"]}
        self.assertEqual(set(necklaces), {(2, 2, 2, 2, 4), (2, 2, 2, 3, 3), (2, 2, 3, 2, 3)})
        self.assertTrue(necklaces[(2, 2, 2, 2, 4)]["has_antipodal_contact"])
        self.assertTrue(necklaces[(2, 2, 2, 3, 3)]["has_antipodal_contact"])
        self.assertFalse(necklaces[(2, 2, 3, 2, 3)]["has_antipodal_contact"])
        self.assertEqual(len(gap["exact_cover_surviving_necklaces"]), 1)
        self.assertEqual(gap["exact_cover_surviving_necklaces"][0]["canonical_gap_necklace"], [2, 2, 3, 2, 3])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("x_0..x_11", proof["boolean_formulation"])
        self.assertIn("C(12,5)=792", proof["enumeration_domain"])
        self.assertIn("Exactly 12", proof["solution_count"])
        self.assertIn("one dihedral orbit", proof["orbit_certificate"])
        self.assertIn("[2,2,3,2,3]", proof["gap_certificate"])
        self.assertIn("three gap necklaces", proof["d1_only_boundary"])
        self.assertIn("conditional selector premise", proof["missing_source"])

        interpretation = self.payload["interpretation"]
        self.assertIn("select exactly the A5/d5", interpretation["direct_result"])
        self.assertIn("Boolean constraint certificate", interpretation["proof_style_upgrade"])
        self.assertIn("antipodal clauses", interpretation["d1_boundary"])
        self.assertIn("conditional/non-strict", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives exact-cover constraints", hard_limits)
        self.assertIn("conditional finite certificate", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
