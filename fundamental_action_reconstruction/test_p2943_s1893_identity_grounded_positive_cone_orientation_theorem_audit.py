import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2943_s1893_identity_grounded_positive_cone_orientation_theorem_audit.py"
OUT = ROOT / "generated" / "p2943_s1893_identity_grounded_positive_cone_orientation_theorem_audit.json"
MD = ROOT / "generated" / "p2943_s1893_identity_grounded_positive_cone_orientation_theorem_audit.md"


class P2943IdentityGroundedPositiveConeOrientationTheoremAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(cls_status := self.payload["status"], "P2943_IDENTITY_GROUNDED_POSITIVE_CONE_ORIENTATION_THEOREM_CANDIDATE_NO_STRICT_SOURCE")
        self.assertIn("POSITIVE_CONE", cls_status)
        self.assertIsNotNone(self.payload["input_hashes"]["P2942"])

    def test_positive_cone_certificate(self):
        cert = self.payload["positive_cone_certificate"]
        self.assertEqual(cert["node_count"], 11)
        self.assertEqual(cert["identity_value"], 0)
        self.assertEqual(cert["positive_prime_coordinate_count"], 5)
        self.assertEqual(cert["positive_nonidentity_node_count"], 10)
        self.assertEqual(cert["product_additivity_defect_count"], 0)
        self.assertTrue(cert["finite_positive_cone_theorem_candidate_verified"])
        self.assertFalse(cert["strict_identity_grounding_theorem_exported"])
        self.assertFalse(cert["strict_provenance_theorem_exported"])
        self.assertFalse(cert["accepted_strict_source"])

    def test_premises_acceptance_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["prime_vector_order_2_3_5_7_11"], [1, 2, 2, 2, 2])
        self.assertTrue(all(row["satisfies_positive_cone"] for row in obj["node_value_rows"]))
        self.assertTrue(all(row["satisfied"] for row in obj["algebraic_theorem_premise_rows"]))
        strict_rows = obj["strict_acceptance_rows"]
        self.assertTrue(strict_rows[0]["satisfied"])
        for row in strict_rows[1:]:
            self.assertFalse(row["satisfied"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_identity_grounding_theorem_exported", "strict_positive_orientation_source_theorem_exported", "strict_selector_source_exported", "strict_aut_breaking_prime_coordinate_source_law_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2943/S1893", MD.read_text(encoding="utf-8"))
        self.assertIn("P2943/S1893", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2943/S1893", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2943", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
