import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2944_s1894_partial_monoid_identity_grounding_theorem_audit.py"
OUT = ROOT / "generated" / "p2944_s1894_partial_monoid_identity_grounding_theorem_audit.json"
MD = ROOT / "generated" / "p2944_s1894_partial_monoid_identity_grounding_theorem_audit.md"


class P2944PartialMonoidIdentityGroundingTheoremAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2944_PARTIAL_MONOID_IDENTITY_GROUNDING_THEOREM_AUDIT_NO_STRICT_SOURCE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2943"])

    def test_identity_grounding_certificate(self):
        cert = self.payload["identity_grounding_certificate"]
        self.assertEqual(cert["node_count"], 11)
        self.assertEqual(cert["partial_product_count"], 29)
        self.assertEqual(cert["two_sided_total_identity_nodes"], [1])
        self.assertEqual(cert["zero_carrier_nodes"], [1])
        self.assertTrue(cert["identity_equals_unique_zero"])
        self.assertTrue(cert["finite_identity_grounding_verified"])
        self.assertFalse(cert["strict_identity_grounding_theorem_exported"])
        self.assertFalse(cert["accepted_strict_source"])

    def test_theorem_rows_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(len(obj["identity_candidate_rows"]), 11)
        self.assertEqual(len(obj["partial_product_rows"]), 29)
        self.assertTrue(all(row["satisfied"] for row in obj["algebraic_theorem_rows"]))
        strict_rows = obj["strict_acceptance_rows"]
        self.assertTrue(strict_rows[0]["satisfied"])
        for row in strict_rows[1:]:
            self.assertFalse(row["satisfied"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_identity_grounding_theorem_exported", "strict_positive_orientation_source_theorem_exported", "strict_selector_source_exported", "strict_aut_breaking_prime_coordinate_source_law_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2944/S1894", MD.read_text(encoding="utf-8"))
        self.assertIn("P2944/S1894", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2944/S1894", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2944", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
