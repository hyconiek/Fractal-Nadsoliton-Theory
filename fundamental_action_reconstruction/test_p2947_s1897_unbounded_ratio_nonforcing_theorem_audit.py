import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2947_s1897_unbounded_ratio_nonforcing_theorem_audit.py"
OUT = ROOT / "generated" / "p2947_s1897_unbounded_ratio_nonforcing_theorem_audit.json"
MD = ROOT / "generated" / "p2947_s1897_unbounded_ratio_nonforcing_theorem_audit.md"


class P2947UnboundedRatioNonforcingTheoremAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2947_UNBOUNDED_RATIO_NONFORCING_THEOREM_AUDIT_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2945"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2946"])

    def test_unbounded_nonforcing_certificate(self):
        cert = self.payload["unbounded_nonforcing_certificate"]
        self.assertTrue(cert["parametric_family_symbolically_verified"])
        self.assertEqual(cert["sample_sum_count"], 11)
        self.assertEqual(cert["sample_non_target_eta_count"], 10)
        self.assertFalse(cert["eta_9_5_forced_by_positive_cone_premises"])
        self.assertFalse(cert["delta_formula_uniquely_selected_by_counts"])
        self.assertFalse(cert["strict_p2938_vector_sum9_source_theorem_exported"])
        self.assertFalse(cert["strict_beta_eta_coupling_theorem_exported"])
        self.assertFalse(cert["accepted_strict_ratio_source_theorem"])

    def test_theorem_rows_and_alias_kernel(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertTrue(obj["p2946_bounded_result_subsumed"])
        theorem_rows = obj["theorem_rows"]
        self.assertTrue(all(row["verified_symbolically"] for row in theorem_rows))
        eta_rows = obj["parametric_eta_witness_rows"]
        self.assertEqual(eta_rows[0]["witness_vector"], [1, 1, 1, 1, 1])
        self.assertEqual(eta_rows[4]["witness_vector"], [5, 1, 1, 1, 1])
        self.assertEqual(eta_rows[4]["eta"]["as_string"], "9/5")
        self.assertTrue(eta_rows[4]["matches_target_eta_9_5"])
        self.assertEqual(sum(row["matches_target_eta_9_5"] for row in eta_rows), 1)
        delta_rows = obj["delta_alias_kernel_rows"]
        self.assertEqual(sum(row["matches_delta_4_5"] for row in delta_rows), len(delta_rows))
        self.assertTrue(all(not row["strict_semantics_selected"] for row in delta_rows))

    def test_acceptance_nonpromotion_and_docs(self):
        obj = self.payload["constructed_theoretical_objects"]
        acceptance = obj["acceptance_rows"]
        self.assertTrue(acceptance[0]["satisfied"])
        self.assertFalse(any(row["satisfied"] for row in acceptance[1:]))
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_ratio_source_theorem_exported", "strict_p2938_vector_source_theorem_exported", "strict_delta_eta_source_law_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2947/S1897", MD.read_text(encoding="utf-8"))
        self.assertIn("P2947/S1897", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2947/S1897", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2947", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
