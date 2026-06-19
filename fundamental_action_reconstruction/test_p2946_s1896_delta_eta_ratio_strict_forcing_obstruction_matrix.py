import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2946_s1896_delta_eta_ratio_strict_forcing_obstruction_matrix.py"
OUT = ROOT / "generated" / "p2946_s1896_delta_eta_ratio_strict_forcing_obstruction_matrix.json"
MD = ROOT / "generated" / "p2946_s1896_delta_eta_ratio_strict_forcing_obstruction_matrix.md"


class P2946DeltaEtaRatioStrictForcingObstructionMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2946_DELTA_ETA_RATIO_STRICT_FORCING_OBSTRUCTION_MATRIX_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2945"])

    def test_strict_forcing_certificate(self):
        cert = self.payload["strict_forcing_certificate"]
        self.assertTrue(cert["finite_model_class_enumerated"])
        self.assertFalse(cert["eta_9_5_forced_by_positive_cone_premises"])
        self.assertFalse(cert["p2938_vector_uniquely_forced"])
        self.assertFalse(cert["delta_formula_uniquely_selected"])
        self.assertFalse(cert["strict_ratio_source_theorem_exported"])
        self.assertFalse(cert["strict_beta_eta_coupling_theorem_exported"])
        self.assertFalse(cert["accepted_strict_delta_eta_source_law"])

    def test_eta_model_class_and_delta_aliases(self):
        obj = self.payload["constructed_theoretical_objects"]
        summary = obj["eta_variant_summary"]
        self.assertEqual(summary["variant_count"], 243)
        self.assertEqual(summary["eta_9_5_match_count"], 45)
        self.assertEqual(summary["eta_9_5_nonmatch_count"], 198)
        self.assertTrue(summary["p2938_vector_present"])
        aliases = obj["delta_formula_alias_rows"]
        self.assertGreaterEqual(sum(row["matches_p2945_delta_4_5"] for row in aliases), 2)
        self.assertTrue(all(not row["selected_by_strict_theorem"] for row in aliases))
        obligations = obj["forcing_obligation_rows"]
        self.assertTrue(obligations[0]["satisfied"])
        self.assertFalse(obligations[1]["satisfied"])
        self.assertFalse(obligations[2]["satisfied"])
        self.assertFalse(obligations[3]["satisfied"])
        self.assertFalse(obligations[4]["satisfied"])

    def test_nonpromotion_and_docs_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_ratio_source_theorem_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_law_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2946/S1896", MD.read_text(encoding="utf-8"))
        self.assertIn("P2946/S1896", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2946/S1896", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2946", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
