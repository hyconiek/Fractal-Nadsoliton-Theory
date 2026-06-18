import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2843_s1793_professorial_swot_foundations_kernel_lagrangian_eom_hamiltonian_closure_path.py"
JSON_PATH = ROOT / "generated" / "p2843_s1793_professorial_swot_foundations_kernel_lagrangian_eom_hamiltonian_closure_path.json"
MD_PATH = ROOT / "generated" / "p2843_s1793_professorial_swot_foundations_kernel_lagrangian_eom_hamiltonian_closure_path.md"


class P2843ProfessorialSWOTClosurePathAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status(self):
        self.assertEqual(self.payload["status"], "P2843_PROFESSORIAL_SWOT_CLOSURE_PATH_AUDIT_NO_NEW_CLOSURE")
        self.assertEqual(
            self.payload["input_statuses_rechecked"]["P2842"],
            "P2842_EXCHANGEABLE_EDGE_PAIR_MEASURE_LOCALIZATION_CANDIDATE_NO_GO_NO_COUPLING_NO_CLOSURE",
        )

    def test_swot_domains_and_content(self):
        swot = self.payload["swot"]
        self.assertEqual(set(swot), {"foundations", "kernel_formulas", "lagrangian_eom_hamiltonian"})
        for domain in swot.values():
            for bucket in ("strengths", "weaknesses", "opportunities", "threats"):
                self.assertGreaterEqual(len(domain[bucket]), 2)
        self.assertIn("K_strict_gate", " ".join(swot["kernel_formulas"]["strengths"]))
        self.assertIn("Hamiltonian", " ".join(swot["lagrangian_eom_hamiltonian"]["weaknesses"]))

    def test_closure_path_boundaries(self):
        rows = self.payload["closure_path"]
        self.assertEqual(len(rows), 8)
        self.assertEqual(rows[0]["name"], "Freeze scope and notation")
        self.assertEqual(rows[1]["name"], "Kernel completion-map theorem")
        self.assertEqual(rows[6]["name"], "Hamiltonian/constraint closure")
        self.assertFalse(rows[6]["admissible_now"])
        self.assertIn("Legendre", rows[6]["professorial_action"])

    def test_acceptance_matrix_and_negative_exports(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["accepted_as_professorial_swot_and_closure_path_audit"])
        self.assertFalse(acceptance["exports_ltotal_or_toe_closure"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["strict_kernel_identity_exported"])
        self.assertFalse(flags["role_bearing_ltotal_promoted"])
        self.assertFalse(flags["eom_closure_exported"])
        self.assertFalse(flags["hamiltonian_closure_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_source_metrics_and_documents_updated(self):
        metrics = self.payload["source_metrics"]
        self.assertIn("K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md", metrics["files_accounted"])
        self.assertGreater(metrics["mention_counts"]["K_strict_gate"], 0)
        self.assertGreater(metrics["mention_counts"]["L_total_or_lagrangian"], 0)
        self.assertIn("P2843/S1793", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2843/S1793", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2843/S1793", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2843", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
