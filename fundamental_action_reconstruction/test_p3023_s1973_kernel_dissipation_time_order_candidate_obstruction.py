import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3023_s1973_kernel_dissipation_time_order_candidate_obstruction.py"
OUT = ROOT / "generated" / "p3023_s1973_kernel_dissipation_time_order_candidate_obstruction.json"
MD = ROOT / "generated" / "p3023_s1973_kernel_dissipation_time_order_candidate_obstruction.md"

class P3023KernelDissipationTimeOrderCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3023_KERNEL_DISSIPATION_TIME_ORDER_CANDIDATE_OBSTRUCTION_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3022"])

    def test_descent_positive_but_equivariance_failure(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["chain_edge_count"], 11)
        self.assertEqual(cert["strict_descent_chain_edges"], 11)
        self.assertFalse(cert["cyclic_reset_strict_descent"])
        self.assertEqual(cert["unit_equivariant_rows"], 1)
        self.assertEqual(cert["unit_row_count"], 4)
        self.assertFalse(cert["accepted_as_strict_time_order_with_physical_unit"])

    def test_obligations_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "KernelDissipationTimeOrderCandidate_EquivarianceUnitObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["new_typed_time_order_object"])
        self.assertTrue(obligations["strict_descent_chain"])
        self.assertFalse(obligations["no_cyclic_reset_obstruction"])
        self.assertFalse(obligations["U12_equivariant_directed_successor"])
        self.assertFalse(obligations["strict_chart_or_selector_source"])
        self.assertFalse(obligations["physical_unit_theorem"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3023/S1973", MD.read_text(encoding="utf-8"))
        self.assertIn("P3023/S1973", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3023/S1973", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3023", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
