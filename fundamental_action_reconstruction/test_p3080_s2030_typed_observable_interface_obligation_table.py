import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3080_s2030_typed_observable_interface_obligation_table.py"
OUT = ROOT / "generated" / "p3080_s2030_typed_observable_interface_obligation_table.json"
MD = ROOT / "generated" / "p3080_s2030_typed_observable_interface_obligation_table.md"

class P3080TypedObservableInterfaceObligationTableTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3080_TYPED_OBSERVABLE_INTERFACE_OBLIGATION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3079"])

    def test_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3079_accepted_internal_causal_order_sources"], 0)
        self.assertEqual(cert["interface_objects"], 4)
        self.assertEqual(cert["standard_physics_obligations"], 6)
        self.assertEqual(cert["obligation_cell_rows"], 24)
        self.assertEqual(cert["strictly_sourced_cells"], 0)
        self.assertEqual(cert["formal_imported_cells"], 10)
        self.assertEqual(cert["absent_cells"], 14)
        self.assertEqual(cert["accepted_standard_physics_interface_objects"], 0)
        self.assertEqual(cert["fully_discharged_obligations"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_object_and_obligation_summaries(self):
        objects = self.payload["constructed_theoretical_objects"]["interface_object_rows"]
        self.assertTrue(all(not row["standard_physics_interface_accepted"] for row in objects))
        by_id = {row["id"]: row for row in objects}
        self.assertEqual(by_id["z12_laplacian_smoothing_flow"]["formal_imported_obligations"], 2)
        self.assertEqual(by_id["imported_standard_physics_readout_template"]["formal_imported_obligations"], 6)
        summary = self.payload["constructed_theoretical_objects"]["obligation_summary_rows"]
        self.assertTrue(all(not row["discharged_by_current_artifacts"] for row in summary))

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["typed_interface_table_constructed"])
        self.assertIn("dimension/action-unit source audit", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3080/S2030", MD.read_text(encoding="utf-8"))
        self.assertIn("P3080/S2030", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3080/S2030", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3080", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
