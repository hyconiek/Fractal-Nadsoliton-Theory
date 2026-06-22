import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3021_s1971_action_quantum_reference_cell_source_obstruction.py"
OUT = ROOT / "generated" / "p3021_s1971_action_quantum_reference_cell_source_obstruction.json"
MD = ROOT / "generated" / "p3021_s1971_action_quantum_reference_cell_source_obstruction.md"

class P3021ActionQuantumReferenceCellSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3021_ACTION_QUANTUM_REFERENCE_CELL_SOURCE_OBSTRUCTION_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3020"])

    def test_positive_partitions_but_rescaling_obstruction(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["reference_cell_candidate_count"], 4)
        self.assertEqual(cert["positive_reference_cell_candidates"], 4)
        self.assertEqual(cert["scale_row_count"], 4)
        self.assertEqual(cert["stable_cell_count_rows"], 4)
        self.assertEqual(cert["quanta_rescale_rows"], 4)
        self.assertFalse(cert["accepted_as_strict_action_quantum_reference_cell_source"])

    def test_obligations_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "ActionQuantumReferenceCellSource_RescalingPartitionObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["typed_formal_action_input"])
        self.assertTrue(obligations["positive_finite_reference_cell_candidates"])
        self.assertTrue(obligations["cell_counts_stable_under_observable_rescaling"])
        self.assertFalse(obligations["action_quantum_invariant_under_observable_rescaling"])
        self.assertFalse(obligations["strict_reference_cell_source_theorem"])
        self.assertFalse(obligations["physical_action_unit_and_hamiltonian_coupling"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3021/S1971", MD.read_text(encoding="utf-8"))
        self.assertIn("P3021/S1971", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3021/S1971", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3021", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
