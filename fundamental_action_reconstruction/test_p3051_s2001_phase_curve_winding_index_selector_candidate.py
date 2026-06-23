import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3051_s2001_phase_curve_winding_index_selector_candidate.py"
OUT = ROOT / "generated" / "p3051_s2001_phase_curve_winding_index_selector_candidate.json"
MD = ROOT / "generated" / "p3051_s2001_phase_curve_winding_index_selector_candidate.md"

class P3051PhaseCurveWindingIndexSelectorCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3051_PHASE_CURVE_WINDING_INDEX_SELECTOR_CANDIDATE_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3050"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["base_integer_winding"], 1)
        self.assertTrue(cert["base_nonzero_integer_winding"])
        self.assertEqual(cert["aut_translation_rows"], 48)
        self.assertEqual(cert["translation_stable_rows"], 48)
        self.assertEqual(cert["orientation_reversing_rows"], 24)
        self.assertEqual(cert["strict_source_exported_rows"], 0)
        self.assertEqual(cert["source_acceptance_criteria"], 8)
        self.assertEqual(cert["satisfied_source_acceptance_criteria"], 4)
        self.assertFalse(cert["strict_winding_selector_source_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "PhaseCurveWindingIndex_SelectorCandidateAudit")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["new_typed_object_not_area_replay"])
        self.assertTrue(obligations["finite_nonzero_topological_hint"])
        self.assertTrue(obligations["aut_translation_witness_matrix"])
        self.assertFalse(obligations["strict_winding_orientation_source"])
        self.assertFalse(obligations["selector_readout_or_p3046_coupling"])
        self.assertFalse(obligations["ltotal_bridge_role_toe_installation"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3051/S2001", MD.read_text(encoding="utf-8"))
        self.assertIn("P3051/S2001", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3051/S2001", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3051", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
