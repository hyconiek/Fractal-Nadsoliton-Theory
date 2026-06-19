import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2915_s1865_gamma_edge_quotient_measure_bridge_candidate_gate.py"
JSON_PATH = ROOT / "generated" / "p2915_s1865_gamma_edge_quotient_measure_bridge_candidate_gate.json"
MD_PATH = ROOT / "generated" / "p2915_s1865_gamma_edge_quotient_measure_bridge_candidate_gate.md"


class P2915GammaEdgeQuotientMeasureBridgeCandidateGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_p2914_input(self):
        self.assertEqual(self.payload["status"], "P2915_GAMMA_EDGE_QUOTIENT_MEASURE_BRIDGE_CANDIDATE_GATE_NO_EXPORT")
        self.assertTrue(self.acceptance["p2914_rechecked_measure_obstruction"])

    def test_quotient_candidates_resolve_arithmetic(self):
        self.assertEqual(self.acceptance["quotient_candidate_count"], 3)
        self.assertEqual(self.acceptance["arithmetic_normalization_resolving_candidate_count"], 3)
        self.assertEqual(self.acceptance["strictly_selected_quotient_count"], 0)
        self.assertEqual(self.acceptance["site_normalized_m"], "1/12")
        self.assertEqual(self.acceptance["per_edge_quotient_weight"], "1/144")
        self.assertEqual(self.acceptance["quotient_edge_total"], "1/1")

    def test_candidate_fibers(self):
        rows = self.objects["quotient_candidate_rows"]
        self.assertEqual(len(rows), 3)
        for row in rows:
            self.assertEqual(row["fiber_count"], 12)
            self.assertTrue(row["all_fibers_size_12"])
            self.assertTrue(row["normalization_mismatch_resolved_arithmetically"])
            self.assertFalse(row["strictly_selected"])

    def test_no_ltotal_or_closure_export(self):
        self.assertFalse(self.acceptance["unique_strict_quotient_exported"])
        self.assertFalse(self.acceptance["strict_continuum_integration_theorem_exported"])
        self.assertFalse(self.acceptance["accepted_as_nonproxy_ltotal_measure_bridge"])
        self.assertFalse(any(self.flags.values()))

    def test_documents_updated(self):
        self.assertIn("P2915/S1865", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2915/S1865", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2915/S1865", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2915", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
