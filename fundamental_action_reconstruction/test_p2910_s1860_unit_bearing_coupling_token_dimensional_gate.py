import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2910_s1860_unit_bearing_coupling_token_dimensional_gate.py"
JSON_PATH = ROOT / "generated" / "p2910_s1860_unit_bearing_coupling_token_dimensional_gate.json"
MD_PATH = ROOT / "generated" / "p2910_s1860_unit_bearing_coupling_token_dimensional_gate.md"


class P2910UnitBearingCouplingTokenDimensionalGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_p2909_input(self):
        self.assertEqual(self.payload["status"], "P2910_UNIT_BEARING_COUPLING_TOKEN_DIMENSIONAL_GATE_READINESS_NO_EXPORT")
        self.assertTrue(self.acceptance["p2909_rechecked_no_live_frontier"])

    def test_new_typed_object_and_dimension_gate(self):
        self.assertTrue(self.acceptance["new_typed_object_constructed"])
        self.assertTrue(self.acceptance["outside_xi_j_defect_placement_lane"])
        self.assertTrue(self.acceptance["dimension_algebra_passed"])
        self.assertEqual(self.acceptance["required_gamma_dimension"], {"Action": 1, "Length": 0, "Time": 0})
        token = self.objects["coupling_token_candidate"]
        self.assertEqual(token["name"], "Gamma_9_5")
        self.assertEqual(token["reconstructed_density_dimension"], {"Action": 1, "Length": -1, "Time": 0})

    def test_missing_ltotal_theorems(self):
        self.assertFalse(self.acceptance["strict_gamma_source_exported"])
        self.assertFalse(self.acceptance["localization_pullback_exported"])
        self.assertFalse(self.acceptance["variational_chain_rule_exported"])
        self.assertFalse(self.acceptance["compatible_with_j_provenance"])
        self.assertFalse(self.acceptance["accepted_as_unit_bearing_ltotal_coupling"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.flags["strict_gamma_9_5_source_exported"])
        self.assertFalse(self.flags["unit_bearing_nonproxy_ltotal_exported"])
        self.assertFalse(self.flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2910/S1860", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2910/S1860", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2910/S1860", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2910", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
