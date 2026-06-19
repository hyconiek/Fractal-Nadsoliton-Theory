import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2911_s1861_gamma_localization_pullback_skeleton_gate.py"
JSON_PATH = ROOT / "generated" / "p2911_s1861_gamma_localization_pullback_skeleton_gate.json"
MD_PATH = ROOT / "generated" / "p2911_s1861_gamma_localization_pullback_skeleton_gate.md"


class P2911GammaLocalizationPullbackSkeletonGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_p2910_input(self):
        self.assertEqual(self.payload["status"], "P2911_GAMMA_LOCALIZATION_PULLBACK_SKELETON_GATE_READINESS_NO_EXPORT")
        self.assertTrue(self.acceptance["p2910_rechecked_dimension_readiness"])

    def test_finite_pullback_counts(self):
        self.assertEqual(self.acceptance["directed_edge_count"], 144)
        self.assertEqual(self.acceptance["site_count"], 12)
        self.assertEqual(self.acceptance["loop_edge_count"], 12)
        self.assertEqual(self.acceptance["nonloop_edge_count"], 132)
        self.assertEqual(len(self.objects["pullback_rows"]), 144)

    def test_pullback_properties(self):
        self.assertTrue(self.acceptance["all_column_sums_one"])
        self.assertTrue(self.acceptance["all_weights_nonnegative"])
        self.assertTrue(self.acceptance["endpoint_support_only"])
        self.assertEqual(self.acceptance["translation_equivariance_failure_count"], 0)
        self.assertTrue(self.acceptance["finite_localization_skeleton_constructed"])

    def test_no_ltotal_export(self):
        self.assertFalse(self.acceptance["strict_pullback_theorem_exported"])
        self.assertFalse(self.acceptance["variational_chain_rule_exported"])
        self.assertFalse(self.acceptance["accepted_as_nonproxy_ltotal_localization"])
        self.assertFalse(any(self.flags.values()))

    def test_documents_updated(self):
        self.assertIn("P2911/S1861", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2911/S1861", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2911/S1861", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2911", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
