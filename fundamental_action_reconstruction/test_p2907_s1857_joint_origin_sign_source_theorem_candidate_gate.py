import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2907_s1857_joint_origin_sign_source_theorem_candidate_gate.py"
JSON_PATH = ROOT / "generated" / "p2907_s1857_joint_origin_sign_source_theorem_candidate_gate.json"
MD_PATH = ROOT / "generated" / "p2907_s1857_joint_origin_sign_source_theorem_candidate_gate.md"


class P2907JointOriginSignSourceTheoremCandidateGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2907_JOINT_ORIGIN_SIGN_SOURCE_THEOREM_CANDIDATE_GATE_NO_STRICT_EXPORT")
        self.assertTrue(self.acceptance["p2906_rechecked_joint_theorem_required"])
        self.assertTrue(self.acceptance["p2906_rechecked_joint_theorem_not_exported"])

    def test_joint_candidate_collapses_xi_family(self):
        self.assertEqual(self.acceptance["xi_family_count_before_joint_postulate"], 24)
        self.assertEqual(self.acceptance["xi_family_count_after_joint_postulate"], 1)
        self.assertEqual(self.acceptance["selected_origin"], 0)
        self.assertEqual(self.acceptance["selected_sign"], 1)
        self.assertEqual(self.acceptance["selected_defect_edge"], [0, 5])
        self.assertEqual(self.objects["selected_xi_target"]["xi"], "Xi_{0,+}")

    def test_symbolic_variational_derivative(self):
        candidate = self.objects["joint_origin_sign_candidate"]
        derivatives = candidate["local_variational_derivative"]
        self.assertEqual(self.acceptance["nonzero_local_derivative_count"], 1)
        self.assertEqual(self.acceptance["zero_local_derivative_count"], 143)
        self.assertEqual(derivatives["q_9_5(0,5)"], "U_9_5")
        self.assertEqual(derivatives["q_9_5(0,6)"], 0)

    def test_no_strict_export_or_ltotal(self):
        self.assertTrue(self.acceptance["passes_p2906_joint_origin_sign_finite_gate"])
        self.assertFalse(self.acceptance["strict_nadsoliton_provenance_exported"])
        self.assertFalse(self.acceptance["unit_bearing_nonproxy_ltotal_coupling_exported"])
        self.assertFalse(self.acceptance["accepted_as_strict_source_theorem"])
        self.assertFalse(any(self.flags.values()))

    def test_documents_updated(self):
        self.assertIn("P2907/S1857", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2907/S1857", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2907/S1857", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2907", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
