import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2948_s1898_torsion_character_ratio_package_theorem_skeleton.py"
OUT = ROOT / "generated" / "p2948_s1898_torsion_character_ratio_package_theorem_skeleton.json"
MD = ROOT / "generated" / "p2948_s1898_torsion_character_ratio_package_theorem_skeleton.md"


class P2948TorsionCharacterRatioPackageTheoremSkeletonTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2948_TORSION_CHARACTER_RATIO_PACKAGE_THEOREM_SKELETON_FINITE_ONLY")
        self.assertIsNotNone(self.payload["input_hashes"]["P2938"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2945"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2947"])

    def test_package_certificate(self):
        cert = self.payload["package_certificate"]
        self.assertTrue(cert["finite_spine_passes"])
        self.assertTrue(cert["exact_p2938_vector_selected_finitely"])
        self.assertTrue(cert["sum9_selected_finitely"])
        self.assertTrue(cert["delta_4_5_constructed_finitely"])
        self.assertTrue(cert["eta_9_5_constructed_finitely"])
        self.assertTrue(cert["eta_equals_one_plus_delta"])
        self.assertFalse(cert["strict_torsion_character_source_theorem_exported"])
        self.assertFalse(cert["strict_delta_numerator_semantics_exported"])
        self.assertFalse(cert["strict_beta_eta_coupling_theorem_exported"])
        self.assertFalse(cert["accepted_strict_damping_source_package"])

    def test_spine_acceptance_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        spine = obj["finite_spine_rows"]
        self.assertTrue(all(row["satisfied_finitely"] for row in spine))
        self.assertEqual(spine[0]["value"], [1, 2, 2, 2, 2])
        self.assertEqual(spine[1]["value"], 9)
        self.assertEqual(spine[2]["value"]["as_string"], "4/5")
        self.assertEqual(spine[3]["value"]["as_string"], "9/5")
        acceptance = obj["acceptance_rows"]
        self.assertTrue(acceptance[0]["satisfied"])
        self.assertTrue(acceptance[1]["satisfied"])
        self.assertFalse(any(row["satisfied"] for row in acceptance[2:]))
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_torsion_character_source_theorem_exported", "strict_p2938_vector_source_theorem_exported", "strict_delta_eta_source_law_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2948/S1898", MD.read_text(encoding="utf-8"))
        self.assertIn("P2948/S1898", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2948/S1898", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2948", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
