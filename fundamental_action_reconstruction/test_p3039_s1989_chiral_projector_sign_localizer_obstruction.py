import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3039_s1989_chiral_projector_sign_localizer_obstruction.py"
OUT = ROOT / "generated" / "p3039_s1989_chiral_projector_sign_localizer_obstruction.json"
MD = ROOT / "generated" / "p3039_s1989_chiral_projector_sign_localizer_obstruction.md"

class P3039ChiralProjectorSignLocalizerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3039_CHIRAL_PROJECTOR_SIGN_LOCALIZER_OBSTRUCTION_NO_SOURCE_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3038"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["receiver_rows"], 4)
        self.assertEqual(cert["accepted_nonpremise_localizer_rows"], 0)
        self.assertEqual(cert["translation_phase_orbit_size"], 12)
        self.assertEqual(cert["phase_polarity_projector_count"], 12)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)
        self.assertFalse(cert["nonpremise_chiral_sign_localizer_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "ChiralProjectorSignLocalizer_ObstructionMatrix")
        self.assertTrue(obj["torsor"]["inversion_unit_11_flips_polarity"])
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["exact_p3038_missing_premise_targeted"])
        self.assertTrue(obligations["finite_chiral_projector_torsor_constructed"])
        self.assertTrue(obligations["inversion_odd_sign_verified"])
        self.assertTrue(obligations["current_hint_receivers_tested"])
        self.assertFalse(obligations["nonpremise_phase_origin_localizer"])
        self.assertFalse(obligations["nonpremise_polarity_sign_source"])
        self.assertFalse(obligations["aut_z12_compatible_selector_source"])
        self.assertFalse(obligations["coupling_back_to_p3038_as_source"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3039/S1989", MD.read_text(encoding="utf-8"))
        self.assertIn("P3039/S1989", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3039/S1989", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3039", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
