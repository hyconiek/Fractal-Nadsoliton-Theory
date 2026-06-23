import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3050_s2000_area_sign_to_p3046_polarity_coupling_audit.py"
OUT = ROOT / "generated" / "p3050_s2000_area_sign_to_p3046_polarity_coupling_audit.json"
MD = ROOT / "generated" / "p3050_s2000_area_sign_to_p3046_polarity_coupling_audit.md"

class P3050AreaSignToP3046PolarityCouplingAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3050_AREA_SIGN_TO_P3046_POLARITY_COUPLING_AUDIT_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3049"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["area_sign_rows"], 6)
        self.assertEqual(cert["nonzero_area_sign_rows"], 5)
        self.assertEqual(cert["neutral_area_rows"], 1)
        self.assertEqual(cert["coupling_rows"], 11)
        self.assertEqual(cert["nonzero_candidate_coupling_rows"], 10)
        self.assertEqual(cert["aut_equivariant_nonzero_coupling_rows"], 10)
        self.assertEqual(cert["polarity_selected_rows"], 0)
        self.assertEqual(cert["accepted_orientation_coupling_rows"], 0)
        self.assertEqual(cert["source_acceptance_criteria"], 7)
        self.assertEqual(cert["satisfied_source_acceptance_criteria"], 3)
        self.assertFalse(cert["p3046_coupling_polarity_selected"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "AreaSignToP3046Polarity_CouplingAudit")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["p3049_remaining_coupling_premise_targeted"])
        self.assertTrue(obligations["area_sign_torsor_constructed"])
        self.assertTrue(obligations["all_nonzero_equivariant_maps_enumerated"])
        self.assertFalse(obligations["unique_p3046_coupling_polarity"])
        self.assertFalse(obligations["nonconventional_orientation_theorem"])
        self.assertFalse(obligations["selector_readout_or_ltotal_installation"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3050/S2000", MD.read_text(encoding="utf-8"))
        self.assertIn("P3050/S2000", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3050/S2000", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3050", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
