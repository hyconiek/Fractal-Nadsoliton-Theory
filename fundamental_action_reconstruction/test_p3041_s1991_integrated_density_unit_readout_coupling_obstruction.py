import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3041_s1991_integrated_density_unit_readout_coupling_obstruction.py"
OUT = ROOT / "generated" / "p3041_s1991_integrated_density_unit_readout_coupling_obstruction.json"
MD = ROOT / "generated" / "p3041_s1991_integrated_density_unit_readout_coupling_obstruction.md"

class P3041IntegratedDensityUnitReadoutCouplingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3041_INTEGRATED_DENSITY_UNIT_READOUT_COUPLING_OBSTRUCTION_NO_SOURCE_EXPORT")
        for key in ["P3036", "P3038", "P3039", "P3040"]:
            self.assertIsNotNone(self.payload["input_hashes"][key])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["receiver_rows"], 4)
        self.assertEqual(cert["finite_nonzero_rows"], 4)
        self.assertEqual(cert["accepted_physical_unit_readout_coupling_rows"], 0)
        self.assertEqual(cert["unit_grid_rows"], 5)
        self.assertEqual(cert["unit_scaling_exact_rows"], 5)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)
        self.assertFalse(cert["physical_unit_readout_coupling_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "IntegratedDensityUnitReadoutCoupling_ObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["exact_p3038_unit_readout_premise_targeted"])
        self.assertTrue(obligations["finite_integrated_density_readout_nonzero"])
        self.assertTrue(obligations["unit_scale_orbit_constructed"])
        self.assertTrue(obligations["p3036_unit_boundary_consulted"])
        self.assertFalse(obligations["external_physical_unit_source"])
        self.assertFalse(obligations["scale_orbit_quotient_or_absolute_calibration"])
        self.assertFalse(obligations["readout_coupling_theorem"])
        self.assertFalse(obligations["nonproxy_variational_or_action_installation"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3041/S1991", MD.read_text(encoding="utf-8"))
        self.assertIn("P3041/S1991", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3041/S1991", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3041", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
