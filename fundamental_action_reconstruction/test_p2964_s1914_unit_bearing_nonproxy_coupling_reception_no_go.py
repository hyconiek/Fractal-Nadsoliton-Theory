import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2964_s1914_unit_bearing_nonproxy_coupling_reception_no_go.py"
OUT = ROOT / "generated" / "p2964_s1914_unit_bearing_nonproxy_coupling_reception_no_go.json"
MD = ROOT / "generated" / "p2964_s1914_unit_bearing_nonproxy_coupling_reception_no_go.md"

class P2964UnitBearingNonproxyCouplingReceptionNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_inputs_and_constants(self):
        self.assertEqual(self.payload["status"], "P2964_UNIT_BEARING_NONPROXY_COUPLING_RECEPTION_NO_GO_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2961"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2963"])
        self.assertEqual(self.payload["package_constants"]["aggregate_vector"], [1, 2, 2, 2, 2])
        self.assertEqual(self.payload["package_constants"]["sum"], 9)
        self.assertAlmostEqual(self.payload["package_constants"]["eta"], 1.8)

    def test_coupling_rows(self):
        cert = self.payload["coupling_certificate"]
        self.assertEqual(cert["template_count"], 6)
        self.assertEqual(cert["accepted_current_strict_couplings"], [])
        self.assertFalse(cert["strict_nonproxy_coupling_exported"])
        rows = {r["template"]: r for r in self.payload["constructed_theoretical_objects"]["coupling_template_rows"]}
        self.assertTrue(rows["aggregate_norm_action_density"]["receives_aggregate"])
        self.assertFalse(rows["aggregate_norm_action_density"]["requires_KC_artifact_exchange"])
        self.assertFalse(rows["aggregate_norm_action_density"]["canonical_scale_quotient"])
        self.assertFalse(rows["strict_completed_coupling_schema"]["current_artifact_available"])
        self.assertIn("current_artifact_export", rows["strict_completed_coupling_schema"]["missing_obligations"])

    def test_obligations_and_matrix(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["aggregate_reception_without_KC_exchange_representable"])
        self.assertFalse(obligations["canonical_scale_quotient_exported"])
        self.assertTrue(obligations["unit_bearing_nonproxy_density_exported"])
        self.assertTrue(obligations["variational_derivative_exported"])
        self.assertFalse(obligations["strict_coefficient_source_law_exported"])
        self.assertFalse(obligations["accepted_current_strict_coupling"])
        matrix = self.payload["constructed_theoretical_objects"]["finite_acceptance_matrix"]
        self.assertEqual(len(matrix), 128)
        self.assertEqual(sum(1 for row in matrix if row["accepts_strict_nonproxy_coupling"]), 1)

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2964/S1914", MD.read_text(encoding="utf-8"))
        self.assertIn("P2964/S1914", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2964/S1914", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2964", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
