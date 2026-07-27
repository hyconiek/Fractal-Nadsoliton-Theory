import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3172_s2122_legacy_star_operator_model_generator_potential_audit.py"
OUT = ROOT / "generated" / "p3172_s2122_legacy_star_operator_model_generator_potential_audit.json"
MD = ROOT / "generated" / "p3172_s2122_legacy_star_operator_model_generator_potential_audit.md"
PACKET = ROOT / "P3172_S2122_LEGACY_STAR_OPERATOR_MODEL_GENERATOR_POTENTIAL_AUDIT.md"


class P3172LegacyStarOperatorModelGeneratorPotentialAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_scope(self):
        self.assertEqual(self.payload["status"], "P3172_LEGACY_STAR_OPERATOR_MODEL_GENERATOR_POTENTIAL_AUDIT")
        scope = self.payload["scope"]
        self.assertTrue(scope["no_strict_bridge"])
        self.assertTrue(scope["no_parameter_fitting_to_strict"])
        self.assertTrue(scope["no_role_transfer"])
        self.assertTrue(scope["treat_as_research_program"])

    def test_operator_class_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["n_operator_classes"], 18)
        self.assertEqual(cert["n_TAK"] + cert["n_NIE"] + cert["n_WARUNKOWO"], 18)
        self.assertGreaterEqual(cert["n_TAK"], 3)
        self.assertGreaterEqual(cert["n_NIE"], 3)
        self.assertGreaterEqual(cert["n_WARUNKOWO"], 8)
        rows = {r["class"]: r["verdict"] for r in self.payload["constructed_theoretical_objects"]["operator_class_matrix"]}
        self.assertEqual(rows["graph_Laplacian"], "TAK")
        self.assertEqual(rows["diffusion"], "TAK")
        self.assertEqual(rows["integral"], "TAK")
        self.assertEqual(rows["spectral"], "TAK")
        self.assertEqual(rows["Dirac"], "NIE")
        self.assertEqual(rows["Maxwell"], "NIE")
        self.assertEqual(rows["Yang_Mills"], "NIE")

    def test_dual_dynamics_and_inverse(self):
        dual = self.payload["constructed_theoretical_objects"]["dual_dynamics_structure"]
        self.assertTrue(dual["common_generator"])
        self.assertTrue(dual["finite_witness"]["shared_spectral_measure"])
        self.assertLess(dual["finite_witness"]["unitarity_residual_fro"], 1e-8)
        inv = self.payload["constructed_theoretical_objects"]["inverse_recovery"]
        self.assertTrue(inv["partially_invertible_on_radial_circulant_class"])
        self.assertFalse(inv["surjective_onto_all_hermitian"])
        self.assertLess(inv["max_abs_recovery_error"], 1e-10)
        self.assertFalse(self.payload["finite_certificate"]["dI_is_metric_beta_0_01"])
        self.assertGreater(self.payload["finite_certificate"]["dI_violations_beta_0_01"], 0)
        self.assertFalse(self.payload["finite_certificate"]["dI_is_metric_beta1"])
        self.assertTrue(self.payload["finite_certificate"]["dI_triangle_holds_beta1"])
        self.assertAlmostEqual(
            self.payload["finite_certificate"]["beta_star_approx"],
            1.0751507762299,
            places=10,
        )

    def test_no_promotion_and_docs(self):
        dec = self.payload["decision"]
        self.assertTrue(dec["accepted_as_mathematical_generator"])
        self.assertFalse(dec["accepted_as_fundamental_physics"])
        self.assertFalse(dec["strict_bridge_performed"])
        self.assertFalse(dec["ToE_promoted"])
        self.assertFalse(dec["L_total_promoted"])
        self.assertTrue(all(v is False for v in dec["negative_export_flags"].values()))
        self.assertIn("P3172/S2122", MD.read_text(encoding="utf-8"))
        self.assertTrue(PACKET.exists())
        self.assertIn("P3172/S2122", (REPO / "AGENTS.md").read_text(encoding="utf-8"))
        self.assertIn("P3172/S2122", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        scores = {s["domain"]: s["score"] for s in self.payload["constructed_theoretical_objects"]["potential_scores"]}
        self.assertEqual(scores["unification_program"], 2)
        self.assertEqual(scores["operator_theory"], 8)


if __name__ == "__main__":
    unittest.main()
