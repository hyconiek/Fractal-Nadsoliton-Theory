import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2919_s1869_gamma_scale_orbit_source_object_intake_gate.py"
OUT = ROOT / "generated" / "p2919_s1869_gamma_scale_orbit_source_object_intake_gate.json"
MD = ROOT / "generated" / "p2919_s1869_gamma_scale_orbit_source_object_intake_gate.md"


class P2919GammaScaleOrbitSourceObjectIntakeGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_p2918_input(self):
        self.assertEqual(self.payload["status"], "P2919_GAMMA_SCALE_ORBIT_SOURCE_OBJECT_INTAKE_GATE_NO_EXPORT")
        self.assertTrue(self.payload["acceptance_matrix"]["p2918_rechecked_gamma_free_scalar"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2918"])

    def test_source_object_schema(self):
        schema = self.payload["constructed_theoretical_objects"]["required_source_object_schema"]
        self.assertEqual(schema["object_name"], "Strict_Gamma_9_5_Action_Unit_Scale_Fixing_Source_Object")
        self.assertIn("sigma_Gamma", schema["source_map"])
        self.assertEqual(len(schema["acceptance_obligations"]), 5)

    def test_scale_orbit_matrix(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["finite_observable_count"], 5)
        self.assertTrue(acc["finite_observables_fix_relative_weights"])
        self.assertFalse(acc["finite_observables_fix_gamma_scale"])
        self.assertFalse(acc["strict_sigma_gamma_source_object_exported"])

    def test_candidate_intake_counts(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["source_object_candidate_count"], 6)
        self.assertEqual(acc["scale_breaking_candidate_count"], 1)
        self.assertEqual(acc["strict_scale_breaking_candidate_count"], 0)
        self.assertEqual(acc["accepted_source_object_count"], 0)

    def test_false_closure_exports_and_docs(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["nonproxy_ltotal_exported", "eom_closure_exported", "hamiltonian_closure_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2919/S1869", MD.read_text(encoding="utf-8"))
        self.assertIn("P2919/S1869", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2919/S1869", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2919", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
