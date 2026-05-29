import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_unit_mirror_tie_obstruction_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_unit_mirror_tie_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_unit_mirror_tie_obstruction_report.md"


class StrictAlphaHebbianUnitMirrorTieObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_UNIT_MIRROR_TIE_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "full-aut-invariance-cannot-break-k1-k5-unit-mirror-tie")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_score_signature_tie(self):
        signature = self.payload["score_signature"]
        self.assertEqual(signature["k1_contiguous"]["power_exact"], {"rational_part": "7", "sqrt3_coefficient": "4", "expression": "7 + 4*sqrt3"})
        self.assertEqual(signature["k5_d5"]["power_exact"], {"rational_part": "7", "sqrt3_coefficient": "4", "expression": "7 + 4*sqrt3"})
        self.assertEqual(signature["k1_contiguous"]["nyquist_power_exact"], "1")
        self.assertEqual(signature["k5_d5"]["nyquist_power_exact"], "1")
        self.assertTrue(signature["same_mode_power"])
        self.assertTrue(signature["same_nyquist_power"])
        self.assertTrue(signature["same_for_any_scalar_anti_nyquist_score"])

    def test_full_aut_action_certificate(self):
        cert = self.payload["full_aut_action_certificate"]
        self.assertEqual(cert["orbit_size"], 2)
        self.assertEqual(cert["k1_to_k5_witness_units"], [5, 7])
        self.assertEqual(cert["k5_to_k1_witness_units"], [5, 7])
        self.assertEqual(cert["k5_stabilizer"], [1, 11])
        self.assertEqual(cert["k1_stabilizer"], [1, 11])
        self.assertFalse(cert["singleton_d5_full_aut_invariant"])
        invariant_sizes = sorted(len(subset) for subset in cert["full_aut_invariant_subsets_on_orbit"])
        self.assertEqual(invariant_sizes, [0, 2])

    def test_subgroup_reduction_audit(self):
        audit = {tuple(row["subgroup"]): row for row in self.payload["subgroup_reduction_audit"]}
        self.assertTrue(audit[(1,)]["k5_singleton_invariant"])
        self.assertFalse(audit[(1, 5)]["k5_singleton_invariant"])
        self.assertFalse(audit[(1, 7)]["k5_singleton_invariant"])
        self.assertTrue(audit[(1, 11)]["k5_singleton_invariant"])
        self.assertFalse(audit[(1, 5, 7, 11)]["k5_singleton_invariant"])
        self.assertEqual(audit[(1, 11)]["orbit_size_from_k5_d5"], 1)
        self.assertEqual(audit[(1, 5, 7, 11)]["orbit_size_from_k5_d5"], 2)

    def test_interpretation_and_guardrails(self):
        interpretation = self.payload["selector_interpretation"]
        self.assertIn("one Aut(Z_12) orbit", interpretation["obstruction"])
        self.assertIn("{1,11}", interpretation["conditional_escape"])
        self.assertIn("Anti-Nyquist filtering", interpretation["relation_to_anti_nyquist_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("Full Aut(Z_12) invariance forbids", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
