import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_recursive_self_information_mirror_selector_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_recursive_self_information_mirror_selector_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_recursive_self_information_mirror_selector_audit_report.md"


class StrictAlphaRecursiveSelfInformationMirrorSelectorAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_RECURSIVE_SELF_INFORMATION_MIRROR_SELECTOR_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "recursive-self-information-mirror-depth-does-not-create-unit-orientation-bit",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_recursive_level_audit(self):
        audit = self.payload["recursive_self_information_level_audit"]
        self.assertEqual(audit["max_level"], 5)
        self.assertEqual(audit["enumerated_max_level"], 3)
        self.assertEqual(audit["singleton_d5_fixed_count_across_levels"], 0)
        self.assertFalse(audit["aut_invariant_singleton_d5_selector_exists_any_level"])
        rows = audit["rows"]
        self.assertEqual([row["level"] for row in rows], list(range(6)))
        self.assertEqual([row["object_count"] for row in rows[:4]], [2, 4, 16, 65536])
        self.assertEqual([row["fixed_record_count"] for row in rows[:4]], [0, 2, 8, 4096])
        self.assertEqual([row["mirror_orbit_count"] for row in rows[:4]], [1, 3, 12, 34816])
        self.assertTrue(all(not row["singleton_d5_tower_fixed_by_mirror"] for row in rows))
        self.assertFalse(rows[4]["enumerated_exactly"])
        self.assertEqual(rows[4]["object_count"], "2^65536")

    def test_reflection_depth_and_proof_certificate(self):
        depth = self.payload["reflection_depth_audit"]
        self.assertFalse(depth["depth_as_scalar_can_break_unit_mirror"])
        self.assertTrue(all(row["depth_is_aut_invariant_scalar"] for row in depth["rows"]))
        self.assertTrue(all(not row["depth_selects_A5_over_A1"] for row in depth["rows"]))

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("swaps A1 and A5", proof["base_action"])
        self.assertIn("powerset", proof["recursive_lift"])
        self.assertIn("whole J_n-orbits", proof["fixed_record_condition"])
        self.assertIn("never fixed", proof["singleton_tower_obstruction"])
        self.assertIn("Aut-invariant scalar", proof["information_reading"])

    def test_interpretation_and_guardrails(self):
        interpretation = self.payload["interpretation"]
        self.assertIn("recursion alone preserves", interpretation["direct_correction"])
        self.assertIn("information-about-information", interpretation["what_it_can_do"])
        self.assertIn("Aut-breaking", interpretation["what_it_cannot_do"])
        self.assertIn("orientation character", interpretation["information_principle"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives recursive self-information", hard_limits)
        self.assertIn("does not discharge QW-2191", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
