from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2691_s1641_alpha_geo_role_safe_amplitude_source_audit.py"
OUT = ROOT / "generated" / "p2691_s1641_alpha_geo_role_safe_amplitude_source_audit.json"
MD = ROOT / "generated" / "p2691_s1641_alpha_geo_role_safe_amplitude_source_audit.md"


class P2691AlphaGeoRoleSafeAmplitudeSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_alpha_geo_frontier(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("alpha_geo role-safe", audit["mode"])
        for key in (
            "p2690_selected_p2691",
            "p2680_amplitude_atom",
            "strict_alpha_source",
            "legacy_role_blockers",
            "forbidden_imports",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_separate_strict_value_from_role_safe_source(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2690_selected_p2691"])
        self.assertTrue(state["alpha_strict_source_exported"])
        self.assertEqual(state["alpha_value"], "4 ln(2)")
        self.assertTrue(state["p2680_scalar_shape_normalization_present"])
        self.assertFalse(state["p2680_scalar_shape_role_safe"])
        self.assertFalse(state["p2680_amplitude_role_safe_source_exported"])
        self.assertTrue(state["role_draft_blocks_scalar_normalization_as_physical_role_proof"])
        self.assertTrue(state["toe_audit_keeps_apd_source_open"])

    def test_scalar_normalization_computation_is_exact_but_not_role_semantic(self) -> None:
        comp = self.payload["amplitude_normalization_computation"]
        self.assertTrue(comp["alpha_matches_4_ln2"])
        self.assertTrue(comp["scalar_shape_normalization_is_exact"])
        self.assertLess(comp["max_abs_residual"], 1e-15)
        self.assertEqual(len(comp["sample_rows"]), 13)
        self.assertFalse(comp["role_semantics_generated_by_this_computation"])

    def test_obligation_matrix_keeps_missing_role_safe_requirements(self) -> None:
        matrix = {row["obligation"]: row for row in self.payload["obligation_matrix"]}
        self.assertTrue(matrix["strict_alpha_geo_value_source"]["satisfied_now"])
        self.assertTrue(matrix["exact_scalar_shape_normalization"]["satisfied_now"])
        self.assertFalse(matrix["amplitude_absorption_as_strict_bridge_source"]["satisfied_now"])
        self.assertFalse(matrix["physical_role_safety_theorem"]["satisfied_now"])
        self.assertFalse(matrix["apd_or_lagrangian_dynamical_source"]["satisfied_now"])

    def test_decision_bounds_alpha_atom_and_updates_docs(self) -> None:
        decision = self.payload["decision"]
        self.assertFalse(decision["amplitude_role_safe_source_exported_now"])
        self.assertTrue(decision["bounded_no_go_for_current_alpha_geo_amplitude_atom"])
        self.assertIn("P2692", decision["next_honest_step"])
        self.assertFalse(decision["role_transfer_started_now"])
        self.assertFalse(decision["ltotal_promoted_now"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2691/S1641", MD.read_text(encoding="utf-8"))
        self.assertIn("P2691/S1641", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2691/S1641", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2691/S1641", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
