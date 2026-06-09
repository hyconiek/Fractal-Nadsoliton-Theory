from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2603NadsolitonFractalCodimensionSlopeSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["nadsoliton_fractal_codimension_slope_source_theorem"]["theorem_export"]
        cls.audit = cls.theorem["fractal_codimension_slope_audit"]
        cls.assignment = cls.theorem["strict_damping_source_assignment_after_p2603"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2603")
        self.assertEqual(self.payload["stage_id"], "S1553")
        self.assertIn("SLOPE_SOURCE", self.payload["status"])
        self.assertTrue(self.payload["gatekeeper_checks"]["p2543_slope_obstruction_inherited"])
        self.assertTrue(self.payload["gatekeeper_checks"]["p2602_prime_log_source_inherited"])

    def test_codimension_selects_delta_four_fifths(self) -> None:
        self.assertEqual(self.audit["df_exact"], "9/5")
        self.assertEqual(self.audit["fractal_codimension_formula"], "delta = D_f - 1")
        self.assertEqual(self.audit["codimension_slope_delta"], "4/5")
        self.assertEqual(self.audit["prime_anchor_formula"], "v_p = (D_f - 1) log(p) = (4/5) log(p)")
        self.assertEqual(self.audit["candidate_slope_row_count"], 4)
        self.assertEqual(self.audit["accepting_slope_row_count"], 1)
        self.assertTrue(self.audit["only_codimension_candidate_satisfies_prime_anchor"])
        self.assertTrue(self.audit["all_candidates_prime_log_proportional"])

    def test_strict_damping_source_discharge(self) -> None:
        self.assertTrue(self.theorem["fractal_codimension_slope_source_exported"])
        self.assertTrue(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertTrue(self.theorem["prime_log_proportionality_source_exported"])
        self.assertTrue(self.theorem["multiplicative_character_law_source_exported"])
        self.assertTrue(self.theorem["m2_operator_signature_source_exported"])
        self.assertTrue(self.theorem["beta_eta_numeric_source_exported"])
        self.assertTrue(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertTrue(self.assignment["strict_damping_beta_eta_source_accepts"])
        self.assertEqual(self.assignment["remaining_missing_source_key_count_after_p2603"], 0)

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["legacy_to_strict_completion_bridge_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_theorem"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2603/S1553", MD.read_text(encoding="utf-8"))
        self.assertIn("P2603/S1553", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2603/S1553", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
