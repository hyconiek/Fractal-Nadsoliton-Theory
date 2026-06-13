from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2690_s1640_selector_free_nonexact_boundary_phase_sector_provider_audit.py"
OUT = ROOT / "generated" / "p2690_s1640_selector_free_nonexact_boundary_phase_sector_provider_audit.json"
MD = ROOT / "generated" / "p2690_s1640_selector_free_nonexact_boundary_phase_sector_provider_audit.md"


class P2690SelectorFreeNonexactBoundaryPhaseSectorProviderAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_provider_frontier(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("selector-free nonexact", audit["mode"])
        for key in (
            "p2689_selected_p2690",
            "chain_level_one_bit_carrier",
            "sector_provider_candidates",
            "pair12_subinterface_blockers",
            "forbidden_imports",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_use_p2689_and_prior_provider_audits(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2689_selected_p2690"])
        self.assertTrue(state["p2663_nonexact_sector_needed"])
        self.assertFalse(state["p2663_unique_bits_without_sector_choice"])
        self.assertEqual(state["p2663_exact_coboundary_square_values"], [0])
        self.assertFalse(state["p2664_local_even_class_selects_nonexact"])
        self.assertTrue(state["p2664_declared_theta_can_select"])
        self.assertEqual(state["p2664_passing_sector_selectors"], [])
        self.assertEqual(state["p2665_passing_selector_transfers"], [])
        self.assertFalse(state["p2672_typed_descent_exported"])
        self.assertFalse(state["p2673_pair12_subinterface_exported"])

    def test_carrier_enumeration_separates_exact_from_nonexact_bit(self) -> None:
        carrier = self.payload["carrier_enumeration"]
        self.assertTrue(carrier["carrier_exists"])
        self.assertGreater(carrier["exact_rows"], 0)
        self.assertGreater(carrier["nonexact_rows"], 0)
        self.assertEqual(carrier["exact_square_bit_values"], [0])
        self.assertEqual(carrier["nonexact_square_bit_values"], [1])
        self.assertTrue(carrier["nonzero_bit_requires_nonexact_sector"])

    def test_provider_matrix_has_no_selector_free_exported_pass(self) -> None:
        providers = self.payload["provider_matrix"]
        self.assertEqual(len(providers), 5)
        self.assertTrue(any(row["provider_id"] == "declared_theta_holonomy_source" and row["prefers_nonexact_sector"] and not row["selector_free"] for row in providers))
        self.assertFalse(any(row["selector_free"] and row["prefers_nonexact_sector"] and row["exported_now"] for row in providers))
        self.assertTrue(all(not row["exported_now"] for row in providers))

    def test_decision_freezes_entropy_uv_route_without_exports_and_updates_docs(self) -> None:
        decision = self.payload["decision"]
        self.assertTrue(decision["carrier_exists"])
        self.assertFalse(decision["selector_free_provider_exported"])
        self.assertTrue(decision["freeze_entropy_uv_unit_route"])
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["uv_unit_or_beta_source_exported_now"])
        self.assertFalse(decision["role_transfer_started_now"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertIn("P2691", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2690/S1640", MD.read_text(encoding="utf-8"))
        self.assertIn("P2690/S1640", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2690/S1640", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2690/S1640", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
