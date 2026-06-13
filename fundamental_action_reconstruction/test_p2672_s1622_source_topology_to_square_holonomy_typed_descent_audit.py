from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2672_s1622_source_topology_to_square_holonomy_typed_descent_audit.py"
OUT = ROOT / "generated" / "p2672_s1622_source_topology_to_square_holonomy_typed_descent_audit.json"
MD = ROOT / "generated" / "p2672_s1622_source_topology_to_square_holonomy_typed_descent_audit.md"


class P2672SourceTopologyToSquareHolonomyTypedDescentAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_descent_not_numbers(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "source_topology_sign_content",
            "square_holonomy_content",
            "typed_descent_content",
            "sign_preservation_content",
            "reversal_guard_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2671_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2671_has_near_miss_pair_sector_subsets"])
        self.assertTrue(upstream["p2671_no_passing_observable_origin"])
        self.assertTrue(upstream["p2671_no_boundary_phase_bit_target"])
        self.assertTrue(upstream["p2671_no_ltotal_reopening"])

    def test_typed_descent_gate_is_real_near_miss_not_source(self) -> None:
        gate = self.payload["typed_descent_gate"]["criteria"]
        self.assertTrue(gate["source_topology_sign_material_present"])
        self.assertTrue(gate["boundary_square_holonomy_material_present"])
        self.assertTrue(gate["typed_descent_material_present"])
        self.assertTrue(gate["sector_swap_reversal_guard_material_present"])
        self.assertFalse(gate["current_strict_typed_descent_exported"])
        self.assertFalse(gate["current_sign_preservation_proof_exported"])
        self.assertFalse(gate["current_sector_swap_reversal_forbidden"])
        self.assertFalse(gate["bridge_completed_boundary_dynamics_provenance"])
        self.assertFalse(gate["passes_typed_descent_gate_now"])

    def test_finite_map_witness_blocks_convention_reversal(self) -> None:
        witness = self.payload["finite_map_witness"]
        self.assertEqual(witness["total_maps_checked"], 4)
        self.assertEqual(witness["maps_with_positive_to_sector1"], 2)
        self.assertEqual(witness["passing_sourced_descent_count"], 0)
        self.assertTrue(all(row["convention_reversal_available"] for row in witness["rows"]))
        self.assertFalse(any(row["passes_as_sourced_descent_now"] for row in witness["rows"]))

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["passes_typed_descent_gate_now"])
        self.assertEqual(decision["passing_sourced_descent_count"], 0)
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2672/S1622", MD.read_text(encoding="utf-8"))
        self.assertIn("P2672/S1622", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2672/S1622", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
