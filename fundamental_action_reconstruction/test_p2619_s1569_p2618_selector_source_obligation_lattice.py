from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2619_s1569_p2618_selector_source_obligation_lattice import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2619SelectorSourceObligationLatticeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["p2618_selector_source_obligation_lattice"]["theorem_export"]

    def test_identity_and_inherited_blocks(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2619")
        self.assertEqual(self.payload["stage_id"], "S1569")
        self.assertIn("SELECTOR_SOURCE_OBLIGATION_LATTICE", self.payload["status"])
        self.assertTrue(self.theorem["inherits_p2618_full_completion_block"])
        self.assertTrue(self.theorem["inherits_p2616_role_block"])

    def test_c2_enumeration_blocks_invariant_legacy_scalar_sources(self) -> None:
        rows = self.theorem["equivariance_theorem"]["finite_enumeration_rows"]
        invariant = [row for row in rows if row["action_kind"] == "invariant"]
        torsor = [row for row in rows if row["action_kind"] == "torsor"]
        self.assertTrue(invariant)
        self.assertTrue(torsor)
        self.assertTrue(all(row["equivariant_function_count"] == 0 for row in invariant))
        self.assertTrue(all(row["equivariant_function_count"] > 0 for row in torsor))
        proof_text = " ".join(self.theorem["equivariance_theorem"]["proof_steps"])
        self.assertIn("f(x)=-f(x)", proof_text)
        self.assertIn("source torsor", proof_text)

    def test_minimal_lattice_obligations(self) -> None:
        minimal = sorted(self.theorem["premise_lattice"]["minimal_accepting_supports"])
        self.assertEqual(minimal, sorted([["orientation_odd_source"], ["symmetry_breaking_boundary"], ["spin_pin_orientation_source"]]))
        rejected = self.theorem["premise_lattice"]["legacy_scalar_or_axis_only_rejected_supports"]
        self.assertIn(["beta_tors_scalar_invariant"], rejected)
        self.assertIn(["axis_only_selector_up_to_Z2"], rejected)
        self.assertIn(["beta_tors_scalar_invariant", "axis_only_selector_up_to_Z2"], rejected)

    def test_scope_guards(self) -> None:
        self.assertFalse(self.theorem["strict_selector_source_exported"])
        self.assertFalse(self.theorem["beta_tors_chi11_route_reopened"])
        self.assertFalse(self.theorem["gf2_bridge_revalidated"])
        self.assertFalse(self.theorem["role_transfer_revalidated"])
        self.assertFalse(self.theorem["role_bearing_ltotal_reenabled"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2619/S1569", MD.read_text(encoding="utf-8"))
        self.assertIn("P2619/S1569", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2619/S1569", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
