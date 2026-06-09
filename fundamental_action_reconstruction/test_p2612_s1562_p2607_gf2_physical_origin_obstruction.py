from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2612_s1562_p2607_gf2_physical_origin_obstruction import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2612P2607Gf2PhysicalOriginObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["p2607_gf2_physical_origin_obstruction"]["theorem_export"]
        cls.obs = cls.theorem["obstruction_certificate"]

    def test_identity_and_recomputed_rank(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2612")
        self.assertEqual(self.payload["stage_id"], "S1562")
        self.assertIn("GF2_PHYSICAL_ORIGIN_OBSTRUCTION", self.payload["status"])
        self.assertTrue(self.theorem["gf2_physical_origin_obstruction_exported"])
        self.assertEqual(self.obs["rank_over_gf2_recomputed"], 11)
        self.assertEqual(self.obs["nullity_over_gf2_recomputed"], 0)

    def test_rank_is_algorithmic_not_physical(self) -> None:
        self.assertTrue(self.obs["lower_triangular_unit_matrix"])
        self.assertTrue(self.obs["rank_explained_by_unit_diagonal"])
        self.assertTrue(self.obs["source_origin_audit"]["all_algorithmic_index_markers_present"])
        self.assertTrue(self.obs["rhs_rank_independence_audit"]["rank_independent_of_rhs_samples"])
        self.assertFalse(self.obs["physical_origin_requirement_audit"]["all_physical_requirements_present"])
        self.assertIn("chiral_current_conservation_equations", self.obs["physical_origin_requirement_audit"]["missing_requirements"])

    def test_quarantine_and_alternative(self) -> None:
        self.assertTrue(self.obs["p2610_quarantine_inherited"]["p2607_quarantined"])
        self.assertTrue(self.obs["p2610_quarantine_inherited"]["p2608_quarantined"])
        self.assertTrue(self.obs["p2611_role_semantics_inherited"]["role_semantics_defined"])
        self.assertFalse(self.obs["p2611_role_semantics_inherited"]["current_ltotal_assignment_accepts"])
        self.assertEqual(self.obs["strict_alternative_proposal"]["preferred_path"], "P2601_monoid_action_uniqueness")

    def test_no_reexports(self) -> None:
        self.assertFalse(self.theorem["p2607_quarantine_lifted"])
        self.assertFalse(self.theorem["gf2_physical_origin_theorem_exported"])
        self.assertFalse(self.theorem["bridge_completion_revalidated"])
        self.assertFalse(self.theorem["p2608_role_bearing_ltotal_reenabled"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2612/S1562", MD.read_text(encoding="utf-8"))
        self.assertIn("P2612/S1562", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2612/S1562", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
