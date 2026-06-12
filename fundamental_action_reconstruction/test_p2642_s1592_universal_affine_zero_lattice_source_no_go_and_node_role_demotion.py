from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2642_s1592_universal_affine_zero_lattice_source_no_go_and_node_role_demotion.py"
OUT = ROOT / "generated" / "p2642_s1592_universal_affine_zero_lattice_source_no_go_and_node_role_demotion.json"
MD = ROOT / "generated" / "p2642_s1592_universal_affine_zero_lattice_source_no_go_and_node_role_demotion.md"


class P2642UniversalAffineZeroLatticeSourceNoGoAndNodeRoleDemotionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("universal_affine_zero_lattice_content", audit["patterns"])
        self.assertIn("source_atom_origin_stride_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_universal_affine_parameterization_is_recorded(self) -> None:
        theorem = self.payload["symbolic_parameterization_theorem"]
        self.assertEqual(theorem["legacy_nodes"], [2, 5, 8, 11])
        self.assertIn("r_{k0,s}(d) = (4s/3)d", theorem["all_exact_affine_lifts"])
        self.assertEqual(theorem["uv_positive_at_d1_iff"], "r(1)>0 iff 1 + 3k0 - s > 0")

    def test_strict_aut_source_scan_has_no_inverse_hierarchy_candidate(self) -> None:
        scan = self.payload["strict_aut_source_scan"]
        self.assertEqual(scan["fixed_origins"], [0, 6])
        self.assertEqual(scan["fixed_nonzero_strides"], [6])
        self.assertEqual(scan["candidate_count"], 2)
        self.assertEqual(scan["uv_safe_count"], 1)
        self.assertEqual(scan["inverse_hierarchy_count"], 0)

    def test_premise_generator_with_canonical_origin_still_has_no_inverse_hierarchy_candidate(self) -> None:
        scan = self.payload["premise_generator_with_canonical_origin_scan"]
        self.assertEqual(scan["fixed_origins"], [0, 6])
        self.assertEqual(scan["candidate_count"], 24)
        self.assertEqual(scan["uv_safe_count"], 12)
        self.assertEqual(scan["inverse_hierarchy_count"], 0)
        self.assertEqual(scan["by_origin"]["0"]["uv_safe_count"], 0)
        self.assertEqual(scan["by_origin"]["6"]["uv_safe_count"], 12)

    def test_p2639_role_like_lifts_fail_source_recheck(self) -> None:
        recheck = self.payload["p2639_role_like_source_recheck"]
        self.assertGreaterEqual(recheck["role_like_candidate_count"], 2)
        self.assertEqual(recheck["strict_aut_source_pass_count"], 0)
        self.assertEqual(recheck["premise_generator_plus_canonical_origin_pass_count"], 0)
        self.assertTrue(all(row["source_failure"] == "origin_not_aut_fixed_hidden_selector_required" for row in recheck["rows"]))

    def test_demotion_decision_blocks_closure_and_updates_docs(self) -> None:
        decision = self.payload["demotion_decision"]
        self.assertTrue(decision["gates"]["universal_affine_family_proven"])
        self.assertFalse(decision["gates"]["strict_aut_source_candidate_with_uv_and_inverse_hierarchy_exists"])
        self.assertFalse(decision["gates"]["premise_generator_canonical_origin_candidate_with_uv_and_inverse_hierarchy_exists"])
        self.assertIn("DEMOTE", decision["legacy_integer_node_gauge_role_status"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2642/S1592", MD.read_text(encoding="utf-8"))
        self.assertIn("P2642/S1592", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2642/S1592", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
