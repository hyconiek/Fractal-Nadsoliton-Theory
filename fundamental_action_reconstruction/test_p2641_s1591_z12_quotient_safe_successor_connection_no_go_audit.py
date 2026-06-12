from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2641_s1591_z12_quotient_safe_successor_connection_no_go_audit.py"
OUT = ROOT / "generated" / "p2641_s1591_z12_quotient_safe_successor_connection_no_go_audit.json"
MD = ROOT / "generated" / "p2641_s1591_z12_quotient_safe_successor_connection_no_go_audit.md"


class P2641Z12QuotientSafeSuccessorConnectionNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("z12_successor_connection_content", audit["patterns"])
        self.assertIn("offset_stride_origin_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_z12_group_invariants_are_computed(self) -> None:
        inv = self.payload["z12_group_invariants"]
        self.assertEqual(inv["aut_group"], [1, 5, 7, 11])
        self.assertEqual(inv["fixed_elements_usable_as_aut_invariant_origins"], [0, 6])
        self.assertEqual(inv["nonzero_aut_invariant_translation_strides"], [6])
        self.assertEqual(inv["fixed_generators"], [])

    def test_strict_aut_invariant_scan_does_not_preserve_inverse_hierarchy(self) -> None:
        scan = self.payload["strict_aut_invariant_lift_scan"]
        self.assertEqual(scan["candidate_count"], 2)
        self.assertEqual(scan["uv_safe_count"], 1)
        self.assertEqual(scan["inverse_hierarchy_count"], 0)
        uv_safe = [row for row in scan["candidates"] if row["uv_positive_at_d1"]]
        self.assertEqual(uv_safe[0]["k0"], 6)
        self.assertLess(uv_safe[0]["strict_lifted_abs_k7_over_k1"], 1.0)

    def test_p2639_role_like_lifts_remain_not_strict_aut_sources(self) -> None:
        rows = self.payload["role_like_compatibility"]
        pairs = {(row["k0"], row["stride"]) for row in rows}
        self.assertIn((4, 3), pairs)
        self.assertIn((10, 6), pairs)
        self.assertFalse(any(row["strict_aut_invariant_origin_and_stride"] for row in rows))
        self.assertTrue(any(row["stride_aut_invariant"] and not row["k0_aut_fixed"] for row in rows))

    def test_premise_based_fixing_is_support_but_not_role_like_source(self) -> None:
        premise = self.payload["premise_based_fixing_audit"]
        self.assertTrue(premise["premise_based_fixing_is_real_support"])
        self.assertFalse(premise["but_exports_zero_lattice_origin"])
        self.assertFalse(premise["but_exports_role_like_stride_3_or_6"])
        self.assertFalse(premise["minimal_identity_origin_stride1_lift"]["uv_positive_at_d1"])
        self.assertLess(premise["nearest_uv_safe_stride1_lift_from_p2639_family"]["strict_lifted_abs_k7_over_k1"], 1.0)

    def test_closure_decision_and_docs_do_not_promote(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["promote_successor_connection_to_bridge_completion"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertIn("FINITE_Z12_SUCCESSOR_CONNECTION_NO_GO", decision["classification"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2641/S1591", MD.read_text(encoding="utf-8"))
        self.assertIn("P2641/S1591", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2641/S1591", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
