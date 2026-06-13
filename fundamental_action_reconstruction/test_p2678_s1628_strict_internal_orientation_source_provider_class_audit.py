from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2678_s1628_strict_internal_orientation_source_provider_class_audit.py"
OUT = ROOT / "generated" / "p2678_s1628_strict_internal_orientation_source_provider_class_audit.json"
MD = ROOT / "generated" / "p2678_s1628_strict_internal_orientation_source_provider_class_audit.md"


class P2678StrictInternalOrientationSourceProviderClassAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_different_provider_class(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "different_provider_class_content",
            "oriented_torsor_content",
            "legacy_scalar_even_content",
            "symmetry_breaking_content",
            "sigma_f301_binding_content",
            "forbidden_collapse_content",
            "closure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2677_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2677_current_route_no_go"])
        self.assertTrue(upstream["p2677_not_global_no_go"])
        self.assertTrue(upstream["p2677_missing_internal_source"])
        self.assertTrue(upstream["p2677_missing_nonfiat_choice"])

    def test_c2_equivariant_provider_enumeration(self) -> None:
        enum = self.payload["c2_equivariant_provider_enumeration"]
        self.assertEqual(enum["sign_maps_checked"], 4)
        self.assertEqual(enum["orientation_odd_map_count"], 2)
        self.assertEqual(enum["orientation_even_map_count"], 2)
        self.assertEqual(
            set(enum["formal_provider_classes_with_odd_domain"]),
            {"oriented_c2_torsor", "spin_pin_orientation_source", "boundary_symmetry_breaking_source"},
        )
        self.assertEqual(enum["passing_provider_class_count"], 0)
        formal_rows = [row for row in enum["provider_rows"] if row["has_formal_orientation_odd_domain"]]
        self.assertTrue(all(row["formal_odd_map_tables"] for row in formal_rows))
        self.assertFalse(any(row["provider_exported_now"] for row in enum["provider_rows"]))

    def test_provider_obligation_lattice(self) -> None:
        lattice = self.payload["provider_obligation_lattice"]
        self.assertEqual(lattice["total_states"], 128)
        self.assertEqual(lattice["passing_states"], 1)
        self.assertEqual(lattice["hamming_distance_to_pass"], 5)
        self.assertTrue(lattice["current_state"]["different_provider_class_selected"])
        self.assertTrue(lattice["current_state"]["closure_guards_preserved"])
        self.assertIn("orientation_odd_domain_exported", lattice["missing_current_obligations"])
        self.assertIn("precollapse_sigma_side_binding_exported", lattice["missing_current_obligations"])

    def test_no_reopening_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["formal_odd_provider_class_count"], 3)
        self.assertEqual(decision["passing_provider_class_count"], 0)
        self.assertFalse(decision["strict_internal_orientation_source_exported_now"])
        self.assertFalse(decision["precollapse_sigma_f301_binding_exported_now"])
        self.assertFalse(decision["s6_reopened_now"])
        self.assertFalse(decision["o3_reopened_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2678/S1628", MD.read_text(encoding="utf-8"))
        self.assertIn("P2678/S1628", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2678/S1628", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
