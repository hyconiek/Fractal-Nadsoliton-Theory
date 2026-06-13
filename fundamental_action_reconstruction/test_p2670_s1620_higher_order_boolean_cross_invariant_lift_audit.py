from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2670_s1620_higher_order_boolean_cross_invariant_lift_audit.py"
OUT = ROOT / "generated" / "p2670_s1620_higher_order_boolean_cross_invariant_lift_audit.json"
MD = ROOT / "generated" / "p2670_s1620_higher_order_boolean_cross_invariant_lift_audit.md"
class P2670HigherOrderBooleanCrossInvariantLiftAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload=json.loads(OUT.read_text(encoding="utf-8"))
    def test_content_first_grep_covers_research_not_numbers(self):
        audit=self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in ("higher_order_boolean_content","physical_origin_content","boundary_cycle_sector_content","selector_content","nonclosure_guard_content"):
            self.assertIn(key, audit["patterns"]); self.assertGreater(audit["patterns"][key]["count"],0)
    def test_upstream_p2669_consistency(self):
        u=self.payload["upstream_consistency"]
        self.assertTrue(u["p2669_mathematical_candidates_exist"])
        self.assertTrue(u["p2669_no_mixed_candidate"])
        self.assertTrue(u["p2669_no_passing_source"])
        self.assertTrue(u["p2669_no_convention_free_source"])
    def test_higher_order_enumeration_no_source(self):
        w=self.payload["higher_order_witness"]
        self.assertEqual(w["total_functions"],256)
        self.assertGreater(w["candidate_count"],0)
        self.assertGreater(w["higher_order_candidate_count"],0)
        self.assertGreater(w["auxiliary_dependent_candidate_count"],0)
        self.assertEqual(w["passing_source_count"],0)
        self.assertFalse(w["convention_free_physical_origin_exported_for_any_candidate"])
    def test_no_closure_and_docs_updated(self):
        d=self.payload["closure_decision"]
        self.assertTrue(d["mathematical_higher_order_candidates_exist"])
        self.assertEqual(d["passing_source_count"],0)
        self.assertFalse(d["boundary_phase_bit_target_exported_now"])
        self.assertFalse(d["beta_source_exported_now"])
        self.assertFalse(d["qw2191_discharged_now"])
        self.assertFalse(d["role_bearing_ltotal_now"])
        self.assertFalse(d["toe_closure_now"])
        self.assertTrue(all(v is False for v in self.payload["negative_export_flags"].values()))
        self.assertIn("P2670/S1620", MD.read_text(encoding="utf-8"))
        self.assertIn("P2670/S1620", (ROOT/"STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2670/S1620", (ROOT/"STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
if __name__ == "__main__": unittest.main()
