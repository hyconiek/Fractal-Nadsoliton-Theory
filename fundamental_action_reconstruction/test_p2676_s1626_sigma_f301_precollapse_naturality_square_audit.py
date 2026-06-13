from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2676_s1626_sigma_f301_precollapse_naturality_square_audit.py"
OUT = ROOT / "generated" / "p2676_s1626_sigma_f301_precollapse_naturality_square_audit.json"
MD = ROOT / "generated" / "p2676_s1626_sigma_f301_precollapse_naturality_square_audit.md"


class P2676SigmaF301PrecollapseNaturalitySquareAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_naturality_route(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "naturality_square_content",
            "sigma_to_f301_content",
            "precollapse_content",
            "chart_label_retaining_content",
            "nonprojector_content",
            "nonconvention_content",
            "selector_source_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2675_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2675_s6_not_passed"])
        self.assertTrue(upstream["p2675_o3_not_passed"])
        self.assertTrue(upstream["p2675_missing_pre_collapse"])
        self.assertTrue(upstream["p2675_missing_nonprojector"])
        self.assertTrue(upstream["p2675_missing_no_fiat"])

    def test_finite_boolean_naturality_witness(self) -> None:
        witness = self.payload["finite_boolean_naturality_witness"]
        self.assertEqual(witness["total_maps_checked"], 16)
        self.assertEqual(witness["formal_square_candidate_count"], 2)
        self.assertEqual(witness["passing_export_gate_count"], 0)
        self.assertEqual(
            set(witness["formal_square_convention_classes"]),
            {"xor_orientation", "xnor_reversal_orientation"},
        )
        formal_rows = [row for row in witness["rows"] if row["passes_formal_square"]]
        self.assertTrue(all(row["source_sensitive"] for row in formal_rows))
        self.assertTrue(all(row["chart_sensitive_nonprojector_form"] for row in formal_rows))
        self.assertFalse(any(row["has_internal_orientation_source"] for row in formal_rows))

    def test_obstruction_certificates_block_export(self) -> None:
        certs = {row["certificate"]: row["value"] for row in self.payload["obstruction_certificates"]}
        self.assertTrue(certs["finite_formal_candidates_exist"])
        self.assertTrue(certs["orientation_reversal_pair_survives"])
        self.assertTrue(certs["internal_orientation_source_absent"])
        self.assertTrue(certs["export_gate_zero"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["formal_square_candidate_count"], 2)
        self.assertEqual(decision["passing_export_gate_count"], 0)
        self.assertFalse(decision["precollapse_naturality_square_exported_now"])
        self.assertFalse(decision["sigma_to_f301_typed_arrow_exported_now"])
        self.assertFalse(decision["s6_exported_now"])
        self.assertFalse(decision["o3_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2676/S1626", MD.read_text(encoding="utf-8"))
        self.assertIn("P2676/S1626", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2676/S1626", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
