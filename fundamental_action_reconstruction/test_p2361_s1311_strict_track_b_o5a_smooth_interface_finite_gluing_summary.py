from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2361_s1311_strict_track_b_o5a_smooth_interface_finite_gluing_summary.py"
OUT = ROOT / "generated" / "p2361_s1311_strict_track_b_o5a_smooth_interface_finite_gluing_summary.json"
MD = ROOT / "generated" / "p2361_s1311_strict_track_b_o5a_smooth_interface_finite_gluing_summary.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2361O5aSmoothInterfaceFiniteGluingSummaryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_o5a_smooth_interface_finite_gluing_summary_probe"]
        cls.summary = cls.probe["track_B_o5a_smooth_interface_finite_gluing_summary"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2361_s1311_v1")
        self.assertEqual(self.payload["packet_id"], "P2361")
        self.assertEqual(self.payload["stage_id"], "S1311")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_local_o5a_identity_and_scope(self) -> None:
        self.assertEqual(self.summary["local_chern_polynomial"], "8*(K*x1 + K*x2 + K*x3 + 2*x1*x2*x3)")
        self.assertEqual(self.summary["local_reversed_normal_polynomial"], "-8*(K*x1 + K*x2 + K*x3 + 2*x1*x2*x3)")
        self.assertEqual(self.summary["local_o5a_residual"], "0")
        self.assertEqual(
            self.summary["o5a_smooth_interface_status"],
            "FORMAL_SUMMARY_EXPORTED_FOR_FINITE_SMOOTH_CLOSED_INTERFACES_ONLY",
        )
        self.assertTrue(self.summary["o5_full_status"].startswith("OPEN"))
        self.assertIn("not move to corners", self.summary["remaining_gap"])

    def test_finite_incidence_case_residuals(self) -> None:
        rows = {row["case_id"]: row for row in self.summary["case_rows"]}
        self.assertEqual(len(rows), 3)
        self.assertTrue(self.summary["all_case_columns_balanced"])
        self.assertTrue(self.summary["all_case_residuals_zero"])
        self.assertEqual(rows["two_region_single_smooth_closed_interface"]["column_sums"], ["0"])
        self.assertEqual(rows["three_region_smooth_chain_two_interfaces"]["column_sums"], ["0", "0"])
        self.assertEqual(rows["four_region_disjoint_smooth_interface_tree"]["column_sums"], ["0", "0", "0"])
        for row in rows.values():
            self.assertTrue(row["all_columns_balanced"], row["case_id"])
            self.assertEqual(row["pre_gluing_interface_total"], "0", row["case_id"])
            self.assertEqual(row["interface_residual"], "0", row["case_id"])
            self.assertEqual(row["boundary_residual"], "0", row["case_id"])
            self.assertEqual(row["pairing_residual"], "0", row["case_id"])
            self.assertEqual(row["pre_gluing_degree"], row["post_gluing_degree"], row["case_id"])

    def test_upstream_replays(self) -> None:
        self.assertEqual(self.summary["p2360_sign_reversal_residual_replayed"], "0")
        self.assertEqual(self.summary["p2360_interface_pair_sum_replayed"], "0")
        self.assertEqual(self.summary["p2358_symbolic_interface_residual_replayed"], "0")
        self.assertTrue(self.summary["p2358_all_case_residuals_replayed"])
        self.assertEqual(self.summary["p2359_interface_cancellation_replayed"], "0")
        self.assertEqual(self.summary["p2359_gluing_consistency_residual_replayed"], "0")
        self.assertIn("O5_regularization_corners_and_gluing", self.summary["p2353_minimal_cut_replayed"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        for key, value in self.probe["current_export_dependencies"].items():
            self.assertTrue(value, key)
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("O5 full regularization theorem", theorem["not_licensed"])
        self.assertIn("corner contribution theorem", theorem["not_licensed"])
        self.assertIn("legacy-to-strict kernel bridge", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
