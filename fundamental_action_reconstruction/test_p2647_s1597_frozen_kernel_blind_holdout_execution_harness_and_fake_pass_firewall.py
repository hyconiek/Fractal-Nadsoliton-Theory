from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2647_s1597_frozen_kernel_blind_holdout_execution_harness_and_fake_pass_firewall.py"
OUT = ROOT / "generated" / "p2647_s1597_frozen_kernel_blind_holdout_execution_harness_and_fake_pass_firewall.json"
MD = ROOT / "generated" / "p2647_s1597_frozen_kernel_blind_holdout_execution_harness_and_fake_pass_firewall.md"


class P2647FrozenKernelBlindHoldoutHarnessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_schema_are_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("blind_holdout_execution_content", audit["patterns"])
        self.assertIn("tail_ratio_measurement_schema_content", audit["patterns"])
        self.assertIn("fake_pass_or_control_firewall_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        schema = self.payload["required_blind_payload_schema"]
        self.assertIn("p2646_payload_sha256", schema["top_level_keys"])
        self.assertIn("measured_tail_ratio", schema["measurement_keys"])
        self.assertIn("uses_holdout_refit", schema["control_baseline_keys"])

    def test_fake_pass_firewall_accepts_strict_and_rejects_legacy_midpoint(self) -> None:
        firewall = self.payload["fake_pass_firewall"]
        self.assertTrue(firewall["firewall_passes"])
        self.assertTrue(firewall["strict_fixture_expected_pass"]["all_pairs_pass"])
        self.assertFalse(firewall["legacy_fixture_expected_fail"]["all_pairs_pass"])
        self.assertFalse(firewall["midpoint_fixture_expected_fail_by_strict_inequality"]["all_pairs_pass"])
        self.assertGreater(firewall["strict_fixture_expected_pass"]["minimum_ratio_margin"], 0.0)
        self.assertGreater(firewall["strict_fixture_expected_pass"]["minimum_slope_margin"], 0.0)
        self.assertLess(firewall["legacy_fixture_expected_fail"]["minimum_ratio_margin"], 0.0)

    def test_empty_payload_is_not_admissible_and_no_data_claim_is_made(self) -> None:
        witness = self.payload["empty_payload_validation_witness"]
        self.assertFalse(witness["admissible_for_execution"])
        self.assertFalse(witness["schema_ready"])
        self.assertIn("measurements", witness["missing_top_level_keys"])
        data_audit = self.payload["blind_data_presence_audit"]
        self.assertFalse(data_audit["real_blind_holdout_payload_found"])
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["real_blind_holdout_executed"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_docs_are_updated(self) -> None:
        self.assertIn("P2647/S1597", MD.read_text(encoding="utf-8"))
        self.assertIn("P2647/S1597", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2647/S1597", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
