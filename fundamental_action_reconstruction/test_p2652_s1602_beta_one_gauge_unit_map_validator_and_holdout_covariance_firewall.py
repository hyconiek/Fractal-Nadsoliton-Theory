from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2652_s1602_beta_one_gauge_unit_map_validator_and_holdout_covariance_firewall.py"
OUT = ROOT / "generated" / "p2652_s1602_beta_one_gauge_unit_map_validator_and_holdout_covariance_firewall.json"
MD = ROOT / "generated" / "p2652_s1602_beta_one_gauge_unit_map_validator_and_holdout_covariance_firewall.md"


class P2652BetaOneGaugeUnitMapValidatorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("gauge_unit_map_content", audit["patterns"])
        self.assertIn("holdout_covariance_content", audit["patterns"])
        self.assertIn("fake_unit_pass_firewall_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2647_schema_ready_but_no_real_holdout"])
        self.assertTrue(upstream["p2648_margin_ready_but_no_real_holdout"])
        self.assertTrue(upstream["p2651_allows_beta_one_only_as_working_gauge"])
        self.assertTrue(upstream["p2651_blocks_beta_source"])
        self.assertTrue(upstream["p2651_requires_unit_map_bookkeeping"])

    def test_required_schema_extends_holdout_with_unit_map(self) -> None:
        schema = self.payload["required_gauge_unit_payload_schema"]
        self.assertIn("gauge_scale_a", schema["gauge_contract_keys"])
        self.assertIn("unit_map_source_precommitted_before_holdout", schema["gauge_contract_keys"])
        self.assertIn("raw_near", schema["unit_aware_measurement_keys"])
        self.assertIn("near_in_beta_one_gauge", schema["unit_aware_measurement_keys"])
        self.assertIn("holdout_fit", schema["forbidden_unit_map_sources"])
        self.assertIn("per_pair_fit", schema["forbidden_unit_map_sources"])

    def test_covariance_firewall_blocks_fake_unit_passes(self) -> None:
        firewall = self.payload["gauge_unit_map_firewall"]
        covariant = firewall["covariant_synthetic_fixture"]
        raw_fake = firewall["raw_distance_beta_one_substitution_fixture"]
        target_fake = firewall["target_fit_unit_map_fixture"]
        self.assertTrue(covariant["unit_map_admissible_for_p2647_p2648_execution"])
        self.assertTrue(covariant["scale_equation_ok"])
        self.assertFalse(covariant["covariance_errors"])
        self.assertFalse(raw_fake["unit_map_admissible_for_p2647_p2648_execution"])
        self.assertTrue(raw_fake["covariance_errors"])
        self.assertFalse(raw_fake["raw_distance_guard_declared"])
        self.assertFalse(target_fake["unit_map_admissible_for_p2647_p2648_execution"])
        self.assertTrue(target_fake["unit_map_source_forbidden_as_target_fit"])
        self.assertFalse(target_fake["unit_map_source_precommitted_before_holdout"])
        self.assertTrue(firewall["firewall_passes"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertTrue(decision["unit_map_firewall_passes"])
        self.assertFalse(decision["real_blind_holdout_payload_loaded"])
        self.assertFalse(decision["unit_map_source_theorem_exported_now"])
        self.assertFalse(decision["empirical_confirmation_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2652/S1602", MD.read_text(encoding="utf-8"))
        self.assertIn("P2652/S1602", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2652/S1602", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
