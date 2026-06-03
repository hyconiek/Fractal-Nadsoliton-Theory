#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import p2473_s1423_strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit as p2473

OUT = ROOT / "generated" / "p2473_s1423_strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit.json"
MD = ROOT / "generated" / "p2473_s1423_strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit.md"


class P2473StrictPointwiseIntervalDecimalP2459ReplayChainFingerprintBindingAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not OUT.exists():
            p2473.main()
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_replay_chain_fingerprint_binding_audit"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2473")
        self.assertEqual(self.payload["stage_id"], "S1423")
        self.assertIn("P2459_REPLAY_CHAIN_FINGERPRINT_BINDING_AUDIT", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])

    def test_artifact_binding_results(self) -> None:
        self.assertEqual(len(self.theorem["artifact_binding_audits"]), 4)
        self.assertTrue(self.theorem["all_artifact_packet_and_stage_ids_match"])
        self.assertTrue(self.theorem["all_theorem_fingerprints_match_declared"])
        self.assertTrue(self.theorem["all_declared_source_fingerprints_match_current_sources"])
        self.assertTrue(self.theorem["all_audited_gatekeepers_pass"])
        for artifact in self.theorem["artifact_binding_audits"]:
            self.assertTrue(artifact["packet_id_matches_expected"])
            self.assertTrue(artifact["stage_id_matches_expected"])
            self.assertTrue(artifact["theorem_fingerprint_matches_declared"])
            self.assertTrue(artifact["all_source_fingerprints_match_declared"])
            self.assertTrue(artifact["all_gatekeepers_pass"])
            self.assertRegex(artifact["artifact_file_sha256"], r"^[0-9a-f]{64}$")

    def test_chain_consistency_counts(self) -> None:
        chain = self.theorem["chain_consistency"]
        self.assertEqual(chain["p2466_descent_tail_count_from_p2471"], 6361)
        self.assertEqual(chain["p2469_full_opposite_tail_count"], 45165)
        self.assertEqual(chain["p2470_remaining_non_tail_count"], 48320)
        self.assertEqual(chain["p2471_p2459_universe_count"], 99846)
        self.assertEqual(chain["finite_replay_partition_sum"], 99846)
        self.assertTrue(chain["finite_replay_partition_sum_equals_p2459_universe"])
        self.assertEqual(chain["p2470_remaining_budget_after_replay"], 0)
        self.assertEqual(chain["p2471_missing_cells"], 0)
        self.assertEqual(chain["p2471_extra_cells"], 0)
        self.assertTrue(chain["p2471_all_family_partitions_are_disjoint_and_complete"])
        self.assertEqual(chain["p2472_transition_pair_count"], 4)
        self.assertEqual(chain["p2472_seam_replay_count"], 16)
        self.assertTrue(self.theorem["finite_replay_chain_count_binding_passes"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        for name, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, name)
        self.assertTrue(self.theorem["finite_replay_chain_fingerprint_binding_audit_exported"])
        self.assertFalse(self.theorem["directed_rounding_interval_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["symbolic_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["analytic_monotonicity_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["global_continuum_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["legacy_role_transfer_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_rg_audit_docs_and_lay_summary(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        md_text = MD.read_text(encoding="utf-8")
        self.assertIn("Plain-language progress note", md_text)
        self.assertIn("serial numbers and seals", self.theorem["lay_summary"])
        self.assertIn("P2473/S1423", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2473/S1423", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
