#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2431_s1381_admissible_next_theorem_target_antichain_certificate.py"
OUT = ROOT / "generated" / "p2431_s1381_admissible_next_theorem_target_antichain_certificate.json"
MD = ROOT / "generated" / "p2431_s1381_admissible_next_theorem_target_antichain_certificate.md"
P2430 = ROOT / "generated" / "p2430_s1380_repair_derivative_witness_cover_minimality_certificate.json"


class P2431AdmissibleNextTheoremTargetAntichainCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2430.exists():
            subprocess.run([sys.executable, str(ROOT / "p2430_s1380_repair_derivative_witness_cover_minimality_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["admissible_next_theorem_target_antichain_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2431")
        self.assertEqual(self.payload["stage_id"], "S1381")
        self.assertEqual(self.payload["status"], "PASS_ADMISSIBLE_NEXT_TARGET_ANTICHAIN_NO_GATE_DISCHARGE_NO_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_candidate_counts_and_singletons(self) -> None:
        self.assertEqual(self.theorem["candidate_row_count_size_le_2"], 15)
        self.assertEqual(self.theorem["admissible_candidate_count_size_le_2"], 6)
        self.assertEqual(self.theorem["inadmissible_candidate_count_size_le_2"], 9)
        self.assertEqual(
            self.theorem["admissible_singletons_from_current"],
            [["source_obligation_discharge"], ["chi11_source_export"], ["qw2191_selector_discharge"]],
        )
        self.assertEqual(
            self.theorem["inadmissible_singletons_from_current"],
            [["role_transfer_audit_license"], ["role_bearing_ltotal_export"]],
        )

    def test_pareto_antichain(self) -> None:
        self.assertEqual(
            self.theorem["minimal_readiness_complete_candidates"],
            [["source_obligation_discharge"], ["chi11_source_export", "qw2191_selector_discharge"]],
        )
        self.assertEqual(
            self.theorem["max_completed_target_coverage_candidates"],
            [["source_obligation_discharge"], ["chi11_source_export", "qw2191_selector_discharge"]],
        )
        self.assertFalse(self.theorem["role_transfer_admissible_as_first_move"])
        self.assertFalse(self.theorem["role_bearing_ltotal_admissible_as_first_move"])

    def test_inherited_hard_limits_and_docs(self) -> None:
        self.assertTrue(self.theorem["p2430_global_minimal_cover_inherited"])
        self.assertTrue(self.theorem["p2430_selector_pair_inherited"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("admissible next theorem-target", MD.read_text(encoding="utf-8"))
        self.assertIn("P2431/S1381", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2431/S1381", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
