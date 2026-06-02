#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2430_s1380_repair_derivative_witness_cover_minimality_certificate.py"
OUT = ROOT / "generated" / "p2430_s1380_repair_derivative_witness_cover_minimality_certificate.json"
MD = ROOT / "generated" / "p2430_s1380_repair_derivative_witness_cover_minimality_certificate.md"
P2429 = ROOT / "generated" / "p2429_s1379_repair_derivative_nearest_miss_witness_table_certificate.json"


class P2430RepairDerivativeWitnessCoverMinimalityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2429.exists():
            subprocess.run([sys.executable, str(ROOT / "p2429_s1379_repair_derivative_nearest_miss_witness_table_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["repair_derivative_witness_cover_minimality_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2430")
        self.assertEqual(self.payload["stage_id"], "S1380")
        self.assertEqual(
            self.payload["status"],
            "PASS_REPAIR_DERIVATIVE_WITNESS_COVER_MINIMALITY_NO_GATE_DISCHARGE_NO_CLOSURE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_minimal_covers(self) -> None:
        gates = [
            "source_obligation_discharge",
            "chi11_source_export",
            "qw2191_selector_discharge",
            "role_transfer_audit_license",
            "role_bearing_ltotal_export",
        ]
        self.assertEqual(self.theorem["global_minimal_covers"], [gates])
        self.assertTrue(self.theorem["unique_global_minimal_cover_is_all_five_gates"])
        self.assertTrue(self.theorem["unique_toe_minimal_cover_is_all_five_gates"])
        self.assertTrue(self.theorem["selector_minimal_cover_is_chi11_qw2191_pair"])
        self.assertEqual(
            self.theorem["minimal_cover_sizes_by_target"],
            {
                "bridge_source_ready": 1,
                "role_bearing_ltotal_ready": 1,
                "role_transfer_ready": 1,
                "selector_source_ready": 2,
                "toe_ready": 5,
            },
        )

    def test_cover_lattice_counts(self) -> None:
        self.assertEqual(self.theorem["global_covering_row_count"], 1)
        self.assertEqual(self.theorem["global_proper_failure_row_count"], 31)
        self.assertEqual(self.theorem["toe_covering_row_count"], 1)
        self.assertEqual(self.theorem["toe_proper_failure_row_count"], 31)
        self.assertEqual(
            self.theorem["global_uncovered_count_distribution"],
            {"0": 1, "1": 5, "2": 10, "3": 10, "4": 5, "5": 1},
        )

    def test_inherited_hard_limits_and_docs(self) -> None:
        self.assertTrue(self.theorem["p2429_witness_row_count_inherited"])
        self.assertTrue(self.theorem["p2429_target_gate_pair_count_inherited"])
        self.assertTrue(self.theorem["p2429_toe_nearest_miss_count_inherited"])
        self.assertEqual(self.theorem["current_non_theorem_evidence_covers_required_witness_gates"], [])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("witness-cover minimality", MD.read_text(encoding="utf-8"))
        self.assertIn("P2430/S1380", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2430/S1380", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
