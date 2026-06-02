#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2417_s1367_apd_witness_to_source_obligation_nonpromotion_matrix_certificate.py"
OUT = ROOT / "generated" / "p2417_s1367_apd_witness_to_source_obligation_nonpromotion_matrix_certificate.json"
MD = ROOT / "generated" / "p2417_s1367_apd_witness_to_source_obligation_nonpromotion_matrix_certificate.md"
P2416 = ROOT / "generated" / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json"


class P2417APDWitnessToSourceObligationNonpromotionMatrixCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2416.exists():
            subprocess.run([sys.executable, str(ROOT / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["apd_witness_to_source_obligation_nonpromotion_matrix_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_witness_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2417")
        self.assertEqual(self.payload["stage_id"], "S1367")
        self.assertEqual(self.payload["status"], "PASS_VALUE_WITNESSES_POSITIVE_SOURCE_DISCHARGE_ZERO_NO_BRIDGE_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_nonpromotion_matrix(self) -> None:
        self.assertEqual(self.theorem["source_obligation_atom_count"], 8)
        self.assertEqual(self.theorem["artifact_column_count"], 4)
        self.assertTrue(self.theorem["all_recent_witnesses_positive"])
        rows = self.finite["nonpromotion_matrix_rows"]
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(not row["discharged"] for row in rows))
        self.assertTrue(all(not any(row["discharged_by_artifacts"].values()) for row in rows))
        self.assertEqual(self.theorem["discharged_atoms"], [])

    def test_source_lattice_counts(self) -> None:
        self.assertEqual(self.theorem["current_source_discharge_mask"], 0)
        self.assertEqual(self.theorem["full_source_discharge_mask"], 255)
        self.assertEqual(self.theorem["source_obligation_assignment_count"], 256)
        self.assertEqual(self.theorem["proper_subset_failure_count"], 255)
        self.assertEqual(self.theorem["nearest_miss_count"], 8)
        self.assertFalse(self.theorem["bridge_source_ready_from_current_artifacts"])
        self.assertTrue(self.theorem["p2411_full_bridge_rule_inherited"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertTrue(self.theorem["no_p2413_to_p2416_artifact_discharges_source_atom"])
        self.assertTrue(self.theorem["value_witness_to_source_nonpromotion_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertIn("No P2413-P2416 value witness", "\n".join(self.theorem["not_licensed"]))
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("zero source discharge", MD.read_text(encoding="utf-8"))
        self.assertIn("P2417/S1367 APD", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2417/S1367 APD", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
