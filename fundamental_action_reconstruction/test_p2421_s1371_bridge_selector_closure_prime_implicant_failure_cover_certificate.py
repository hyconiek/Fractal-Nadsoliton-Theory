#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2421_s1371_bridge_selector_closure_prime_implicant_failure_cover_certificate.py"
OUT = ROOT / "generated" / "p2421_s1371_bridge_selector_closure_prime_implicant_failure_cover_certificate.json"
MD = ROOT / "generated" / "p2421_s1371_bridge_selector_closure_prime_implicant_failure_cover_certificate.md"
P2420 = ROOT / "generated" / "p2420_s1370_bridge_selector_nonclosure_reason_matrix_certificate.json"


class P2421BridgeSelectorClosurePrimeImplicantFailureCoverCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2420.exists():
            subprocess.run([sys.executable, str(ROOT / "p2420_s1370_bridge_selector_nonclosure_reason_matrix_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["bridge_selector_closure_prime_implicant_failure_cover_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_witness_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2421")
        self.assertEqual(self.payload["stage_id"], "S1371")
        self.assertEqual(self.payload["status"], "PASS_UNIQUE_PRIME_IMPLICANT_FAILURE_COVER_NO_GATE_DISCHARGE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_truth_table_and_anf(self) -> None:
        self.assertEqual(self.theorem["closure_gate_count"], 7)
        self.assertEqual(self.theorem["total_assignment_count"], 128)
        self.assertEqual(self.theorem["truth_vector_true_count"], 1)
        self.assertEqual(self.theorem["full_mask"], 127)
        self.assertEqual(self.theorem["anf_term_count"], 1)
        self.assertEqual(self.theorem["anf_degree"], 7)
        self.assertEqual(self.finite["anf_terms"][0]["mask"], 127)

    def test_prime_implicant_and_failure_cover(self) -> None:
        self.assertEqual(self.theorem["prime_implicant_count"], 1)
        self.assertEqual(self.theorem["prime_implicant_masks"], [127])
        self.assertEqual(self.theorem["success_cnf_unit_clause_count"], 7)
        self.assertEqual(self.theorem["failure_cover_literal_count"], 7)
        self.assertEqual(self.theorem["failure_cover_rows_count"], 7)
        self.assertEqual(self.theorem["proper_failure_count"], 127)
        self.assertEqual(self.theorem["nearest_miss_count"], 7)
        self.assertTrue(self.theorem["all_derivative_edges_decisive"])

    def test_current_repair_and_inheritance(self) -> None:
        self.assertEqual(self.theorem["current_missing_gate_count"], 5)
        self.assertEqual(self.theorem["current_repair_distance"], 5)
        self.assertIn("source_obligation_discharge", self.theorem["current_missing_gates"])
        self.assertIn("chi11_source_export", self.theorem["current_missing_gates"])
        self.assertIn("qw2191_selector_discharge", self.theorem["current_missing_gates"])
        self.assertTrue(self.theorem["p2420_minimal_toe_masks_inherited"])
        self.assertTrue(self.theorem["p2420_repair_distance_inherited"])
        self.assertTrue(self.theorem["p2420_subcube_failure_count_inherited"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("prime-implicant", MD.read_text(encoding="utf-8"))
        self.assertIn("P2421/S1371", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2421/S1371", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
