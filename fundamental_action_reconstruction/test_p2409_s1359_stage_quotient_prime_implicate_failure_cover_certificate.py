#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2409_s1359_stage_quotient_prime_implicate_failure_cover_certificate.py"
OUT = ROOT / "generated" / "p2409_s1359_stage_quotient_prime_implicate_failure_cover_certificate.json"
MD = ROOT / "generated" / "p2409_s1359_stage_quotient_prime_implicate_failure_cover_certificate.md"
P2408 = ROOT / "generated" / "p2408_s1358_stage_quotient_prime_implicant_derivative_certificate.json"


class P2409StageQuotientPrimeImplicateFailureCoverCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2408.exists():
            subprocess.run([sys.executable, str(ROOT / "p2408_s1358_stage_quotient_prime_implicant_derivative_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["stage_quotient_prime_implicate_failure_cover_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_dual_boolean_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2409")
        self.assertEqual(self.payload["stage_id"], "S1359")
        self.assertEqual(self.payload["status"], "PASS_PRIME_IMPLICATE_FAILURE_COVER_NO_SHORTCUT_NO_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_success_prime_implicate_cnf(self) -> None:
        self.assertEqual(self.theorem["truth_vector_by_mask_0_to_7"], [0, 0, 0, 0, 0, 0, 0, 1])
        self.assertEqual(
            self.theorem["success_cnf"],
            "O_ontology_guard_package AND S_strict_internal_completion_package AND R_role_successor_projection_package",
        )
        self.assertEqual(self.theorem["success_prime_implicate_count"], 3)
        self.assertTrue(self.finite["prime_implicates_are_positive_stage_units"])

    def test_failure_cover(self) -> None:
        self.assertEqual(
            self.theorem["failure_dnf"],
            "not O_ontology_guard_package OR not S_strict_internal_completion_package OR not R_role_successor_projection_package",
        )
        self.assertEqual(self.theorem["failure_cover_term_count"], 3)
        self.assertEqual(self.theorem["shortcut_failure_mask_count"], 7)
        self.assertEqual(self.theorem["shortcut_masks"], [0, 1, 2, 3, 4, 5, 6])
        self.assertTrue(self.finite["all_shortcuts_hit_by_failure_cover"])
        self.assertTrue(self.finite["full_mask_avoids_failure_cover"])

    def test_repair_distances_and_p2408_inheritance(self) -> None:
        self.assertEqual(self.theorem["one_stage_missing_masks"], [3, 5, 6])
        self.assertEqual(self.theorem["max_nearest_repair_distance"], 3)
        self.assertTrue(self.theorem["p2408_single_success_prime_implicant_inherited"])
        self.assertTrue(self.theorem["p2408_derivative_nearest_misses_inherited"])

    def test_shortcut_witness_rows_are_explicit(self) -> None:
        table = self.finite["shortcut_witness_table"]
        self.assertEqual(len(table), 7)
        empty = table[0]
        self.assertEqual(empty["mask"], 0)
        self.assertEqual(len(empty["minimal_repair_options"]), 3)
        nearest = [row for row in table if row["nearest_repair_distance"] == 1]
        self.assertEqual([row["mask"] for row in nearest], [3, 5, 6])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("prime-implicate failure-cover certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2409/S1359 stage-quotient", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2409/S1359 prime-implicate", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
