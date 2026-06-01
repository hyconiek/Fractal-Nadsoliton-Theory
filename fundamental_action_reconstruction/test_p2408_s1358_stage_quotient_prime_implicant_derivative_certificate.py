#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2408_s1358_stage_quotient_prime_implicant_derivative_certificate.py"
OUT = ROOT / "generated" / "p2408_s1358_stage_quotient_prime_implicant_derivative_certificate.json"
MD = ROOT / "generated" / "p2408_s1358_stage_quotient_prime_implicant_derivative_certificate.md"
P2407 = ROOT / "generated" / "p2407_s1357_stage_quotient_projection_barrier_certificate.json"


class P2408StageQuotientPrimeImplicantDerivativeCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2407.exists():
            subprocess.run([sys.executable, str(ROOT / "p2407_s1357_stage_quotient_projection_barrier_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["stage_quotient_prime_implicant_derivative_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_boolean_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2408")
        self.assertEqual(self.payload["stage_id"], "S1358")
        self.assertEqual(self.payload["status"], "PASS_UNIQUE_PRIME_IMPLICANT_AND_SINGLE_EDGE_DERIVATIVES_NO_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_truth_vector_and_anf(self) -> None:
        self.assertEqual(self.theorem["truth_vector_by_mask_0_to_7"], [0, 0, 0, 0, 0, 0, 0, 1])
        self.assertEqual(self.theorem["true_masks"], [7])
        self.assertEqual(
            self.theorem["anf_polynomial"],
            "O_ontology_guard_package * S_strict_internal_completion_package * R_role_successor_projection_package",
        )
        self.assertEqual(self.theorem["anf_degree"], 3)

    def test_unique_prime_implicant(self) -> None:
        self.assertEqual(self.theorem["prime_implicant_count"], 1)
        prime = self.theorem["prime_implicants"][0]
        self.assertEqual(prime["required_ones"], self.finite["stage_atoms"])
        self.assertEqual(prime["required_zeros"], [])
        self.assertEqual(prime["support_masks"], [7])
        self.assertEqual(prime["literal_count"], 3)

    def test_boolean_derivative_edges(self) -> None:
        self.assertEqual(
            self.theorem["boolean_derivative_edge_counts"],
            {
                "O_ontology_guard_package": 1,
                "S_strict_internal_completion_package": 1,
                "R_role_successor_projection_package": 1,
            },
        )
        self.assertEqual(sorted(self.theorem["nearest_miss_masks_by_stage"].values()), [3, 5, 6])
        self.assertTrue(self.finite["all_derivatives_have_single_edge_support"])

    def test_inherits_p2407_limits(self) -> None:
        self.assertTrue(self.theorem["p2407_full_mask_only_inherited"])
        self.assertTrue(self.theorem["p2407_stage_degree_three_inherited"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("prime-implicant derivative certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2408/S1358 stage-quotient", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2408/S1358 prime-implicant", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
