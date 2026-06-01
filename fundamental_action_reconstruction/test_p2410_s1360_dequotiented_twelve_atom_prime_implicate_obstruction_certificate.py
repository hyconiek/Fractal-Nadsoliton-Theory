#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2410_s1360_dequotiented_twelve_atom_prime_implicate_obstruction_certificate.py"
OUT = ROOT / "generated" / "p2410_s1360_dequotiented_twelve_atom_prime_implicate_obstruction_certificate.json"
MD = ROOT / "generated" / "p2410_s1360_dequotiented_twelve_atom_prime_implicate_obstruction_certificate.md"
P2406 = ROOT / "generated" / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.json"
P2409 = ROOT / "generated" / "p2409_s1359_stage_quotient_prime_implicate_failure_cover_certificate.json"


class P2410DequotientedTwelveAtomPrimeImplicateObstructionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2406.exists():
            subprocess.run([sys.executable, str(ROOT / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.py")], check=True)
        if not P2409.exists():
            subprocess.run([sys.executable, str(ROOT / "p2409_s1359_stage_quotient_prime_implicate_failure_cover_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["dequotiented_twelve_atom_prime_implicate_obstruction_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_obstruction_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2410")
        self.assertEqual(self.payload["stage_id"], "S1360")
        self.assertEqual(self.payload["status"], "PASS_DEQUOTIENTED_12_ATOM_OBSTRUCTION_LEDGER_NO_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_assignment_and_clause_counts(self) -> None:
        self.assertEqual(self.theorem["atom_count"], 12)
        self.assertEqual(self.theorem["total_assignment_count"], 4096)
        self.assertEqual(self.theorem["success_assignment_count"], 1)
        self.assertEqual(self.theorem["failure_assignment_count"], 4095)
        self.assertEqual(self.theorem["enumerated_nonempty_clause_count"], (3**12) - 1)
        self.assertEqual(self.theorem["valid_implicate_count"], (3**12) - (2**12))

    def test_prime_implicates_and_failure_cover(self) -> None:
        self.assertEqual(self.theorem["success_cnf_unit_count"], 12)
        self.assertEqual(self.theorem["failure_cover_term_count"], 12)
        self.assertEqual(self.theorem["prime_implicate_count"], 12)
        self.assertEqual(self.theorem["prime_implicate_atoms"], self.finite["stage_packages"]["O_ontology_guard_package"] + self.finite["stage_packages"]["S_strict_internal_completion_package"] + self.finite["stage_packages"]["R_role_successor_projection_package"])
        self.assertEqual(len(self.finite["failure_cover"]), 12)
        self.assertTrue(all(row["covered_failure_mask_count"] == 2048 for row in self.finite["failure_cover"]))

    def test_repair_spectrum(self) -> None:
        spectrum = self.theorem["repair_spectrum"]
        self.assertEqual(len(spectrum), 12)
        self.assertEqual(sum(row["failure_mask_count"] for row in spectrum), 4095)
        self.assertEqual(spectrum[0]["missing_atom_count"], 1)
        self.assertEqual(spectrum[0]["failure_mask_count"], 12)
        self.assertEqual(spectrum[-1]["missing_atom_count"], 12)
        self.assertEqual(spectrum[-1]["failure_mask_count"], 1)
        self.assertEqual(self.theorem["nearest_one_atom_missing_failure_count"], 12)

    def test_stage_preimages_and_inheritance(self) -> None:
        self.assertEqual(self.theorem["stage_quotient_full_mask_preimage_count"], 1)
        self.assertEqual(self.theorem["stage_quotient_preimage_counts"]["0"], (2**5 - 1) * (2**4 - 1) * (2**3 - 1))
        self.assertTrue(self.theorem["p2406_degree_twelve_inherited"])
        self.assertTrue(self.theorem["p2409_three_stage_cnf_inherited"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("dequotiented twelve-atom", MD.read_text(encoding="utf-8"))
        self.assertIn("P2410/S1360 dequotiented", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2410/S1360 dequotiented", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
