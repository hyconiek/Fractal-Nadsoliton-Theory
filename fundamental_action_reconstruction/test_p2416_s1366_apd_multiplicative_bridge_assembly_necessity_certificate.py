#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.py"
OUT = ROOT / "generated" / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json"
MD = ROOT / "generated" / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.md"
P2415 = ROOT / "generated" / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.json"


class P2416APDMultiplicativeBridgeAssemblyNecessityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2415.exists():
            subprocess.run([sys.executable, str(ROOT / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["apd_multiplicative_bridge_assembly_necessity_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_witness_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2416")
        self.assertEqual(self.payload["stage_id"], "S1366")
        self.assertEqual(self.payload["status"], "PASS_APD_VALUE_ASSEMBLY_EXACT_NO_SOURCE_NO_ROLE_TRANSFER")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_pointwise_factorization(self) -> None:
        self.assertEqual(self.theorem["domain_size"], 12)
        self.assertEqual(self.theorem["factor_names"], ["alpha_normalization", "phase_frequency_transport", "damping_compression"])
        self.assertLessEqual(self.theorem["max_abs_full_factorization_residual"], 1e-12)
        rows = self.finite["pointwise_factor_rows"]
        self.assertEqual([row["d"] for row in rows], list(range(12)))
        self.assertTrue(all(abs(row["quotient_minus_factor_product"]) <= 1e-12 for row in rows))

    def test_subset_lattice_necessity(self) -> None:
        self.assertEqual(self.theorem["subset_count"], 8)
        self.assertEqual(
            self.theorem["exact_subsets_without_extra_scalar"],
            ["alpha_normalization+phase_frequency_transport+damping_compression"],
        )
        self.assertIn(
            "phase_frequency_transport+damping_compression",
            self.theorem["exact_subsets_up_to_one_global_scalar"],
        )
        self.assertTrue(self.theorem["unique_exact_without_extra_scalar"])
        self.assertTrue(self.theorem["alpha_missing_repairable_only_by_global_scalar"])
        self.assertTrue(self.theorem["phase_missing_not_scalar_repairable"])
        self.assertTrue(self.theorem["damping_missing_not_scalar_repairable"])

    def test_hard_limits_and_inheritance(self) -> None:
        self.assertTrue(self.theorem["scratch_necessity_inherited"])
        self.assertTrue(self.theorem["p2411_full_bridge_still_requires_all_obligations"])
        self.assertTrue(self.theorem["p2413_amplitude_witness_inherited"])
        self.assertTrue(self.theorem["p2414_damping_nonabsorption_inherited"])
        self.assertTrue(self.theorem["p2415_phase_transport_inherited"])
        self.assertTrue(self.theorem["apd_value_assembly_witness_ready"])
        self.assertFalse(self.theorem["strict_dynamic_source_exported"])
        self.assertFalse(self.theorem["selector_source_exported"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertIn("No strict dynamic source", "\n".join(self.theorem["not_licensed"]))

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("APD multiplicative", MD.read_text(encoding="utf-8"))
        self.assertIn("P2416/S1366 APD", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2416/S1366 APD", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
