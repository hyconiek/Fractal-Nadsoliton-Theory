#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.py"
OUT = ROOT / "generated" / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.json"
MD = ROOT / "generated" / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.md"
P2412 = ROOT / "generated" / "p2412_s1362_chi11_selector_scope_separation_certificate.json"


class P2413AmplitudeScalarNormalizationBridgeWitnessCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2412.exists():
            subprocess.run([sys.executable, str(ROOT / "p2412_s1362_chi11_selector_scope_separation_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["amplitude_scalar_normalization_bridge_witness_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_witness_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2413")
        self.assertEqual(self.payload["stage_id"], "S1363")
        self.assertEqual(self.payload["status"], "PASS_AMPLITUDE_SCALAR_WITNESS_NO_ROLE_TRANSFER_NO_FULL_BRIDGE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_exact_scalar_witness(self) -> None:
        self.assertEqual(self.theorem["domain_size"], 12)
        self.assertEqual(self.theorem["normalization_map"], "N_alpha: K -> alpha_geo^{-1} K")
        self.assertEqual(self.theorem["cos_zero_solution_count"], 0)
        self.assertTrue(self.theorem["all_denominators_positive_exact"])
        self.assertTrue(self.theorem["all_signs_preserved_by_positive_alpha"])
        self.assertLess(self.theorem["max_float_residual"], 1e-15)
        self.assertTrue(self.theorem["scalar_normalization_witness_ready"])

    def test_rows_and_proof_certificate(self) -> None:
        rows = self.finite["witness_rows"]
        self.assertEqual([row["d"] for row in rows], list(range(12)))
        self.assertTrue(all(not row["cos_zero_congruence_3d_eq_4_mod_12"] for row in rows))
        self.assertTrue(all(row["denominator_positive_exact"] for row in rows))
        proof = self.finite["proof_certificate"]
        self.assertIn("K_legacy_ont(d)=alpha_geo*L_shape(d)", proof["factorization_step"])
        self.assertIn("3d == 4 mod 12", proof["nonzero_cos_step"])

    def test_hard_limits_and_inheritance(self) -> None:
        self.assertFalse(self.theorem["full_amplitude_source_theorem_ready"])
        self.assertFalse(self.theorem["role_safe_amplitude_absorption_ready"])
        self.assertFalse(self.theorem["bridge_completion_ready"])
        self.assertTrue(self.theorem["p2411_full_bridge_still_requires_all_obligations"])
        self.assertTrue(self.theorem["p2412_chi11_scope_separation_inherited"])
        self.assertTrue(self.theorem["scratch_amplitude_witness_inherited"])
        self.assertTrue(self.theorem["scratch_blocks_role_transfer_inherited"])
        self.assertIn("No role-safe alpha_geo", "\n".join(self.theorem["not_licensed"]))

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("amplitude scalar-normalization", MD.read_text(encoding="utf-8"))
        self.assertIn("P2413/S1363 amplitude", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2413/S1363 amplitude", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
