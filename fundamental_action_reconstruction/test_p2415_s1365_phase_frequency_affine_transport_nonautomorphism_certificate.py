#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.py"
OUT = ROOT / "generated" / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.json"
MD = ROOT / "generated" / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.md"
P2414 = ROOT / "generated" / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.json"


class P2415PhaseFrequencyAffineTransportNonautomorphismCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2414.exists():
            subprocess.run([sys.executable, str(ROOT / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["phase_frequency_affine_transport_nonautomorphism_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_witness_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2415")
        self.assertEqual(self.payload["stage_id"], "S1365")
        self.assertEqual(self.payload["status"], "PASS_PHASE_FREQUENCY_AFFINE_TRANSPORT_WITNESS_NO_SOURCE_NO_SELECTOR_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_affine_transport_witness(self) -> None:
        self.assertEqual(self.theorem["domain_size"], 12)
        self.assertTrue(self.theorem["continuous_affine_phase_transport_exact_float"])
        self.assertLessEqual(self.theorem["max_abs_affine_transport_residual"], 1e-14)
        rows = self.finite["phase_transport_rows"]
        self.assertEqual([row["d"] for row in rows], list(range(12)))
        self.assertTrue(any(row["distance_x_d_to_nearest_integer"] > 1e-14 for row in rows))

    def test_nonautomorphism_and_nonscalar_checks(self) -> None:
        self.assertTrue(self.theorem["affine_coordinates_not_all_integers"])
        self.assertEqual(self.theorem["z12_unit_offset_automorphism_count_checked"], 48)
        self.assertTrue(self.theorem["no_z12_unit_offset_reindex_matches_strict_sign_pattern"])
        self.assertGreater(self.theorem["best_z12_unit_offset_mismatch_count"], 0)
        self.assertTrue(self.theorem["scalar_phase_replacement_fails"])
        self.assertGreater(self.theorem["scalar_phase_best_fit_max_abs_residual"], 1e-3)

    def test_phase_bits_and_inherited_limits(self) -> None:
        self.assertEqual(len(self.theorem["phase_factor_bits"]), 12)
        self.assertTrue(self.theorem["phase_factor_bits_match_z2_node_bits"])
        self.assertTrue(self.theorem["gf2_solution_inherited_unique"])
        self.assertTrue(self.theorem["p2411_full_bridge_still_requires_all_obligations"])
        self.assertTrue(self.theorem["p2412_chi11_scope_separation_inherited"])
        self.assertTrue(self.theorem["p2413_amplitude_witness_inherited"])
        self.assertTrue(self.theorem["p2414_damping_nonabsorption_inherited"])
        self.assertTrue(self.theorem["scratch_affine_transport_inherited"])
        self.assertTrue(self.theorem["scratch_nonautomorphism_inherited"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["orientation_selector_source_exported"])
        self.assertFalse(self.theorem["phase_frequency_bridge_component_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertIn("No strict dynamic source", "\n".join(self.theorem["not_licensed"]))
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("phase/frequency affine", MD.read_text(encoding="utf-8"))
        self.assertIn("P2415/S1365 phase/frequency", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2415/S1365 phase/frequency", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
