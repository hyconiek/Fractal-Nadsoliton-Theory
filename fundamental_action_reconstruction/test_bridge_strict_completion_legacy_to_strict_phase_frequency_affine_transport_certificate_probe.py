from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_phase_frequency_affine_transport_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_phase_frequency_affine_transport_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_phase_frequency_affine_transport_certificate_report.md"


class LegacyToStrictPhaseFrequencyAffineTransportCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_PHASE_FREQUENCY_AFFINE_TRANSPORT_CERTIFICATE__NO_SELECTOR_SOURCE",
        )
        self.assertIn("continuous-affine-phase-transport-exact", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "necessity",
            "phase_sign_z2_coboundary",
            "gf2_linear_system",
            "component_gap_matrix",
            "finite_diagonal_completion_map",
            "legacy_bridge_guardrail",
        })
        self.assertIn("phase affine", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("Z12 automorphism", payload["grep_disambiguation"]["finding"])

        definition = payload["phase_affine_definition"]
        self.assertIn("theta_L", definition["legacy_argument"])
        self.assertIn("theta_S", definition["strict_argument"])
        self.assertIn("omega_S/omega_L", definition["affine_transport"])
        self.assertEqual(definition["z12_units_checked"], [1, 5, 7, 11])

    def test_rows_summary_and_cross_checks(self):
        payload = self.payload
        rows = payload["phase_transport_rows"]
        self.assertEqual(len(rows), 12)
        self.assertEqual([row["d"] for row in rows], list(range(12)))
        self.assertTrue(all(abs(row["affine_transport_residual"]) <= 1e-14 for row in rows))
        self.assertTrue(all(row["matches_z2_node_bit"] for row in rows))
        self.assertEqual([row["phase_factor_gf2_bit"] for row in rows], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])

        automorphism_rows = payload["z12_unit_offset_automorphism_rows"]
        self.assertEqual(len(automorphism_rows), 48)
        self.assertTrue(all(not row["matches_strict_sign_pattern"] for row in automorphism_rows))
        self.assertEqual(min(row["mismatch_count"] for row in automorphism_rows), 2)

        summary = payload["phase_frequency_affine_transport_summary"]
        self.assertEqual(summary["domain_size"], 12)
        self.assertTrue(summary["continuous_affine_phase_transport_exact"])
        self.assertLessEqual(summary["max_abs_affine_transport_residual"], 1e-14)
        self.assertTrue(summary["affine_map_is_not_z12_automorphism"])
        self.assertTrue(summary["affine_coordinates_not_all_integers"])
        self.assertEqual(summary["z12_unit_offset_automorphism_count_checked"], 48)
        self.assertTrue(summary["no_z12_unit_offset_reindex_matches_strict_sign_pattern"])
        self.assertEqual(summary["best_z12_unit_offset_mismatch_count"], 2)
        self.assertTrue(summary["phase_factor_bits_match_z2_node_bits"])
        self.assertTrue(summary["gf2_solution_inherited_unique"])
        self.assertTrue(summary["scalar_phase_replacement_fails"])
        self.assertGreater(summary["scalar_phase_best_fit_max_abs_residual"], 0.1)
        self.assertTrue(summary["finite_diagonal_completion_map_inherited"])
        self.assertTrue(summary["component_gap_phase_row_still_source_open"])
        self.assertFalse(summary["strict_phase_frequency_source_exported"])
        self.assertFalse(summary["orientation_selector_source_exported"])
        self.assertFalse(summary["raw_kernel_identity_claimed"])
        self.assertFalse(summary["full_bridge_theorem_exported"])

        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("theta_L(x(d))", proof["affine_step"])
        self.assertIn("not a Z12 unit", proof["non_automorphism_step"])
        self.assertIn("P(d)=cos(theta_S(d))/cos(theta_L(d))", proof["phase_factor_step"])
        self.assertIn("not scalar normalization", proof["non_scalar_step"])
        self.assertIn("not a strict derivation of omega/phi", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No raw identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No strict dynamical source for omega/phi", hard_limits)
        self.assertIn("No orientation/selector source", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
