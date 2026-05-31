from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_decoder_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_decoder_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_decoder_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientCoordinateSyndromeDecoderCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COORDINATE_SYNDROME_DECODER_CERTIFICATE__FULL_DECODER_TABLE",
        )
        self.assertIn("syndrome-decoder", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_component_quotient_coordinate_isomorphism_certificate",
            "phase_sign_component_quotient_coordinate_residual_syndrome_certificate",
        })
        self.assertIn("coordinate syndrome decoder", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["coordinate_syndrome_decoder_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["coordinate_count"], 12)
        self.assertEqual(definition["node_count"], 12)
        self.assertEqual(definition["expected_coordinate_vector"], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertIn("c+delta(c)", definition["decoder_formula"])
        self.assertIn("U*T*s=s", definition["residual_decoder_formula"])

    def test_decoder_certificate(self):
        payload = self.payload
        certificate = payload["decoder_certificate"]
        self.assertEqual(certificate["coordinate_rows_checked"], 4096)
        self.assertEqual(certificate["residual_rows_checked"], 4096)
        self.assertEqual(certificate["coordinate_decode_failures"], [])
        self.assertEqual(certificate["residual_decode_failures"], [])
        self.assertEqual(len(certificate["zero_residual_decoder_rows"]), 1)
        zero_row = certificate["zero_residual_decoder_rows"][0]
        self.assertEqual(zero_row["residual_syndrome"], [0] * 12)
        self.assertEqual(zero_row["correction_delta_T_residual"], [0] * 12)
        self.assertEqual(zero_row["reencoded_residual_U_T_residual"], [0] * 12)
        self.assertTrue(certificate["decoded_coordinate_examples"])
        self.assertTrue(certificate["residual_weight_to_correction_weight_examples"])
        self.assertEqual(sum(certificate["correction_weight_histogram"].values()), 4096)

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["coordinate_syndrome_decoder_summary"]
        self.assertTrue(summary["enumerated_all_2^12_coordinate_vectors"])
        self.assertTrue(summary["enumerated_all_2^12_residual_syndromes"])
        self.assertTrue(summary["all_coordinate_vectors_decode_to_audited_coordinate"])
        self.assertTrue(summary["all_residual_syndromes_reencode_correctly"])
        self.assertTrue(summary["zero_residual_decodes_to_zero_correction"])
        self.assertTrue(summary["correction_weight_histogram_sums_to_coordinate_space"])
        self.assertTrue(summary["matches_residual_syndrome_certificate"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("c+T*r(c)=c_target", proof["coordinate_decoder_step"])
        self.assertIn("U*T*s=s", proof["residual_decoder_step"])
        self.assertIn("zero residual syndrome", proof["zero_decoder_step"])
        self.assertIn("4096 coordinate vectors", proof["histogram_step"])
        self.assertIn("does not derive omega/phi", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
