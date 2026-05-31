from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_generator_basis_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_generator_basis_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_generator_basis_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientCoordinateSyndromeGeneratorBasisCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COORDINATE_SYNDROME_GENERATOR_BASIS_CERTIFICATE__HYPERCUBE_EDGE_CHECKS",
        )
        self.assertIn("generator-basis", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_component_quotient_coordinate_isomorphism_certificate",
            "phase_sign_component_quotient_coordinate_syndrome_decoder_certificate",
        })
        self.assertIn("coordinate syndrome generator basis", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["generator_basis_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["coordinate_count"], 12)
        self.assertEqual(definition["node_count"], 12)
        self.assertIn("r(c+e_i)+r(c)=U_i", definition["generator_formula"])
        self.assertIn("T*U_i=e_i", definition["decoder_generator_formula"])

    def test_generator_basis_and_hypercube_edges(self):
        payload = self.payload
        certificate = payload["generator_basis_certificate"]
        self.assertEqual(len(certificate["generator_rows"]), 12)
        self.assertEqual(certificate["residual_generator_column_rank"], 12)
        self.assertEqual(certificate["decoded_generator_row_rank"], 12)
        self.assertEqual(len(certificate["residual_generator_weights"]), 12)
        for index, row in enumerate(certificate["generator_rows"]):
            self.assertEqual(row["coordinate_index"], index)
            self.assertEqual(row["decoded_T_U_i"], row["coordinate_unit_e_i"])
            self.assertTrue(row["decoded_matches_unit"])
            self.assertGreater(row["residual_generator_weight"], 0)

        edge_certificate = payload["hypercube_edge_certificate"]
        self.assertEqual(edge_certificate["hypercube_edge_checks"], 4096 * 12)
        self.assertEqual(edge_certificate["hypercube_edge_failures"], [])
        self.assertTrue(edge_certificate["hypercube_edge_examples"])

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["coordinate_syndrome_generator_basis_summary"]
        self.assertTrue(summary["all_12_generators_checked"])
        self.assertTrue(summary["all_generators_decode_to_coordinate_units"])
        self.assertTrue(summary["residual_generators_have_full_rank_12"])
        self.assertTrue(summary["decoded_generators_have_full_rank_12"])
        self.assertTrue(summary["checked_all_4096_times_12_hypercube_edges"])
        self.assertTrue(summary["all_hypercube_edge_deltas_match_generators"])
        self.assertTrue(summary["matches_syndrome_decoder_report"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("T*U_i=e_i", proof["generator_step"])
        self.assertIn("rank 12", proof["rank_step"])
        self.assertIn("4096*12 hypercube edges", proof["edge_step"])
        self.assertIn("c+T*r(c)=c_target", proof["decoder_link_step"])
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
