from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_positive_factor_sign_separation_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_positive_factor_sign_separation_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_positive_factor_sign_separation_certificate_report.md"


class StrictCompletionPositiveFactorSignSeparationCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_rule(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_POSITIVE_FACTOR_SIGN_SEPARATION_CERTIFICATE__PHASE_ONLY_Z2_SIGN",
        )
        self.assertIn("phase-only", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "necessity_certificate",
            "phase_sign_z2_coboundary_certificate",
            "damping_exact_rational_certificate",
        })
        rule = payload["sign_separation_rule"]
        self.assertEqual(rule["factorization"], "T(d)=A(d)*P(d)*D(d)")
        self.assertIn("A(d)>0", rule["positive_factor_statement"])
        self.assertIn("sign(T(d))=sign(P(d))", rule["sign_consequence"])
        self.assertIn("zero node and edge sign bits", rule["z2_consequence"])

    def test_node_and_edge_rows(self):
        payload = self.payload
        node_rows = payload["node_sign_separation_rows"]
        self.assertEqual(len(node_rows), 12)
        self.assertTrue(all(row["alpha_normalization_positive"] for row in node_rows))
        self.assertTrue(all(row["damping_compression_positive"] for row in node_rows))
        self.assertTrue(all(row["positive_factor_z2_bit"] == 0 for row in node_rows))
        self.assertTrue(all(row["factor_sign_equals_phase_sign"] for row in node_rows))
        self.assertTrue(all(row["quotient_sign_equals_phase_sign"] for row in node_rows))
        self.assertTrue(all(row["z2_bit_equals_phase_sign_bit"] for row in node_rows))
        self.assertEqual([row["factor_product_sign"] for row in node_rows], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])

        edge_rows = payload["edge_sign_separation_rows"]
        self.assertEqual(len(edge_rows), 11)
        self.assertTrue(all(row["positive_factor_edge_bit"] == 0 for row in edge_rows))
        self.assertTrue(all(not row["positive_factors_change_edge_sign"] for row in edge_rows))
        self.assertTrue(all(row["completion_flip_equals_phase_flip"] for row in edge_rows))
        self.assertEqual([row["edge"] for row in edge_rows if row["completion_edge_bit"] == 1], ["1->2", "5->6", "7->8", "9->10"])

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["positive_factor_sign_summary"]
        self.assertTrue(summary["all_alpha_normalization_factors_positive"])
        self.assertTrue(summary["all_damping_compression_factors_positive"])
        self.assertTrue(summary["all_positive_factor_z2_bits_zero"])
        self.assertTrue(summary["all_factor_signs_equal_phase_signs"])
        self.assertTrue(summary["all_quotient_signs_equal_phase_signs"])
        self.assertTrue(summary["all_z2_bits_equal_phase_sign_bits"])
        self.assertTrue(summary["all_positive_factor_edge_bits_zero"])
        self.assertTrue(summary["all_completion_flips_equal_phase_flips"])
        self.assertEqual(summary["derived_completion_sign_pattern"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertEqual(summary["derived_completion_flip_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertTrue(summary["matches_z2_sign_pattern"])
        self.assertTrue(summary["matches_z2_flip_edges"])
        self.assertTrue(summary["exact_damping_positive_and_decreasing"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("strict_damping_formula_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("strict_transport_derivation_from_nadsoliton_dynamics", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("T(d)=A(d)*P(d)*D(d)", proof["factorization_step"])
        self.assertIn("positive at every audited node", proof["positivity_step"])
        self.assertIn("sign(T(d))=sign(P(d))", proof["node_sign_step"])
        self.assertIn("zero Z2 edge bits", proof["edge_z2_step"])
        self.assertIn("does not derive A(d)", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No proof derives A(d), P(d), D(d)", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
