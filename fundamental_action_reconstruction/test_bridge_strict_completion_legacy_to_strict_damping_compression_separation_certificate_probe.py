from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_report.md"


class LegacyToStrictDampingCompressionSeparationCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_nonduplication(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_DAMPING_COMPRESSION_SEPARATION__LINEAR_TORSION_NO_GO",
        )
        self.assertIn("linear-torsion-denominator-separated", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "component_gap_matrix",
            "necessity",
            "damping_exact",
            "legacy_bridge_guardrail",
            "compression_ontology_audit",
            "t15_bridge_theorem_spec",
        })
        searched = "\n".join(payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("1+gamma*d vs 1+d^eta", searched)
        self.assertIn("beta_tors*d to beta*d^eta", searched)
        self.assertIn("no single linear gamma", payload["grep_disambiguation"]["finding"])

    def test_node_rows_and_separation_summary(self):
        payload = self.payload
        rows = payload["node_rows"]
        self.assertEqual([row["d"] for row in rows], list(range(1, 12)))
        self.assertEqual(rows[0]["required_gamma_for_exact_single_node_match"], 1.0)
        self.assertGreater(rows[-1]["required_gamma_for_exact_single_node_match"], rows[0]["required_gamma_for_exact_single_node_match"])
        self.assertFalse(any(row["legacy_beta_tors_matches_this_node"] for row in rows))
        self.assertTrue(all(row["strict_minus_legacy_denominator_residual"] > 0 for row in rows))

        summary = payload["separation_summary"]
        self.assertEqual(summary["domain"], list(range(1, 12)))
        self.assertEqual(summary["eta"], 1.8)
        self.assertAlmostEqual(summary["eta_minus_one"], 0.8)
        self.assertTrue(summary["required_gamma_values_strictly_increase"])
        self.assertTrue(summary["no_single_linear_gamma_matches_two_distinct_positive_nodes"])
        self.assertTrue(summary["legacy_beta_tors_matches_no_positive_strict_node"])
        self.assertGreater(summary["best_l2_residual_l2_norm"], 0.0)
        self.assertGreater(summary["best_l2_residual_max_abs"], 0.0)
        self.assertGreater(summary["required_gamma_spread_d11_minus_d1"], 0.0)
        self.assertGreater(summary["d11_strict_over_legacy_denominator_ratio"], 1.0)
        self.assertTrue(summary["component_gap_records_compression_missing"])
        self.assertTrue(summary["guardrail_records_legacy_incomplete"])
        self.assertFalse(summary["full_bridge_theorem_exported"])
        self.assertTrue(payload["all_cross_checks_pass"])
        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        proof = self.payload["proof_certificate"]
        self.assertIn("rg was used", proof["nonduplication_step"])
        self.assertIn("gamma=d^(4/5)", proof["algebraic_step"])
        self.assertIn("no constant gamma", proof["algebraic_step"])
        self.assertIn("beta_tors=0.01", proof["legacy_beta_tors_step"])
        self.assertIn("best finite L2 linear denominator", proof["least_squares_step"])
        self.assertIn("not hidden inside beta_tors*d", proof["bridge_meaning_step"])
        self.assertIn("not a strict-source derivation", proof["theoretical_limit"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No raw identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No derivation of eta=9/5", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
