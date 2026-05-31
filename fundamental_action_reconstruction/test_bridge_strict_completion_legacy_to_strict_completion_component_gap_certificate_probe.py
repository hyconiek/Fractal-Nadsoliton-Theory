from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.md"


class LegacyToStrictCompletionComponentGapCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_nonduplication(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_COMPONENT_GAP_MATRIX__NO_FULL_BRIDGE_THEOREM",
        )
        self.assertIn("component-gap-matrix", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "agents_guardrail",
            "s2_priority_packet",
            "necessity",
            "cocycle",
            "gf2_linear_system",
            "gf2_commutative_diagram",
            "path_cohomology",
            "cycle_closure",
            "damping_exact",
            "legacy_bridge_guardrail",
            "torsion_chi11_audit",
            "s1_obstruction",
            "closure_plan",
        })
        searched = "\n".join(payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("completion component gap", searched)
        self.assertIn("damping/compression passage", searched)
        self.assertIn("role-transfer audit", searched)
        self.assertIn("no single component-level legacy->strict bridge gap matrix", payload["grep_disambiguation"]["finding"])

    def test_matrix_rows_summary_and_cross_checks(self):
        payload = self.payload
        rows = payload["component_gap_rows"]
        matrix = payload["component_gap_matrix"]
        columns = payload["matrix_definition"]["columns"]
        self.assertEqual(len(rows), 5)
        self.assertEqual(len(matrix), 5)
        self.assertEqual(columns, [
            "legacy_input_visible",
            "strict_completion_visible",
            "finite_certificate_exported",
            "strict_dynamic_source_exported",
            "selector_or_source_exported",
            "role_transfer_allowed_now",
        ])
        self.assertTrue(all(len(row) == len(columns) for row in matrix))
        self.assertEqual([row["component"] for row in rows], [
            "amplitude_normalization",
            "phase_frequency_transport",
            "damping_compression",
            "topological_phase_bit_chi11",
            "legacy_physical_role_transfer",
        ])
        summary = payload["completion_gap_summary"]
        self.assertEqual(summary["component_count"], 5)
        self.assertEqual(summary["rows_with_finite_certificates"], 5)
        self.assertEqual(summary["rows_with_strict_dynamic_sources"], 0)
        self.assertEqual(summary["rows_with_role_transfer_allowed_now"], 0)
        self.assertEqual(summary["component_matrix_rank_mod2"], 2)
        self.assertTrue(summary["all_rows_have_finite_certificates"])
        self.assertTrue(summary["strict_dynamic_sources_missing"])
        self.assertTrue(summary["selector_source_gap_remains"])
        self.assertTrue(summary["role_transfer_blocked_until_full_bridge"])
        self.assertFalse(summary["bridge_ready_for_role_transfer"])
        self.assertFalse(summary["raw_identity_claimed"])
        self.assertTrue(summary["completion_map_partial_not_full"])
        self.assertTrue(summary["strict_compression_recorded_as_missing_from_legacy"])
        self.assertTrue(summary["beta_tors_to_chi11_candidate_not_theorem"])
        self.assertTrue(payload["all_cross_checks_pass"])
        self.assertTrue(all(payload["cross_checks"].values()))

    def test_component_gaps_and_proof_guardrails(self):
        rows_by_name = {row["component"]: row for row in self.payload["component_gap_rows"]}
        self.assertIn("strict nonlinear compression", rows_by_name["damping_compression"]["bridge_obligation"])
        self.assertIn("d^eta compression", rows_by_name["damping_compression"]["open_gap"])
        self.assertIn("beta_tors", rows_by_name["topological_phase_bit_chi11"]["bridge_obligation"])
        self.assertIn("candidate bridge hypothesis", "\n".join(rows_by_name["topological_phase_bit_chi11"]["finite_evidence"]))
        self.assertIn("sin^2(theta_W)", rows_by_name["legacy_physical_role_transfer"]["open_gap"])
        self.assertFalse(rows_by_name["legacy_physical_role_transfer"]["matrix_bits"]["role_transfer_allowed_now"])

        proof = self.payload["proof_certificate"]
        self.assertIn("rg was used", proof["nonduplication_step"])
        self.assertIn("rank 2", proof["matrix_step"])
        self.assertIn("strict nonlinear d^eta compression", proof["compression_step"])
        self.assertIn("beta_tors->chi_11 remains only a candidate", proof["selector_step"])
        self.assertIn("No legacy physical role is transferable now", proof["role_transfer_step"])
        self.assertIn("not ToE closure", proof["theoretical_limit"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No raw identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No strict dynamical derivation", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
