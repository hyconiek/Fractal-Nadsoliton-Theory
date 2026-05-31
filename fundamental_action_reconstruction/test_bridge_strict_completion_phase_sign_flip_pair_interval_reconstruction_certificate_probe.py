from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_flip_pair_interval_reconstruction_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_flip_pair_interval_reconstruction_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_flip_pair_interval_reconstruction_certificate_report.md"


class StrictCompletionPhaseSignFlipPairIntervalReconstructionCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_FLIP_PAIR_INTERVAL_RECONSTRUCTION_CERTIFICATE__BOUNDARY_INVERSE_ON_PATH",
        )
        self.assertIn("ordered-flip-pairs", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_edge_support_minimality_certificate",
            "phase_sign_node_support_interval_boundary_certificate",
            "phase_sign_reduced_coboundary_inverse_certificate",
        })
        self.assertIn("flip pair", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["flip_pair_reconstruction_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["anchor_bit_b0"], 0)
        self.assertEqual(definition["flip_indices"], [1, 5, 7, 9])
        self.assertEqual(definition["flip_edges"], ["1->2", "5->6", "7->8", "9->10"])

    def test_scan_rows_pair_rows_and_reconstruction(self):
        payload = self.payload
        scan_rows = payload["flip_scan_rows"]
        self.assertEqual([row["role"] for row in scan_rows], ["entry", "exit", "entry", "exit"])
        self.assertEqual([row["state_before"] for row in scan_rows], [0, 1, 0, 1])
        self.assertEqual([row["state_after"] for row in scan_rows], [1, 0, 1, 0])

        pair_rows = payload["paired_interval_rows"]
        self.assertEqual(len(pair_rows), 2)
        self.assertEqual(pair_rows[0]["interval"], {"start": 2, "end": 5})
        self.assertEqual(pair_rows[1]["interval"], {"start": 8, "end": 9})
        self.assertEqual(pair_rows[0]["boundary_edges"], ["1->2", "5->6"])
        self.assertEqual(pair_rows[1]["boundary_edges"], ["7->8", "9->10"])
        self.assertTrue(all(row["boundary_edges_match_edge_bits"] for row in pair_rows))

        reconstruction = payload["reconstruction_certificate"]
        self.assertEqual(reconstruction["reconstructed_intervals_from_flip_pairs"], [{"start": 2, "end": 5}, {"start": 8, "end": 9}])
        self.assertEqual(reconstruction["node_bits_from_anchor_scan"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(reconstruction["edge_bits_from_reconstructed_interval_boundaries"], [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])
        self.assertFalse(reconstruction["has_unclosed_interval"])

    def test_parity_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        parity = payload["parity_and_endpoint_certificate"]
        self.assertEqual(parity["flip_count"], 4)
        self.assertTrue(parity["flip_count_even"])
        self.assertEqual(parity["component_count_from_pairs"], 2)
        self.assertTrue(parity["flip_count_equals_two_components"])
        self.assertTrue(parity["no_endpoint_support"])
        self.assertTrue(parity["scan_returns_to_anchor_state"])

        summary = payload["flip_pair_interval_reconstruction_summary"]
        self.assertTrue(summary["matches_expected_flip_indices"])
        self.assertTrue(summary["matches_expected_intervals"])
        self.assertTrue(summary["node_bits_from_scan_match_z2"])
        self.assertTrue(summary["edge_bits_from_interval_boundaries_match_z2"])
        self.assertTrue(summary["matches_interval_boundary_report"])
        self.assertTrue(summary["matches_edge_support_minimality"])
        self.assertTrue(summary["matches_reduced_coboundary_nodes"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("[1,5,7,9]", proof["scan_step"])
        self.assertIn("[2,5] and [8,9]", proof["pairing_step"])
        self.assertIn("GF(2) boundaries", proof["boundary_inverse_step"])
        self.assertIn("2*component_count", proof["parity_step"])
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
