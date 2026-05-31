from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.md"


class StrictCompletionLegacyKernelIntermediateBridgeGuardrailCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_grep_disambiguation(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_KERNEL_INTERMEDIATE_BRIDGE_GUARDRAIL_CERTIFICATE__TEXT_LOGIC_AUDIT",
        )
        self.assertIn("intermediate-bridge", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "AGENTS_guardrails",
            "K1_kernel_split_note",
            "K2_strict_gate_derivation_chain_note",
            "F2_provenance_packet",
            "F3_frontier_classification_packet",
            "S2_strategic_priority_reorientation_packet",
            "legacy_torsion_chi11_opinion_audit",
            "strict_completion_chain_integrity_report",
        })
        self.assertIn("intermediate bridge kernel", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("restoring K_legacy_ont", payload["grep_disambiguation"]["finding"])
        self.assertIn("no-silent-role-transfer", payload["grep_disambiguation"]["finding"])

    def test_evidence_rows_and_summary(self):
        payload = self.payload
        rows = payload["evidence_rows"]
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(row["passes"] for row in rows))
        self.assertEqual(
            [row["source"] for row in rows],
            [
                "AGENTS",
                "K1/K2",
                "S2",
                "S2_compression_clause",
                "F2/F3",
                "legacy_torsion_chi11_opinion_audit",
                "strict_completion_chain_integrity",
            ],
        )
        summary = payload["legacy_kernel_intermediate_bridge_summary"]
        self.assertTrue(summary["legacy_kernel_restored_as_intermediate"])
        self.assertTrue(summary["strict_kernel_treated_as_completed_legacy_continuation"])
        self.assertTrue(summary["raw_identity_bridge_still_not_silent"])
        self.assertTrue(summary["legacy_kernel_incomplete_for_strict_characteristics"])
        self.assertTrue(summary["strict_compression_missing_from_legacy_recorded"])
        self.assertTrue(summary["role_transfer_audit_required_after_full_bridge"])
        self.assertTrue(summary["beta_tors_to_chi11_remains_candidate_not_theorem"])
        self.assertTrue(summary["q2191_remains_open"])
        self.assertTrue(summary["toe_closure_not_claimed"])
        self.assertTrue(summary["intermediate_bridge_guardrail_certificate_passes"])

    def test_implications_proof_and_hard_limits(self):
        payload = self.payload
        implications = "\n".join(payload["closure_plan_implications"])
        self.assertIn("explicit completion map K_legacy_ont -> K_strict_gate", implications)
        self.assertIn("strict compression", implications)
        self.assertIn("role-transfer audit", implications)
        self.assertIn("beta_tors -> chi_11", implications)

        proof = payload["proof_certificate"]
        self.assertIn("Repo grep", proof["grep_step"])
        self.assertIn("intermediate bridge kernel", proof["guardrail_step"])
        self.assertIn("nonlinear d^eta compression", proof["compression_step"])
        self.assertIn("role-transfer audit", proof["role_transfer_step"])
        self.assertIn("QW-2191 remains open", proof["selector_step"])
        self.assertIn("bridge-completion", proof["verdict_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No raw identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No claim that the legacy kernel already contains all strict nadsoliton characteristics", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No beta_tors -> chi_11", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
