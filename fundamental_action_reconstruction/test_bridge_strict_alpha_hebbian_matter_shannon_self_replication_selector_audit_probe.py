import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_matter_shannon_self_replication_selector_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_matter_shannon_self_replication_selector_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_matter_shannon_self_replication_selector_audit_report.md"


class StrictAlphaHebbianMatterShannonSelfReplicationSelectorAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_model_and_opinion_audit(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_MATTER_SHANNON_SELF_REPLICATION_SELECTOR_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "matter-like-self-replication-can-carry-not-originate-unit-orientation-bit")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

        opinion = payload["opinion_audit"]
        self.assertEqual(
            opinion["overall_verdict"],
            "reject-bridge-closure-claim; accept-one-bit-measurement-and-conditional-carrier-claim",
        )
        verdicts = {row["claim"]: row["verdict"] for row in opinion["claims"]}
        self.assertIn("rejected_by_current_repo_guardrails", verdicts["The Legacy-Strict bridge is fully formed / 99 percent closed."])
        self.assertEqual(verdicts["Matter-like self-replicating information may explain the bit source."], "new_probe_result_conditional_carrier_not_origin")

    def test_replication_channel_audit(self):
        audit = self.payload["matter_replication_channel_audit"]
        self.assertEqual(audit["max_record_length"], 8)
        self.assertEqual(audit["error_grid_denominator"], 8)
        self.assertEqual(audit["channel_case_count"], 40)
        self.assertEqual(audit["total_channel_equivariance_violations"], 0)
        self.assertEqual(audit["total_symmetric_prior_marginal_violations"], 0)
        self.assertEqual(audit["error_free_odd_lengths_with_perfect_seed_recovery"], [1, 3, 5, 7])
        for row in audit["rows"]:
            self.assertEqual(row["channel_equivariance_violations"], 0)
            self.assertEqual(row["symmetric_prior_marginal_violations"], 0)
            self.assertIn(row["majority_decoder_interpretation"], ["tie-prone-even-length-record", "amplifies/preserves supplied seed bit"])

    def test_invariant_decoder_no_go(self):
        audit = self.payload["invariant_decoder_event_audit"]
        self.assertEqual(audit["enumerated_decoder_max_length"], 4)
        self.assertEqual(audit["enumerated_d5_majority_dominant_total"], 0)
        self.assertEqual(audit["enumerated_nonzero_majority_bias_total"], 0)
        expected_counts = [2, 4, 16, 256]
        for row, expected_count in zip(audit["rows"][:4], expected_counts):
            self.assertEqual(row["invariant_decoder_event_count"], expected_count)
            self.assertEqual(row["d5_majority_dominant_invariant_event_count"], 0)
            self.assertEqual(row["nonzero_majority_bias_invariant_event_count"], 0)
        self.assertEqual(audit["rows"][4]["invariant_decoder_event_count_formula"], "2^16")

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("B in {0,1}", proof["orientation_variable"])
        self.assertIn("P(R=r|B=b)", proof["equivariant_replication_channel"])
        self.assertIn("P(R=r)=P(R=Jr)", proof["symmetric_prior_consequence"])
        self.assertIn("Shannon information is correlation with a seed", proof["shannon_reading"])
        self.assertIn("non-invariant decoder", proof["selector_consequence"])

        interpretation = self.payload["selector_interpretation"]
        self.assertIn("self-replicating Shannon information", interpretation["question_tested"])
        self.assertIn("No mirror-equivariant replication", interpretation["negative_result"])
        self.assertIn("carry and stabilize", interpretation["conditional_positive_result"])
        self.assertIn("does not close", interpretation["legacy_strict_bridge_warning"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("downstream", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No legacy physical-role or matter-generation claim", hard_limits)
        self.assertIn("matter-like replication still forbids singleton d5", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
