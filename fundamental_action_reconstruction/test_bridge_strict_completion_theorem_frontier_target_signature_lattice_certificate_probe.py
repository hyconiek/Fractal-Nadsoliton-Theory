import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_target_signature_lattice_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_target_signature_lattice_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_target_signature_lattice_certificate_report.md"


class TheoremFrontierTargetSignatureLatticeCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_targets(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_THEOREM_FRONTIER_TARGET_SIGNATURE_LATTICE_CERTIFICATE__NO_CLOSURE_THEOREM",
        )
        self.assertIn("target-signature-lattice", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "theorem_frontier_truth_table",
            "theorem_frontier_atom_influence",
            "theorem_frontier_cut",
        })
        self.assertEqual(payload["targets"], [
            "bridge_theorem_level_closure",
            "role_transfer_theorem_level_closure",
            "selector_qw2191_closure",
            "toe_closure",
        ])
        self.assertIn("target signature", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("target-lattice audit", payload["grep_disambiguation"]["finding"])

    def test_signature_lattice_summary_and_rows(self):
        payload = self.payload
        summary = payload["theorem_frontier_target_signature_lattice_summary"]
        self.assertEqual(summary["target_count"], 4)
        self.assertEqual(summary["all_target_signature_count"], 16)
        self.assertEqual(summary["reachable_target_signature_count"], 6)
        self.assertEqual(summary["unreachable_target_signature_count"], 10)
        self.assertEqual(summary["reachable_signatures"], ["0000", "0010", "0110", "1000", "1010", "1111"])
        self.assertEqual(summary["counts_by_reachable_signature"], {
            "0000": 56,
            "0010": 49,
            "0110": 7,
            "1000": 8,
            "1010": 7,
            "1111": 1,
        })
        self.assertEqual(summary["minimal_weights_by_reachable_signature"], {
            "0000": 0,
            "0010": 1,
            "0110": 4,
            "1000": 3,
            "1010": 4,
            "1111": 7,
        })
        self.assertTrue(summary["role_implies_selector"])
        self.assertTrue(summary["toe_implies_all_targets"])
        self.assertTrue(summary["only_full_signature_has_toe_closure"])
        self.assertTrue(summary["current_signature_all_false"])
        self.assertTrue(summary["atom_influence_top_atom_inherited"])
        self.assertFalse(summary["bridge_theorem_exported"])
        self.assertFalse(summary["role_transfer_theorem_exported"])
        self.assertFalse(summary["selector_closure_exported"])
        self.assertFalse(summary["toe_closure_claimed"])

        rows = {row["signature"]: row for row in payload["signature_rows"]}
        self.assertEqual(len(rows), 16)
        self.assertEqual(rows["0000"]["minimal_true_atom_count"], 0)
        self.assertEqual(rows["0010"]["minimal_true_atom_count"], 1)
        self.assertEqual(rows["1000"]["minimal_true_atom_count"], 3)
        self.assertEqual(rows["1111"]["minimal_true_atom_count"], 7)
        self.assertFalse(rows["0001"]["reachable"])
        self.assertFalse(rows["1110"]["reachable"])
        self.assertTrue(all(payload["cross_checks"].values()))

    def test_implications_proof_and_hard_limits(self):
        payload = self.payload
        implications = {row["implication"]: row for row in payload["implication_rows"]}
        self.assertTrue(implications["role_transfer_theorem_level_closure => selector_qw2191_closure"]["holds_on_all_assignments"])
        self.assertTrue(implications["toe_closure => bridge_theorem_level_closure"]["holds_on_all_assignments"])
        self.assertTrue(implications["toe_closure => role_transfer_theorem_level_closure"]["holds_on_all_assignments"])
        self.assertTrue(implications["toe_closure => selector_qw2191_closure"]["holds_on_all_assignments"])

        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("4-bit target signature", proof["projection_step"])
        self.assertIn("Exactly 6 of 16", proof["reachability_step"])
        self.assertIn("summing to 128", proof["count_step"])
        self.assertIn("role-transfer closure implies selector closure", proof["implication_step"])
        self.assertIn("all-open-atom frontier cut", proof["minimal_step"])
        self.assertIn("exports no theorem atom", proof["scope_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No reachable target signature", hard_limits)
        self.assertIn("No missing theorem atom", hard_limits)
        self.assertIn("No bridge theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
