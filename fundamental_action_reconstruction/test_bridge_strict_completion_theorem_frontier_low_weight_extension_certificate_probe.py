import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_report.md"


class TheoremFrontierLowWeightExtensionCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_scope(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_THEOREM_FRONTIER_LOW_WEIGHT_EXTENSION_CERTIFICATE__NO_THEOREM_EXPORT",
        )
        self.assertIn("singleton-and-pair-frontier-extensions", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "theorem_frontier_truth_table",
            "theorem_frontier_target_signature_lattice",
            "theorem_frontier_atom_influence",
        })
        self.assertIn("singleton extension", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("low-weight extension audit", payload["grep_disambiguation"]["finding"])

    def test_singleton_pair_summary_and_rows(self):
        payload = self.payload
        summary = payload["theorem_frontier_low_weight_extension_summary"]
        self.assertEqual(summary["open_atom_count"], 7)
        self.assertEqual(summary["singleton_extension_count"], 7)
        self.assertEqual(summary["pair_extension_count"], 21)
        self.assertEqual(summary["low_weight_extension_count"], 28)
        self.assertEqual(summary["singleton_unlock_count"], 1)
        self.assertEqual(summary["singleton_unlock_atoms"], ["chi11_selector_source"])
        self.assertEqual(summary["singleton_unlock_signatures"], ["0010"])
        self.assertTrue(summary["chi11_is_only_singleton_unlock"])
        self.assertEqual(summary["pair_unlock_count"], 6)
        self.assertTrue(summary["pair_unlocks_are_selector_only"])
        self.assertTrue(summary["pair_unlocks_all_contain_chi11"])
        self.assertTrue(summary["no_singleton_closes_bridge_role_or_toe"])
        self.assertTrue(summary["no_pair_closes_bridge"])
        self.assertTrue(summary["no_pair_closes_role_transfer"])
        self.assertTrue(summary["no_pair_closes_toe"])
        self.assertTrue(summary["target_lattice_min_weights_inherited"])
        self.assertTrue(summary["atom_influence_top_atom_inherited"])
        self.assertTrue(summary["current_signature_all_false"])
        self.assertFalse(summary["bridge_theorem_exported"])
        self.assertFalse(summary["role_transfer_theorem_exported"])
        self.assertFalse(summary["selector_closure_exported"])
        self.assertFalse(summary["toe_closure_claimed"])

        singleton_rows = payload["singleton_extension_rows"]
        self.assertEqual(len(singleton_rows), 7)
        singleton_unlocks = [row for row in singleton_rows if row["newly_closed_targets"]]
        self.assertEqual(len(singleton_unlocks), 1)
        self.assertEqual(singleton_unlocks[0]["true_atoms"], ["chi11_selector_source"])
        self.assertEqual(singleton_unlocks[0]["target_signature"], "0010")
        pair_rows = payload["pair_extension_rows"]
        self.assertEqual(len(pair_rows), 21)
        pair_unlocks = [row for row in pair_rows if row["newly_closed_targets"]]
        self.assertEqual(len(pair_unlocks), 6)
        self.assertTrue(all("chi11_selector_source" in row["true_atoms"] for row in pair_unlocks))
        self.assertTrue(all(row["target_signature"] == "0010" for row in pair_unlocks))
        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("All 7 singleton and all 21 pair", proof["enumeration_step"])
        self.assertIn("only singleton extension", proof["singleton_step"])
        self.assertIn("Exactly six pair extensions", proof["pair_step"])
        self.assertIn("No singleton or pair extension closes bridge", proof["blocker_step"])
        self.assertIn("exports no theorem atom", proof["scope_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No singleton or pair extension", hard_limits)
        self.assertIn("No missing theorem atom", hard_limits)
        self.assertIn("No bridge theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
