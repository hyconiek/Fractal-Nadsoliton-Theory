import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_adjacency_unit_label_source_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_adjacency_unit_label_source_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_adjacency_unit_label_source_audit_report.md"


class StrictAlphaCyclicAdjacencyUnitLabelSourceAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_CYCLIC_ADJACENCY_UNIT_LABEL_SOURCE_AUDIT_PROBE__CONDITIONAL_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "cyclic-adjacency-conditionally-selects-chi_11-kernel-but-source-premise-not-derived",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_character_table_and_adjacency_kernel(self):
        characters = {row["name"]: row for row in self.payload["character_table"]}
        self.assertEqual(len(characters), 4)
        required = characters["chi_11_required_d5_unit_axis"]
        self.assertEqual(required["kernel"], [1, 11])
        self.assertEqual(required["nonzero_coset"], [5, 7])

        audit = self.payload["cyclic_adjacency_source_audit"]
        self.assertEqual(audit["unit_image_of_generator_plus_one"], {"1": 1, "5": 5, "7": 7, "11": 11})
        self.assertEqual(audit["unit_image_folded_shell"], {"1": 1, "5": 5, "7": 5, "11": 1})
        self.assertEqual(audit["nearest_neighbor_shell_1_preserving_units"], [1, 11])
        self.assertEqual(audit["nearest_neighbor_shell_1_preserving_character"], "chi_11_required_d5_unit_axis")
        self.assertTrue(audit["exact_positive_if_adjacency_admitted"])
        self.assertIn("abstract V4 meta-automorphisms", audit["why_this_breaks_meta_symmetry"])

    def test_shell_rows_and_selection_source_audit(self):
        shell_rows = {row["cyclic_shell"]: row for row in self.payload["cyclic_adjacency_source_audit"]["shell_rows"]}
        self.assertEqual(shell_rows[1]["preserving_units"], [1, 11])
        self.assertTrue(shell_rows[1]["is_required_kernel_{1,11}"])
        self.assertEqual(shell_rows[1]["selected_character"], "chi_11_required_d5_unit_axis")
        self.assertEqual(shell_rows[5]["preserving_units"], [5, 7])
        self.assertFalse(shell_rows[5]["is_required_kernel_{1,11}"])
        self.assertIn("no_Z2_character", shell_rows[5]["selected_character"])

        by_source = {row["candidate_source"]: row for row in self.payload["selection_source_audit"]}
        self.assertFalse(by_source["abstract_unit_group_V4_only"]["selects_required_chi_11"])
        self.assertTrue(by_source["cyclic_adjacency_shell_{+1,-1}"]["selects_required_chi_11"])
        self.assertEqual(by_source["cyclic_adjacency_shell_{+1,-1}"]["selects_kernel"], [1, 11])
        self.assertIn("if", by_source["cyclic_adjacency_shell_{+1,-1}"]["strict_status"])
        self.assertFalse(by_source["choose_long_step_shell_{+5,-5}"]["selects_required_chi_11"])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("V4 alone does not distinguish 11", proof["finite_group"])
        self.assertIn("{+1,-1}", proof["cyclic_adjacency_datum"])
        self.assertIn("u=1 and u=11", proof["stabilizer_computation"])
        self.assertIn("kernel of chi_11", proof["character_identification"])
        self.assertIn("conditional", proof["conditionality"])

        interpretation = self.payload["interpretation"]
        self.assertIn("would name the {1,11} kernel", interpretation["direct_result"])
        self.assertIn("If strict nadsoliton geometry exports", interpretation["positive_conditional"])
        self.assertIn("extra structure", interpretation["negative_result"])
        self.assertIn("QW-2191 is not discharged", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives cyclic adjacency", hard_limits)
        self.assertIn("conditional source audit", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
