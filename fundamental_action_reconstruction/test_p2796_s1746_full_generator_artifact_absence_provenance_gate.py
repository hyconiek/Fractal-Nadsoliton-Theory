import json
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2796_s1746_full_generator_artifact_absence_provenance_gate.json"
MD_PATH = GEN / "p2796_s1746_full_generator_artifact_absence_provenance_gate.md"


class P2796FullGeneratorArtifactAbsenceProvenanceGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2796_s1746_full_generator_artifact_absence_provenance_gate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["full_generator_artifact_absence_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2796_FULL_GENERATOR_ARTIFACT_ABSENCE_PROVENANCE_GATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2795"], "P2795_POST_SUBCLASS_SATURATION_NO_NEW_CLASS_FRONTIER_CERTIFICATE_NO_CLOSURE")

    def test_repository_scan_finds_no_true_full_generator_certificate(self):
        self.assertGreater(self.witness["generated_json_file_count"], 0)
        self.assertGreater(self.witness["status_row_count"], 0)
        self.assertGreater(self.witness["status_rows_with_no_closure"], 0)
        self.assertGreater(self.witness["graph6_payload_row_count"], 0)
        self.assertEqual(self.witness["full_generator_truth_rows"], [])
        self.assertEqual(self.witness["full_generator_truth_row_count"], 0)
        self.assertGreater(self.witness["false_generator_guard_row_count"], 0)
        self.assertFalse(self.witness["required_artifact_present"])

    def test_external_generator_probe_is_recorded(self):
        commands = [row["command"] for row in self.witness["external_generator_command_rows"]]
        self.assertEqual(commands, ["genreg", "geng", "shortg", "labelg", "dreadnaut", "plantri"])
        self.assertEqual(self.witness["external_generator_command_available_count"], sum(1 for row in self.witness["external_generator_command_rows"] if row["available"]))

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_full_generator_absence_provenance_gate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2796", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2796/S1746", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2796/S1746", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2796/S1746", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2796", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
