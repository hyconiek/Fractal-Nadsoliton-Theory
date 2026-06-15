import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2766_s1716_post_moment_provenance_state_reconciliation.py"
OUT = ROOT / "generated" / "p2766_s1716_post_moment_provenance_state_reconciliation.json"
MD = ROOT / "generated" / "p2766_s1716_post_moment_provenance_state_reconciliation.md"


class P2766PostMomentProvenanceStateReconciliationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2766_POST_MOMENT_PROVENANCE_STATE_RECONCILIATION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2761"], "P2761_KERNEL_MOMENT_COUPLING_PROVENANCE_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2762"], "P2762_REFERENCE_CELL_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2763"], "P2763_MOMENT_COUPLING_SIGN_CONVENTION_CONDITIONAL_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2764"], "P2764_FIELD_CURVATURE_NORMALIZATION_COMPATIBILITY_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2765"], "P2765_NONPROXY_VARIATIONAL_INSERTION_RESIDUAL_OBSTRUCTION_NO_CLOSURE")

    def test_content_scan_and_atom_matrix(self):
        self.assertTrue(self.payload["content_evidence_scan"]["all_patterns_have_hits"])
        matrix = self.payload["provenance_atom_matrix"]
        self.assertEqual(matrix["row_count"], 4)
        self.assertEqual(matrix["closed_atom_count"], 0)
        self.assertEqual(matrix["open_atom_count"], 4)
        self.assertEqual(set(matrix["open_atoms"]), {"reference_cell_action_density", "sign_convention", "field_curvature_normalization", "nonproxy_variational_insertion"})
        self.assertFalse(matrix["all_named_atoms_closed"])

    def test_finite_closure_lattice_blocks_current_profile(self):
        lattice = self.payload["finite_closure_lattice"]
        self.assertEqual(lattice["profile_count"], 16)
        self.assertEqual(lattice["license_profile_count"], 1)
        self.assertEqual(lattice["current_profile"]["closed_count"], 0)
        self.assertFalse(lattice["current_profile_is_license_profile"])
        self.assertFalse(lattice["current_profile"]["would_license_physical_coupling_provenance"])

    def test_acceptance_blocks_ltotal_promotion(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertFalse(acceptance["accepted_as_physical_coupling_provenance_theorem"])
        self.assertFalse(acceptance["accepted_as_ltotal_promotion"])
        self.assertTrue(acceptance["facts"]["no_named_atom_closed"])
        self.assertTrue(acceptance["facts"]["current_profile_does_not_license_physical_provenance"])
        self.assertIn("new_theorem_supplied_in_this_reconciliation", acceptance["missing_criteria"])
        self.assertIn("physical_coupling_provenance_theorem_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("P2697-P2766", self.payload["decision"]["next_honest_step"])
        self.assertIn("fresh broad state-map", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2766/S1716", MD.read_text(encoding="utf-8"))
        self.assertIn("P2766/S1716", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2766/S1716", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2766", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
