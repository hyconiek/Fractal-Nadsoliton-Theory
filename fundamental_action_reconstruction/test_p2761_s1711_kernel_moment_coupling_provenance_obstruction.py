import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2761_s1711_kernel_moment_coupling_provenance_obstruction.py"
OUT = ROOT / "generated" / "p2761_s1711_kernel_moment_coupling_provenance_obstruction.json"
MD = ROOT / "generated" / "p2761_s1711_kernel_moment_coupling_provenance_obstruction.md"


class P2761KernelMomentCouplingProvenanceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=REPO, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2761_KERNEL_MOMENT_COUPLING_PROVENANCE_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P1562"], "PASS_STRICT_KERNEL_COEFFICIENTS_DERIVED")
        self.assertEqual(self.payload["input_statuses"]["P1563"], "PASS_STRICT_KERNEL_TO_EOM_CHAIN_EXPORTED")
        self.assertEqual(self.payload["input_statuses"]["P1866"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(self.payload["input_statuses"]["P2760"], "P2760_FOUNDATION_KERNEL_LAGRANGIAN_GAP_MATRIX_NO_CLOSURE")

    def test_content_scan_and_rows(self):
        self.assertTrue(self.payload["content_evidence_scan"]["all_patterns_have_hits"])
        matrix = self.payload["dimension_obligation_matrix"]
        self.assertEqual(matrix["row_count"], 3)
        self.assertEqual(matrix["accepted_physical_coupling_count"], 0)
        names = {row["coupling"] for row in matrix["rows"]}
        self.assertEqual(names, {"lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"})

    def test_missing_global_provenance_and_stale_flags(self):
        matrix = self.payload["dimension_obligation_matrix"]
        self.assertEqual(matrix["missing_global_provenance_atom_count"], 5)
        self.assertIn("canonical physical length/reference cell mapping for moment powers", matrix["missing_global_provenance_atoms"])
        self.assertTrue(matrix["p1562_stale_closure_flag_detected"])
        self.assertTrue(matrix["later_artifacts_block_closure"])
        self.assertTrue(matrix["stale_flags_quarantined"])

    def test_acceptance_blocks_g5_theorem(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertFalse(acceptance["accepted_as_g5_physical_coupling_provenance_theorem"])
        self.assertIn("unit_reference_source_exported", acceptance["missing_criteria"])
        self.assertIn("sign_convention_theorem_exported", acceptance["missing_criteria"])
        self.assertIn("variational_normalization_exported", acceptance["missing_criteria"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertTrue(all(value is False for value in flags.values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("canonical physical length/reference-cell", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2761", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2761/S1711", MD.read_text(encoding="utf-8"))
        self.assertIn("P2761/S1711", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2761/S1711", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2761", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
