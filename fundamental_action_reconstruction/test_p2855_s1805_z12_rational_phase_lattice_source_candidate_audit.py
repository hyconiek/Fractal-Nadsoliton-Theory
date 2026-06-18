import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2853_SCRIPT = ROOT / "p2853_s1803_phase_frequency_bridge_source_audit.py"
P2854_SCRIPT = ROOT / "p2854_s1804_post_p2853_professorial_state_map_no_new_live_frontier.py"
SCRIPT = ROOT / "p2855_s1805_z12_rational_phase_lattice_source_candidate_audit.py"
JSON_PATH = ROOT / "generated" / "p2855_s1805_z12_rational_phase_lattice_source_candidate_audit.json"
MD_PATH = ROOT / "generated" / "p2855_s1805_z12_rational_phase_lattice_source_candidate_audit.md"


class P2855Z12RationalPhaseLatticeSourceCandidateAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2853_SCRIPT, P2854_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["z12_rational_phase_lattice_source_candidate_audit"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2855_Z12_RATIONAL_PHASE_LATTICE_SOURCE_CANDIDATE_AUDIT_NO_CLOSURE")
        statuses = self.audit["input_statuses_rechecked"]
        self.assertEqual(statuses["P2853"], "P2853_PHASE_FREQUENCY_BRIDGE_SOURCE_AUDIT_NO_CLOSURE")
        self.assertEqual(statuses["P2854"], "P2854_POST_P2853_PROFESSORIAL_STATE_MAP_NO_NEW_LIVE_FRONTIER")

    def test_prime_five_denominator_obstruction(self):
        omega = self.audit["strict_omega"]
        phi = self.audit["strict_phi"]
        self.assertEqual(omega["fraction"], "743/4000")
        self.assertEqual(phi["fraction"], "13/80")
        self.assertIn("5", {str(key) for key in omega["denominator_prime_factors"]})
        self.assertIn("5", {str(key) for key in phi["denominator_prime_factors"]})
        self.assertFalse(omega["z12_compatible_denominator"])
        self.assertFalse(phi["z12_compatible_denominator"])

    def test_z12_approximations_are_nonexact(self):
        omega_approx = self.audit["omega_best_z12_approximation"]
        phi_approx = self.audit["phi_best_z12_approximation"]
        common = self.audit["best_common_z12_pair"]
        self.assertFalse(omega_approx["exact_match_found"])
        self.assertFalse(phi_approx["exact_match_found"])
        self.assertFalse(common["exact_pair_match_found"])
        self.assertGreater(omega_approx["absolute_error_float"], 0.0)
        self.assertGreater(phi_approx["absolute_error_float"], 0.0)

    def test_candidate_matrix_no_source_export(self):
        rows = {row["candidate"]: row for row in self.audit["candidate_matrix"]}
        self.assertFalse(rows["pure_z12_denominator_lattice_exact_source"]["finite_witness_passes"])
        self.assertFalse(rows["pure_z12_denominator_lattice_exact_source"]["exports_strict_source_law"])
        self.assertTrue(rows["bounded_z12_lattice_best_approximation"]["finite_witness_passes"])
        self.assertFalse(rows["bounded_z12_lattice_best_approximation"]["exports_strict_source_law"])
        self.assertTrue(rows["import_prime5_phase_unit_extension"]["finite_witness_passes"])
        self.assertFalse(rows["import_prime5_phase_unit_extension"]["exports_strict_source_law"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)

    def test_phase_bits_and_no_false_closure_documents(self):
        self.assertEqual(self.audit["exact_strict_phase_bits"], self.audit["p2853_phase_bits"])
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_z12_rational_phase_lattice_obstruction_audit"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_strict_phase_frequency_source_law"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["z12_phase_lattice_source_exported"])
        self.assertFalse(flags["prime5_phase_unit_source_exported"])
        self.assertFalse(flags["strict_phase_frequency_source_law_exported"])
        self.assertFalse(flags["full_kernel_bridge_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2855/S1805", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2855/S1805", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2855/S1805", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2855", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
