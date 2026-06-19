import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2937_s1887_z12_multiplication_kernel_prime_coordinate_source_candidate.py"
OUT = ROOT / "generated" / "p2937_s1887_z12_multiplication_kernel_prime_coordinate_source_candidate.json"
MD = ROOT / "generated" / "p2937_s1887_z12_multiplication_kernel_prime_coordinate_source_candidate.md"


class P2937Z12MultiplicationKernelPrimeCoordinateSourceCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2937_Z12_MULTIPLICATION_KERNEL_PRIME_COORDINATE_SOURCE_CANDIDATE_REJECTED_PARTIAL")
        self.assertIsNotNone(self.payload["input_hashes"]["P2936"])

    def test_candidate_certificate(self):
        cert = self.payload["candidate_certificate"]
        self.assertEqual(cert["prime_coordinate_vector_order_2_3_5_7_11"], [1, 2, 0, 0, 0])
        self.assertEqual(cert["nonzero_prime_coordinate_count"], 2)
        self.assertEqual(cert["zero_prime_coordinates"], [5, 7, 11])
        self.assertEqual(cert["product_pair_count_de_le_11"], 29)
        self.assertEqual(cert["product_additivity_defect_count"], 0)
        self.assertFalse(cert["accepted_strict_prime_log_source"])

    def test_prime_rows_and_acceptance_boundary(self):
        rows = self.payload["constructed_theoretical_objects"]["prime_kernel_rows"]
        by_prime = {row["prime"]: row for row in rows}
        self.assertEqual(by_prime[2]["kernel_size"], 2)
        self.assertEqual(by_prime[3]["kernel_size"], 3)
        self.assertEqual(by_prime[5]["kernel_size"], 1)
        self.assertTrue(all(row["kernel_times_image_equals_12"] for row in rows))
        acceptance = self.payload["constructed_theoretical_objects"]["acceptance_rows"]
        self.assertFalse(next(row for row in acceptance if row["criterion"] == "all_five_prime_coordinates_nonzero")["satisfied"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_aut_breaking_prime_coordinate_source_law_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2937/S1887", MD.read_text(encoding="utf-8"))
        self.assertIn("P2937/S1887", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2937/S1887", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2937", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
