import json
import math
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2790_s1740_eighth_16node_witness_no_exhaustion_certificate.json"
MD_PATH = GEN / "p2790_s1740_eighth_16node_witness_no_exhaustion_certificate.md"


class P2790Eighth16NodeWitnessNoExhaustionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2790_s1740_eighth_16node_witness_no_exhaustion_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["eighth_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2790_EIGHTH_16NODE_WITNESS_NO_EXHAUSTION_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2786"], "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2789"], "P2789_ORBIT_STABILIZER_EXACT_QUOTIENT_CERTIFICATE_NO_CLOSURE")

    def test_eighth_witness_is_valid_16node_4regular_connected_graph(self):
        self.assertEqual(self.witness["edge_count"], 32)
        self.assertEqual(self.witness["degree_sequence"], [4] * 16)
        self.assertTrue(self.witness["is_connected"])
        self.assertEqual(self.witness["automorphism_group_size"], 2)
        self.assertEqual(math.factorial(16) // 2, self.witness["orbit_size_by_orbit_stabilizer"])
        self.assertEqual(len(self.witness["exact_charpoly_coefficients"]["adjacency"]), 17)
        self.assertEqual(len(self.witness["exact_charpoly_coefficients"]["laplacian"]), 17)
        self.assertEqual(len(self.witness["exact_charpoly_coefficients"]["signless_laplacian"]), 17)

    def test_eighth_witness_is_nonisomorphic_to_all_seven_local_representatives(self):
        self.assertEqual(len(self.witness["local_pair_isomorphism_checks"]), 7)
        self.assertTrue(self.witness["is_nonisomorphic_to_all_seven_local_representatives"])
        for row in self.witness["local_pair_isomorphism_checks"]:
            self.assertFalse(row["isomorphic_to_eighth_witness"])

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_eighth_witness_no_exhaustion_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2790", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2790/S1740", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2790/S1740", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2790/S1740", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2790", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
