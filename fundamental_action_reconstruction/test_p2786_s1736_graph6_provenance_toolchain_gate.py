import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
JSON = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.json"
MD = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.md"


class P2786Graph6ProvenanceToolchainGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON.exists():
            import p2786_s1736_graph6_provenance_toolchain_gate as p2786
            p2786.main()
        cls.payload = json.loads(JSON.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2785"], "P2785_EXACT_CHARACTERISTIC_POLYNOMIAL_CERTIFICATE_NO_CLOSURE")
        self.assertIn("graph6 provenance", self.payload["audited_question"])

    def test_graph6_provenance_counts(self):
        witness = self.payload["graph6_provenance_witness"]
        self.assertEqual(witness["representative_count"], 7)
        self.assertTrue(witness["all_graph6_roundtrips_exact"])
        self.assertTrue(witness["all_rows_have_16node_4regular_edge_count"])
        self.assertEqual(witness["distinct_graph6_hash_count"], 7)
        self.assertEqual(len(witness["provenance_rows"]), 7)
        for row in witness["provenance_rows"]:
            self.assertTrue(row["roundtrip_edges_equal"])
            self.assertEqual(row["edge_count"], 32)
            self.assertEqual(len(row["graph6_sha256"]), 64)
            self.assertEqual(len(row["exact_payload_sha256"]), 64)
            self.assertTrue(row["graph6"].startswith("O"))

    def test_toolchain_gate_blocks_full_generator_claim(self):
        witness = self.payload["graph6_provenance_witness"]
        self.assertFalse(witness["canonical_generation_tool_available"])
        self.assertFalse(witness["common_python_graph_generator_available"])
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["accepted_as_local_graph6_provenance_certificate"])
        self.assertFalse(acceptance["accepted_as_full_canonical_generation_tool_certificate"])
        self.assertFalse(acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2786", self.payload["decision"]["next_honest_step"])

    def test_documentation_updated(self):
        self.assertIn("P2786/S1736", MD.read_text(encoding="utf-8"))
        self.assertIn("P2786/S1736", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2786/S1736", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2786", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
