import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2961_s1911_exchange_balanced_scale_quotient_source_candidate.py"
OUT = ROOT / "generated" / "p2961_s1911_exchange_balanced_scale_quotient_source_candidate.json"
MD = ROOT / "generated" / "p2961_s1911_exchange_balanced_scale_quotient_source_candidate.md"

class P2961ExchangeBalancedScaleQuotientSourceCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2961_EXCHANGE_BALANCED_SCALE_QUOTIENT_SOURCE_CANDIDATE_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2960"])

    def test_unique_fixed_quotient_target(self):
        cert = self.payload["quotient_certificate"]
        self.assertEqual(cert["bounded_pair_rows"], 144)
        self.assertTrue(cert["unique_fixed_orbit_is_target"])
        self.assertEqual(cert["target_sum_from_fixed_orbit"], 9)
        self.assertTrue(cert["developmental_source_candidate_exported"])
        self.assertFalse(cert["strict_nadsoliton_exchange_symmetry_source_exported"])
        self.assertFalse(cert["strict_ratio_package_source_exported"])

    def test_obligations_and_matrix(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertIn("gcd(a,b)=1", objs["unbounded_lemma"])
        obligations = {row["obligation"]: row["satisfied"] for row in objs["proof_obligation_rows"]}
        self.assertTrue(obligations["canonical_gcd_scale_quotient_constructed"])
        self.assertTrue(obligations["exchange_action_computed_on_quotient"])
        self.assertTrue(obligations["unique_exchange_fixed_primitive_orbit"])
        self.assertTrue(obligations["target_vector_emerges_from_fixed_orbit"])
        self.assertFalse(obligations["strict_nadsoliton_exchange_symmetry_source_exported"])
        self.assertFalse(obligations["unit_bearing_nonproxy_coupling_exported"])
        matrix = objs["finite_acceptance_matrix"]
        self.assertEqual(len(matrix), 64)
        self.assertEqual(sum(1 for row in matrix if row["accepts_as_developmental_source_theorem_candidate"]), 4)
        self.assertEqual(sum(1 for row in matrix if row["accepts_as_strict_ratio_package_source"]), 1)

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2961/S1911", MD.read_text(encoding="utf-8"))
        self.assertIn("P2961/S1911", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2961/S1911", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2961", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
