import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2866_SCRIPT = ROOT / "p2866_s1816_product_theorem_candidate_pair_no_source_audit.py"
SCRIPT = ROOT / "p2867_s1817_coupled_z12_smooth_weighted_log_functional_no_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2867_s1817_coupled_z12_smooth_weighted_log_functional_no_source_audit.json"
MD_PATH = ROOT / "generated" / "p2867_s1817_coupled_z12_smooth_weighted_log_functional_no_source_audit.md"


class P2867CoupledZ12SmoothWeightedLogFunctionalNoSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2866_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["coupled_z12_smooth_weighted_log_functional_no_source_audit"]

    def test_status_and_p2866_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2867_COUPLED_Z12_SMOOTH_WEIGHTED_LOG_FUNCTIONAL_NO_SOURCE_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2866_PRODUCT_THEOREM_CANDIDATE_PAIR_NO_SOURCE_AUDIT_NO_CLOSURE",
        )

    def test_prime11_coordinate_forces_imported_9_over_5(self):
        coordinate = self.audit["coordinate_forcing_certificate"]
        self.assertEqual(coordinate["prime11_contributors"], [11])
        self.assertEqual(coordinate["forced_w_11"]["fraction"], "9/5")
        self.assertEqual(coordinate["forced_w_11_denominator_primes"], [5])
        self.assertFalse(coordinate["forced_w_11_is_z12_smooth"])

    def test_bounded_z12_smooth_scan_has_no_exact_hit(self):
        scan = self.audit["bounded_w11_scan"]
        self.assertTrue(scan["exact_w11_absent"])
        self.assertEqual(scan["exact_w11_hits"], [])
        self.assertGreater(scan["candidate_count"], 0)

    def test_imported_witness_exact_but_not_sourced(self):
        witness = self.audit["imported_exact_witness"]
        self.assertEqual(witness["weights"]["11"]["fraction"], "9/5")
        self.assertEqual(witness["weights"]["11"]["denominator_prime_factors"], {"5": 1})
        self.assertTrue(witness["exact_representation"])
        self.assertFalse(witness["exports_sourcehood"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["eta_source_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2867/S1817", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2867/S1817", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2867/S1817", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2867", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
