import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2867_SCRIPT = ROOT / "p2867_s1817_coupled_z12_smooth_weighted_log_functional_no_source_audit.py"
SCRIPT = ROOT / "p2868_s1818_prime5_extended_aut_invariant_weighted_log_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2868_s1818_prime5_extended_aut_invariant_weighted_log_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2868_s1818_prime5_extended_aut_invariant_weighted_log_no_go_audit.md"


class P2868Prime5ExtendedAutInvariantWeightedLogNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2867_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["prime5_extended_aut_invariant_weighted_log_no_go_audit"]

    def test_status_and_p2867_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2868_PRIME5_EXTENDED_AUT_INVARIANT_WEIGHTED_LOG_NO_GO_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2867_COUPLED_Z12_SMOOTH_WEIGHTED_LOG_FUNCTIONAL_NO_SOURCE_AUDIT_NO_CLOSURE",
        )

    def test_prime5_is_granted_but_rank_system_inconsistent(self):
        self.assertIn(5, self.audit["denominator_primes_allowed"])
        system = self.audit["linear_system_certificate"]
        self.assertFalse(system["consistent"])
        self.assertLess(system["rank_A"], system["rank_augmented"])
        self.assertIn([1, 5, 7, 11], system["orbits"])

    def test_unit_orbit_vector_ties_5_7_11(self):
        system = self.audit["linear_system_certificate"]
        vector = system["orbit_prime_vectors"]["[1, 5, 7, 11]"]
        self.assertEqual(vector["5"]["fraction"], "1/1")
        self.assertEqual(vector["7"]["fraction"], "1/1")
        self.assertEqual(vector["11"]["fraction"], "1/1")

    def test_bounded_shared_unit_orbit_scan_has_no_exact_weight(self):
        scan = self.audit["unit_orbit_bounded_scan"]
        self.assertTrue(scan["exact_shared_weight_absent"])
        self.assertEqual(scan["exact_shared_unit_orbit_weight_hits"], [])
        self.assertGreater(scan["candidate_count"], 0)

    def test_no_false_exports_and_documents_updated(self):
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["eta_source_exported"])
        self.assertFalse(flags["selector_or_localizer_source_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2868/S1818", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2868/S1818", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2868/S1818", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2868", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
