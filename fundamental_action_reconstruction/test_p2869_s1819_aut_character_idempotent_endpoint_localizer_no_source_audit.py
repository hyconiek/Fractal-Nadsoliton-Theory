import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2868_SCRIPT = ROOT / "p2868_s1818_prime5_extended_aut_invariant_weighted_log_no_go_audit.py"
SCRIPT = ROOT / "p2869_s1819_aut_character_idempotent_endpoint_localizer_no_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2869_s1819_aut_character_idempotent_endpoint_localizer_no_source_audit.json"
MD_PATH = ROOT / "generated" / "p2869_s1819_aut_character_idempotent_endpoint_localizer_no_source_audit.md"


class P2869AutCharacterIdempotentEndpointLocalizerNoSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2868_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["aut_character_idempotent_endpoint_localizer_no_source_audit"]

    def test_status_and_p2868_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2869_AUT_CHARACTER_IDEMPOTENT_ENDPOINT_LOCALIZER_NO_SOURCE_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2868_PRIME5_EXTENDED_AUT_INVARIANT_WEIGHTED_LOG_NO_GO_AUDIT_NO_CLOSURE",
        )

    def test_character_idempotent_exactly_isolates_endpoint(self):
        weights = self.audit["endpoint_localizer_weights_on_unit_orbit"]
        self.assertEqual(weights["1"]["fraction"], "0/1")
        self.assertEqual(weights["5"]["fraction"], "0/1")
        self.assertEqual(weights["7"]["fraction"], "0/1")
        self.assertEqual(weights["11"]["fraction"], "1/1")

    def test_scaled_idempotent_exactly_matches_target_vector(self):
        self.assertTrue(self.audit["exact_target_representation"])
        vector = self.audit["scaled_prime_vector"]
        self.assertEqual(vector["5"]["fraction"], "0/1")
        self.assertEqual(vector["7"]["fraction"], "0/1")
        self.assertEqual(vector["11"]["fraction"], "9/5")

    def test_idempotent_remains_imported_non_source(self):
        requirements = self.audit["imported_requirements"]
        self.assertTrue(requirements["chosen_endpoint_label_d_11"])
        self.assertTrue(requirements["chosen_aut_character_polarity_projector"])
        self.assertTrue(requirements["denominator_prime_5_in_scaled_coefficients"])
        self.assertTrue(requirements["unit_bearing_coupling_theorem_missing"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["selector_or_localizer_source_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["eom_closure_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2869/S1819", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2869/S1819", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2869/S1819", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2869", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
