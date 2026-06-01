from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2390_s1340_selector_qualified_beta_tors_chi11_role_audit.py"
OUT = ROOT / "generated" / "p2390_s1340_selector_qualified_beta_tors_chi11_role_audit.json"
MD = ROOT / "generated" / "p2390_s1340_selector_qualified_beta_tors_chi11_role_audit.md"
COMPONENT_GAP_SCRIPT = ROOT / "scratch" / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_probe.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2390SelectorQualifiedBetaTorsChi11RoleAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(COMPONENT_GAP_SCRIPT)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["selector_qualified_beta_tors_chi11_role_audit"]
        cls.cert = cls.probe["implication_certificate"]
        cls.grep = cls.probe["grep_nonduplication_audit"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2390_s1340_v1")
        self.assertEqual(self.payload["packet_id"], "P2390")
        self.assertEqual(self.payload["stage_id"], "S1340")
        self.assertEqual(self.payload["result_kind"], "SELECTOR_QUALIFIED_BETA_TORS_CHI11_ROLE_AUDIT")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_grep_nonduplication_audit(self) -> None:
        self.assertGreater(self.grep["searched_file_count"], 100)
        self.assertGreater(self.grep["pattern_match_counts"]["selector_source"], 0)
        self.assertGreater(self.grep["pattern_match_counts"]["role_transfer"], 0)
        self.assertGreater(self.grep["pattern_match_counts"]["beta_tors_to_chi11_literal"], 0)
        self.assertIn("does not try to re-prove the selector", self.grep["nonduplication_finding"])

    def test_selector_present_but_transfer_not_licensed(self) -> None:
        bits = self.cert["audited_bits"]
        self.assertTrue(bits["strict_selector_exported_P1343"])
        self.assertTrue(bits["declared_scope_global_closure_P1348"])
        self.assertTrue(bits["component_gap_says_beta_tors_to_chi11_candidate_not_theorem"])
        self.assertFalse(bits["explicit_beta_tors_to_chi11_transport_theorem_found"])
        self.assertFalse(bits["component_bridge_ready_for_role_transfer"])
        self.assertFalse(bits["legacy_role_transfer_allowed_now"])
        self.assertFalse(self.cert["transfer_licensed"])
        self.assertEqual(self.cert["status"], "SELECTOR_PRESENT_BETA_TORS_CHI11_ROLE_TRANSFER_STILL_UNLICENSED")

    def test_boolean_implication_selector_not_sufficient(self) -> None:
        row = self.cert["current_state_row"]
        self.assertTrue(row["selector_source"])
        self.assertFalse(row["beta_tors_to_chi11_transport_theorem"])
        self.assertFalse(row["full_bridge_ready"])
        self.assertFalse(row["role_transfer_theorem"])
        self.assertFalse(row["licensed"])
        self.assertIn("explicit_beta_tors_to_chi11_transport_theorem", self.cert["single_missing_unlock_candidates"])
        self.assertTrue(
            any(
                truth_row["selector_source"]
                and not truth_row["beta_tors_to_chi11_transport_theorem"]
                and not truth_row["licensed"]
                for truth_row in self.cert["truth_table_slice_selector_not_sufficient"]
            )
        )

    def test_fingerprint_and_gatekeepers(self) -> None:
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("No beta_tors -> chi11 theorem is exported.", theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))


if __name__ == "__main__":
    unittest.main()
