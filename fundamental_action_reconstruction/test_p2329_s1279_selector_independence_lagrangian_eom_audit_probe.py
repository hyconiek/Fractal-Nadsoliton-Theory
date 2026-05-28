from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2329_s1279_selector_independence_lagrangian_eom_audit_probe.json"
MD = ROOT / "generated" / "p2329_s1279_selector_independence_lagrangian_eom_audit_probe.md"
SCRIPT = ROOT / "p2329_s1279_selector_independence_lagrangian_eom_audit_probe.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2329SelectorIndependenceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_selector_independence_lagrangian_eom_audit_probe"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2329_s1279_v1")
        self.assertEqual(self.payload["packet_id"], "P2329")
        self.assertEqual(self.payload["stage_id"], "S1279")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_selector_independent_terms_variations_and_residuals(self) -> None:
        cert = self.probe["selector_independence_certificate"]
        self.assertTrue(cert["p2086_loaded"])
        self.assertEqual(cert["term_count"], 11)
        self.assertEqual(cert["selector_independent_term_count"], 11)
        self.assertEqual(cert["variation_field_count"], 3)
        self.assertEqual(cert["selector_independent_variation_field_count"], 3)
        self.assertTrue(cert["all_recomposition_residuals_selector_independent"])
        self.assertEqual(
            sorted(row["field"] for row in self.probe["selector_independent_variations"]),
            ["A", "h", "psi"],
        )
        self.assertTrue(all(not row["requires_S_future_for_termwise_variation"] for row in self.probe["selector_independent_terms"]))
        self.assertTrue(all(not row["requires_S_future_for_reduced_variation"] for row in self.probe["selector_independent_variations"]))

    def test_observable_claim_boundary_and_gatekeepers(self) -> None:
        self.assertGreaterEqual(len(self.probe["safe_observable_claims"]), 4)
        self.assertGreaterEqual(len(self.probe["blocked_selector_dependent_claims"]), 3)
        blocked_ids = {row["claim_id"] for row in self.probe["blocked_selector_dependent_claims"]}
        self.assertIn("exact_future_branch_actualization", blocked_ids)
        self.assertIn("provider_lift_per_step_policy_margin", blocked_ids)
        checks = self.payload["gatekeeper_checks"]
        for key in [
            "all_terms_selector_independent",
            "all_variations_selector_independent",
            "all_residuals_selector_independent",
            "p2325_no_admissible_signed_source_preserved",
            "no_qw2191_discharge_claimed",
            "no_g1_g3_update_claimed",
            "no_toe_closure_claimed",
        ]:
            self.assertTrue(checks[key], key)

    def test_theorem_fingerprint(self) -> None:
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("selector-independent", theorem["claim"])
        self.assertTrue(theorem["strict_guardrails"]["no_selector_premise_added"])


if __name__ == "__main__":
    unittest.main()
