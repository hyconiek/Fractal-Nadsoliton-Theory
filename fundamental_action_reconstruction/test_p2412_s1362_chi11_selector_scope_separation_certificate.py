#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2412_s1362_chi11_selector_scope_separation_certificate.py"
OUT = ROOT / "generated" / "p2412_s1362_chi11_selector_scope_separation_certificate.json"
MD = ROOT / "generated" / "p2412_s1362_chi11_selector_scope_separation_certificate.md"
P2411 = ROOT / "generated" / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json"


class P2412Chi11SelectorScopeSeparationCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2411.exists():
            subprocess.run([sys.executable, str(ROOT / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["chi11_selector_scope_separation_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_scope_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2412")
        self.assertEqual(self.payload["stage_id"], "S1362")
        self.assertEqual(self.payload["status"], "PASS_CHI11_SCOPE_SEPARATION_DECLARED_SELECTOR_NOT_BRIDGE_SOURCE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_scope_lattice(self) -> None:
        self.assertEqual(self.theorem["truth_row_count"], 32)
        self.assertEqual(self.theorem["current_mask"], 7)
        self.assertEqual(self.theorem["current_scope_separation_true_masks"], [7])
        self.assertTrue(self.theorem["current_scope_separation_true"])
        self.assertEqual(self.theorem["declared_selector_lane_true_count"], 16)
        self.assertEqual(self.theorem["bridge_selector_closure_lane_true_count"], 8)

    def test_current_source_facts_and_prime_implicant(self) -> None:
        facts = self.finite["current_source_facts"]
        self.assertTrue(facts["declared_scope_strict_selector_available"])
        self.assertTrue(facts["beta_tors_selector_route_retired"])
        self.assertTrue(facts["phase_origin_candidate_finite_audit_passed"])
        self.assertFalse(facts["bridge_chi11_source_theorem_exported"])
        self.assertFalse(facts["qw2191_discharge_exported"])
        prime = self.theorem["scope_separation_prime_implicant"]
        self.assertEqual(prime["literal_count"], 5)
        self.assertIn("bridge_chi11_source_theorem_exported", prime["required_false"])

    def test_candidate_classification_blocks_overread(self) -> None:
        rows = {row["candidate"]: row for row in self.finite["candidate_classification"]}
        self.assertTrue(rows["P1343_P1348_declared_scope_selector"]["licenses_declared_selector_lane"])
        self.assertFalse(rows["P1343_P1348_declared_scope_selector"]["licenses_bridge_chi11_source"])
        self.assertEqual(rows["beta_tors_to_chi11"]["current_status"], "retired_as_selector_search_route")
        self.assertEqual(rows["phase_origin_chiral_bispectrum_coprime_phase"]["current_status"], "finite_candidate_audited_not_strict_core_source")

    def test_inherited_open_bridge_limits(self) -> None:
        self.assertTrue(self.theorem["bridge_chi11_source_is_p2411_top_priority"])
        self.assertTrue(self.theorem["p2366_still_blocks_qw2191_discharge"])
        self.assertTrue(self.theorem["p2411_still_blocks_chi11_source"])
        self.assertIn("not a bridge-level chi11 source theorem", "\n".join(self.theorem["not_licensed"]))

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("scope-separation certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2412/S1362 chi11", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2412/S1362 chi11", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
