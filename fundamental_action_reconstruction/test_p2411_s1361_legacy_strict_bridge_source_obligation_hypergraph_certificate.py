#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.py"
OUT = ROOT / "generated" / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json"
MD = ROOT / "generated" / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.md"
P2410 = ROOT / "generated" / "p2410_s1360_dequotiented_twelve_atom_prime_implicate_obstruction_certificate.json"


class P2411LegacyStrictBridgeSourceObligationHypergraphCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2410.exists():
            subprocess.run([sys.executable, str(ROOT / "p2410_s1360_dequotiented_twelve_atom_prime_implicate_obstruction_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["legacy_strict_bridge_source_obligation_hypergraph_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_hypergraph_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2411")
        self.assertEqual(self.payload["stage_id"], "S1361")
        self.assertEqual(self.payload["status"], "PASS_BRIDGE_SOURCE_OBLIGATION_HYPERGRAPH_FULL_MASK_ONLY_NO_ROLE_TRANSFER")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_bridge_lattice_counts(self) -> None:
        self.assertEqual(self.theorem["bridge_obligation_atom_count"], 8)
        self.assertEqual(self.theorem["bridge_component_count"], 5)
        self.assertEqual(self.theorem["total_bridge_obligation_assignment_count"], 256)
        self.assertEqual(self.theorem["bridge_ready_true_masks"], [255])
        self.assertEqual(self.theorem["bridge_ready_anf_degree"], 8)
        self.assertEqual(self.theorem["proper_subset_failure_count"], 255)

    def test_current_empty_mask_and_nearest_misses(self) -> None:
        self.assertEqual(self.theorem["empty_mask_ready_components"], ["residual_strict_additions_inventory"])
        self.assertFalse(self.theorem["empty_mask_bridge_ready"])
        self.assertEqual(self.theorem["nearest_miss_count"], 8)
        nearest = self.finite["nearest_miss_rows"]
        self.assertEqual(len(nearest), 8)
        self.assertTrue(all(len(row["blocked_components_before_repair"]) >= 1 for row in nearest))

    def test_component_hypergraph_and_priority(self) -> None:
        components = self.finite["bridge_components"]
        self.assertIn("amplitude_normalization_passage", components)
        self.assertIn("selector_source_premise", components)
        self.assertEqual(components["residual_strict_additions_inventory"], [])
        top_atoms = [row["atom"] for row in self.theorem["top_priority_atoms"]]
        self.assertIn("chi11_selector_source_theorem", top_atoms)
        chi11 = next(row for row in self.finite["atom_component_incidence"] if row["atom"] == "chi11_selector_source_theorem")
        self.assertEqual(chi11["component_count"], 2)

    def test_inherited_limits(self) -> None:
        self.assertTrue(self.theorem["p2410_no_ltotal_or_toe_closure_inherited"])
        self.assertTrue(self.theorem["scratch_component_gap_sources_missing_inherited"])
        self.assertTrue(self.theorem["scratch_role_transfer_blocked_inherited"])
        self.assertIn("No chi11 selector-source theorem", "\n".join(self.theorem["not_licensed"]))

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("source-obligation hypergraph", MD.read_text(encoding="utf-8"))
        self.assertIn("P2411/S1361 legacy-to-strict", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2411/S1361 bridge-source", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
