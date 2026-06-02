#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2419_s1369_chi11_phase_selector_coupling_cut_certificate.py"
OUT = ROOT / "generated" / "p2419_s1369_chi11_phase_selector_coupling_cut_certificate.json"
MD = ROOT / "generated" / "p2419_s1369_chi11_phase_selector_coupling_cut_certificate.md"
P2418 = ROOT / "generated" / "p2418_s1368_bridge_source_marginal_unlock_lattice_certificate.json"


class P2419Chi11PhaseSelectorCouplingCutCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2418.exists():
            subprocess.run([sys.executable, str(ROOT / "p2418_s1368_bridge_source_marginal_unlock_lattice_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["chi11_phase_selector_coupling_cut_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_witness_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2419")
        self.assertEqual(self.payload["stage_id"], "S1369")
        self.assertEqual(self.payload["status"], "PASS_CHI11_PHASE_SELECTOR_COUPLING_CUT_NO_SOURCE_OR_SELECTOR_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_phase_selector_cut_counts(self) -> None:
        self.assertEqual(self.theorem["source_obligation_atom_count"], 8)
        self.assertEqual(self.theorem["total_source_assignment_count"], 256)
        self.assertEqual(self.theorem["phase_required_atom_count"], 3)
        self.assertEqual(self.theorem["selector_required_atom_count"], 2)
        self.assertEqual(self.theorem["shared_phase_selector_atoms"], ["chi11_selector_source_theorem"])
        self.assertEqual(self.theorem["shared_phase_selector_atom_count"], 1)
        self.assertEqual(self.theorem["phase_selector_union_atom_count"], 4)

    def test_minimal_coreadiness_and_quadrants(self) -> None:
        self.assertEqual(self.theorem["phase_minimal_readiness_mask_count"], 1)
        self.assertEqual(self.theorem["selector_minimal_readiness_mask_count"], 1)
        self.assertEqual(self.theorem["phase_selector_coreadiness_minimal_mask_count"], 1)
        self.assertEqual(self.theorem["phase_selector_coreadiness_minimal_size"], 4)
        self.assertEqual(self.theorem["phase_selector_coreadiness_minimal_masks"], [204])
        self.assertEqual(
            self.theorem["quadrant_counts"],
            {
                "neither_phase_nor_selector_ready": 176,
                "phase_and_selector_ready": 16,
                "phase_only_ready": 16,
                "selector_only_ready": 48,
            },
        )

    def test_chi11_necessary_but_not_sufficient(self) -> None:
        self.assertTrue(self.theorem["no_phase_readiness_without_chi11"])
        self.assertTrue(self.theorem["no_selector_readiness_without_chi11"])
        self.assertTrue(self.theorem["chi11_is_common_necessary_cut"])
        self.assertTrue(self.theorem["qw2191_not_sufficient_for_selector"])
        self.assertTrue(self.theorem["chi11_not_sufficient_for_selector"])
        self.assertTrue(self.theorem["chi11_not_sufficient_for_phase"])
        deletion = {row["removed_atom"]: row["blocked_components_after_deletion"] for row in self.finite["phase_selector_coreadiness_deletion_sensitivity_rows"]}
        self.assertEqual(
            deletion["chi11_selector_source_theorem"],
            ["phase_frequency_topological_bit_passage", "selector_source_premise"],
        )
        self.assertEqual(deletion["qw2191_symmetry_breaking_or_internal_source_theorem"], ["selector_source_premise"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertTrue(self.theorem["p2418_chi11_top_incidence_inherited"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("chi11 phase-selector", MD.read_text(encoding="utf-8"))
        self.assertIn("P2419/S1369 chi11", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2419/S1369 chi11", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
