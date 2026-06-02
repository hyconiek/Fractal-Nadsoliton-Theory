#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2418_s1368_bridge_source_marginal_unlock_lattice_certificate.py"
OUT = ROOT / "generated" / "p2418_s1368_bridge_source_marginal_unlock_lattice_certificate.json"
MD = ROOT / "generated" / "p2418_s1368_bridge_source_marginal_unlock_lattice_certificate.md"
P2417 = ROOT / "generated" / "p2417_s1367_apd_witness_to_source_obligation_nonpromotion_matrix_certificate.json"


class P2418BridgeSourceMarginalUnlockLatticeCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2417.exists():
            subprocess.run([sys.executable, str(ROOT / "p2417_s1367_apd_witness_to_source_obligation_nonpromotion_matrix_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["bridge_source_marginal_unlock_lattice_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_witness_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2418")
        self.assertEqual(self.payload["stage_id"], "S1368")
        self.assertEqual(self.payload["status"], "PASS_SOURCE_MARGINAL_UNLOCK_LATTICE_NO_SINGLETON_NO_BRIDGE_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_lattice_counts(self) -> None:
        self.assertEqual(self.theorem["source_obligation_atom_count"], 8)
        self.assertEqual(self.theorem["bridge_component_count"], 5)
        self.assertEqual(self.theorem["total_source_assignment_count"], 256)
        self.assertEqual(self.theorem["current_source_mask_from_p2417"], 0)
        self.assertEqual(self.theorem["full_source_mask"], 255)
        self.assertEqual(self.theorem["bridge_ready_masks"], [255])
        self.assertEqual(self.theorem["current_mask_ready_components"], ["residual_strict_additions_inventory"])

    def test_marginal_unlock_structure(self) -> None:
        self.assertEqual(self.theorem["singleton_unlock_count"], 0)
        self.assertEqual(self.theorem["pair_unlock_count"], 3)
        self.assertEqual(self.theorem["minimum_nonresidual_component_unlock_size"], 2)
        self.assertEqual(self.theorem["phase_component_minimal_size"], 3)
        self.assertEqual(self.theorem["selector_component_minimal_size"], 2)
        pair_components = [tuple(row["unlocked_nonresidual_components"]) for row in self.finite["pair_unlock_rows"]]
        self.assertIn(("amplitude_normalization_passage",), pair_components)
        self.assertIn(("damping_compression_passage",), pair_components)
        self.assertIn(("selector_source_premise",), pair_components)

    def test_priority_and_limits(self) -> None:
        self.assertTrue(self.theorem["chi11_top_incidence_inherited"])
        self.assertTrue(self.theorem["p2417_zero_source_mask_inherited"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertIn("No singleton source atom", "\n".join(self.theorem["not_licensed"]))

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("marginal-unlock", MD.read_text(encoding="utf-8"))
        self.assertIn("P2418/S1368 bridge", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2418/S1368 bridge", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
