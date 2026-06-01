#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2407_s1357_stage_quotient_projection_barrier_certificate.py"
OUT = ROOT / "generated" / "p2407_s1357_stage_quotient_projection_barrier_certificate.json"
MD = ROOT / "generated" / "p2407_s1357_stage_quotient_projection_barrier_certificate.md"
P2406 = ROOT / "generated" / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.json"


class P2407StageQuotientProjectionBarrierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2406.exists():
            subprocess.run([sys.executable, str(ROOT / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["stage_quotient_projection_barrier_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.lattice = cls.cert["quotient_lattice"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2407")
        self.assertEqual(self.payload["stage_id"], "S1357")
        self.assertEqual(self.payload["status"], "PASS_STAGE_QUOTIENT_BARRIER_FULL_MASK_ONLY_NO_TOE_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_eight_row_quotient_lattice(self) -> None:
        self.assertEqual(self.theorem["row_count"], 8)
        self.assertEqual(len(self.lattice["rows"]), 8)
        self.assertEqual(self.theorem["stage_atoms"], [
            "O_ontology_guard_package",
            "S_strict_internal_completion_package",
            "R_role_successor_projection_package",
        ])

    def test_only_full_stage_mask_readies_ltotal_and_toe(self) -> None:
        self.assertEqual(self.theorem["ltotal_true_masks"], [7])
        self.assertEqual(self.theorem["toe_true_masks"], [7])
        self.assertEqual(self.theorem["proper_subset_fail_count"], 7)
        self.assertEqual(self.theorem["stage_projection_anf"], "O_ontology_guard_package * S_strict_internal_completion_package * R_role_successor_projection_package")
        self.assertEqual(self.theorem["stage_projection_anf_degree"], 3)

    def test_ontology_and_strict_prefixes_remain_nonphysical(self) -> None:
        self.assertEqual(self.theorem["ontology_only_ready_lanes"], ["typed_information_root_ready"])
        self.assertEqual(
            self.theorem["ontology_plus_strict_ready_lanes"],
            ["typed_information_root_ready", "strict_internal_completion_ready"],
        )
        self.assertNotIn("role_bearing_ltotal_projection_ready", self.theorem["ontology_plus_strict_ready_lanes"])

    def test_p2406_degree_twelve_expansion_inherited(self) -> None:
        self.assertTrue(self.theorem["p2406_degree_twelve_inherited"])
        self.assertTrue(self.theorem["quotient_expansion_matches_p2406_degree"])
        self.assertEqual(self.theorem["expanded_stage_atom_count"], 12)

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("stage-quotient projection barrier certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2407/S1357 stage-quotient", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2407/S1357 stage-quotient", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
