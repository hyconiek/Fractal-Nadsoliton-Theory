#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.py"
OUT = ROOT / "generated" / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.json"
MD = ROOT / "generated" / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.md"
P2405 = ROOT / "generated" / "p2405_s1355_nadsoliton_information_ontology_projection_certificate.json"


class P2406InformationToPhysicsStagedProjectionBarrierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2405.exists():
            subprocess.run([sys.executable, str(ROOT / "p2405_s1355_nadsoliton_information_ontology_projection_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["information_to_physics_staged_projection_barrier_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.summary = cls.cert["staged_projection_summary"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2406")
        self.assertEqual(self.payload["stage_id"], "S1356")
        self.assertEqual(self.payload["status"], "PASS_STAGED_PROJECTION_BARRIER_NO_DIRECT_ROLE_EXPORT_NO_TOE_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_compact_4096_assignment_space(self) -> None:
        self.assertEqual(self.theorem["total_boolean_assignment_count"], 4096)
        self.assertTrue(self.theorem["summary_only_no_full_truth_table_emitted"])
        self.assertEqual(len(self.summary["atom_order"]), 12)

    def test_staged_readiness_barrier(self) -> None:
        self.assertEqual(self.theorem["ontology_only_ready_lanes"], ["ontology_typed_information_root"])
        self.assertEqual(
            self.theorem["ontology_plus_strict_ready_lanes"],
            ["ontology_typed_information_root", "strict_internal_information_completion_ready"],
        )
        self.assertEqual(self.theorem["proper_role_prefix_fail_count_after_ontology_plus_strict"], 7)
        self.assertIn("toe_downstream_projection_candidate", self.theorem["full_projection_ready_lanes"])

    def test_ltotal_and_toe_degree_twelve_singleton(self) -> None:
        self.assertEqual(self.theorem["ltotal_projection_anf_degree"], 12)
        self.assertEqual(self.theorem["toe_projection_anf_degree"], 12)
        self.assertEqual(self.theorem["ltotal_projection_true_assignment_count"], 1)
        self.assertEqual(self.theorem["toe_projection_true_assignment_count"], 1)
        self.assertIn("no_separate_information_layer_under_nadsoliton", self.theorem["ltotal_projection_anf"])
        self.assertIn("beta_power_hierarchy_successor_theorem", self.theorem["toe_projection_anf"])

    def test_inherited_p2404_p2405_guards(self) -> None:
        self.assertTrue(self.theorem["p2404_degree_seven_inherited"])
        self.assertTrue(self.theorem["p2405_degree_five_ontology_guard_inherited"])
        self.assertTrue(self.theorem["p2405_no_underlayer_inherited"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("staged projection barrier certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2406/S1356 information-to-physics", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2406/S1356 staged projection", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
