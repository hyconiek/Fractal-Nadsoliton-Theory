#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2404_s1354_strict_addition_physics_lane_dependency_cut_certificate.py"
OUT = ROOT / "generated" / "p2404_s1354_strict_addition_physics_lane_dependency_cut_certificate.json"
MD = ROOT / "generated" / "p2404_s1354_strict_addition_physics_lane_dependency_cut_certificate.md"
PREREQ = ROOT / "generated" / "p2403_s1353_strict_kernel_primary_physics_generation_rebase_certificate.json"


class P2404StrictAdditionPhysicsLaneDependencyCutCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not PREREQ.exists():
            subprocess.run([sys.executable, str(ROOT / "p2403_s1353_strict_kernel_primary_physics_generation_rebase_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_addition_physics_lane_dependency_cut_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.lattice = cls.cert["dependency_lattice"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2404")
        self.assertEqual(self.payload["stage_id"], "S1354")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.payload["status"], "PASS_STRICT_ADDITION_DEPENDENCY_CUT_NO_ROLE_TRANSFER_NO_TOE_CLOSURE")

    def test_exact_128_row_dependency_lattice(self) -> None:
        self.assertEqual(self.theorem["truth_table_row_count"], 128)
        self.assertEqual(len(self.lattice["rows"]), 128)
        self.assertEqual(self.lattice["full_mask"], 127)
        self.assertEqual(self.theorem["full_mask_ready_lane_count"], 7)

    def test_common_strict_addition_cut(self) -> None:
        expected = [
            "apd_completion_structure",
            "gf2_phase_topological_data",
            "nonlinear_damping_compression",
            "chi11_selector_source_declared",
        ]
        self.assertEqual(self.theorem["common_strict_addition_cut_support"], expected)
        self.assertEqual(self.theorem["common_strict_addition_cut_mask"], 15)
        self.assertTrue(self.theorem["p2403_strict_additions_matched"])

    def test_legacy_and_strict_additions_only_slices(self) -> None:
        self.assertEqual(self.theorem["legacy_only_ready_lanes"], [])
        self.assertEqual(
            self.theorem["strict_additions_only_ready_lanes"],
            ["strict_kernel_structural_candidate_test_readiness", "strict_mass_generation_candidate_test_readiness"],
        )
        self.assertEqual(self.theorem["strict_additions_only_physical_role_lanes_ready"], [])

    def test_ltotal_and_toe_degree_seven_dependencies(self) -> None:
        monomial = "apd_completion_structure * gf2_phase_topological_data * nonlinear_damping_compression * chi11_selector_source_declared * alpha_geo_electroweak_role_theorem * beta_tors_strict_role_theorem * beta_power_hierarchy_successor_theorem"
        self.assertEqual(self.theorem["ltotal_dependency_anf"], monomial)
        self.assertEqual(self.theorem["toe_package_dependency_anf"], monomial)
        self.assertEqual(self.theorem["ltotal_dependency_degree"], 7)
        self.assertEqual(self.theorem["toe_package_dependency_degree"], 7)

    def test_distance_from_legacy_only(self) -> None:
        distance = self.theorem["distance_from_legacy_only_by_lane"]
        self.assertEqual(distance["strict_kernel_structural_candidate_test_readiness"], 4)
        self.assertEqual(distance["strict_mass_generation_candidate_test_readiness"], 4)
        self.assertEqual(distance["legacy_weinberg_role_transfer_to_strict_successor"], 5)
        self.assertEqual(distance["legacy_alpha_em_role_transfer_to_strict_successor"], 6)
        self.assertEqual(distance["legacy_gravity_hierarchy_strict_successor"], 6)
        self.assertEqual(distance["strict_role_bearing_ltotal_promotion_candidate"], 7)
        self.assertEqual(distance["strict_toe_physics_generation_package_candidate"], 7)

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertTrue(self.theorem["p2400_full_role_mask_requirement_inherited"])
        self.assertIn("dependency-cut certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2404/S1354 strict-addition", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2404/S1354 dependency-cut", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
