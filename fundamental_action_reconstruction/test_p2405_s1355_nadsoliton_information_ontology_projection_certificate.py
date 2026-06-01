#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2405_s1355_nadsoliton_information_ontology_projection_certificate.py"
OUT = ROOT / "generated" / "p2405_s1355_nadsoliton_information_ontology_projection_certificate.json"
MD = ROOT / "generated" / "p2405_s1355_nadsoliton_information_ontology_projection_certificate.md"
P2404 = ROOT / "generated" / "p2404_s1354_strict_addition_physics_lane_dependency_cut_certificate.json"


class P2405NadsolitonInformationOntologyProjectionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2404.exists():
            subprocess.run([sys.executable, str(ROOT / "p2404_s1354_strict_addition_physics_lane_dependency_cut_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["nadsoliton_information_ontology_projection_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.poset = cls.cert["ontology_poset"]
        cls.lattice = cls.cert["ontology_guard_lattice"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2405")
        self.assertEqual(self.payload["stage_id"], "S1355")
        self.assertEqual(self.payload["status"], "PASS_INFORMATION_ONTOLOGY_RETYPE_NO_UNDERLAYER_NO_ROLE_EXPORT")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_ax9_pure_information_inherited(self) -> None:
        self.assertEqual(self.theorem["single_primitive"], "informational_nadsoliton_only")
        self.assertFalse(self.theorem["separate_information_layer_allowed"])
        self.assertIn("the nadsoliton is itself the primordial information", self.theorem["canonical_ontology_statement"])

    def test_ontology_guard_lattice(self) -> None:
        self.assertEqual(self.lattice["row_count"], 32)
        self.assertEqual(self.theorem["ontology_guard_true_masks"], [31])
        self.assertEqual(self.theorem["ontology_guard_proper_subset_fail_count"], 31)
        self.assertEqual(self.theorem["ontology_guard_anf_degree"], 5)
        self.assertIn("no_separate_information_layer_under_nadsoliton", self.theorem["ontology_guard_anf"])

    def test_poset_keeps_nadsoliton_as_unique_root(self) -> None:
        self.assertEqual(self.poset["root_nodes"], ["nadsoliton_pure_information_root"])
        self.assertFalse(self.poset["has_cycle"])
        self.assertTrue(self.theorem["unique_information_root"])
        self.assertTrue(self.theorem["all_nodes_reachable_from_information_root"])
        self.assertEqual(self.poset["preferred_internal_order"], ["nadsoliton", "light", "matter", "emergent_observer"])

    def test_interpretation_matrix_blocks_underlayer_and_role_export(self) -> None:
        invalid = set(self.theorem["invalid_interpretations"])
        valid = set(self.theorem["valid_interpretations"])
        self.assertIn("separate_information_substrate_below_nadsoliton", invalid)
        self.assertIn("strict_additions_as_immediate_physical_role_exports", invalid)
        self.assertIn("strict_additions_as_internal_nadsoliton_information_constraints", valid)
        self.assertIn("physical_roles_as_downstream_projection_theorems", valid)

    def test_p2404_dependency_inherited_without_role_export(self) -> None:
        self.assertEqual(self.theorem["p2404_strict_additions_only_physical_role_lanes_ready"], [])
        self.assertEqual(self.theorem["p2404_ltotal_dependency_degree"], 7)
        self.assertEqual(
            self.theorem["strict_information_atoms"],
            [
                "apd_completion_structure",
                "gf2_phase_topological_data",
                "nonlinear_damping_compression",
                "chi11_selector_source_declared",
            ],
        )

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("information-ontology projection certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2405/S1355 nadsoliton", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2405/S1355 information-ontology", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
