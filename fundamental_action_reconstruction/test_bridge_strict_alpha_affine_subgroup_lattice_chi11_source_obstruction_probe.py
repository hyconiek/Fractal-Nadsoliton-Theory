import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_affine_subgroup_lattice_chi11_source_obstruction_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_affine_subgroup_lattice_chi11_source_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_affine_subgroup_lattice_chi11_source_obstruction_report.md"


class StrictAlphaAffineSubgroupLatticeChi11SourceObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_AFFINE_SUBGROUP_LATTICE_CHI11_SOURCE_OBSTRUCTION__NO_FALSE_PASS",
        )
        self.assertEqual(
            payload["status"],
            "only-cyclic-and-D12-quotients-admit-chi11-unit5-or-unit7-kill-it",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["all_automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["affine_unit_subgroup_count"], 5)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_lattice_summary(self):
        summary = self.payload["lattice_summary"]
        self.assertEqual(summary["admitting_subgroup_count"], 2)
        self.assertEqual(summary["killing_subgroup_count"], 3)
        self.assertEqual(
            summary["admitting_subgroups"],
            ["T_semidirect_{1}_cyclic_translation", "T_semidirect_{1,11}_D12"],
        )
        self.assertEqual(
            summary["killing_subgroups"],
            [
                "T_semidirect_{1,5}_unit5_axis",
                "T_semidirect_{1,7}_unit7_axis",
                "T_semidirect_{1,5,7,11}_full_affine_Aut",
            ],
        )
        self.assertEqual(summary["minimal_nontrivial_admitting_reflection_subgroup"], "T_semidirect_{1,11}_D12")
        self.assertEqual(summary["minimal_killing_unit_additions"], ["unit5", "unit7"])
        self.assertFalse(summary["full_aut_admits_chi11"])

    def test_subgroup_rows_and_branch_partitions(self):
        rows = {row["subgroup_name"]: row for row in self.payload["subgroup_rows"]}
        self.assertEqual(set(rows), {
            "T_semidirect_{1}_cyclic_translation",
            "T_semidirect_{1,11}_D12",
            "T_semidirect_{1,5}_unit5_axis",
            "T_semidirect_{1,7}_unit7_axis",
            "T_semidirect_{1,5,7,11}_full_affine_Aut",
        })
        self.assertEqual(rows["T_semidirect_{1}_cyclic_translation"]["orbit_count_on_all_supports"], 66)
        self.assertEqual(rows["T_semidirect_{1,11}_D12"]["orbit_count_on_all_supports"], 38)
        self.assertEqual(rows["T_semidirect_{1,5}_unit5_axis"]["orbit_count_on_all_supports"], 38)
        self.assertEqual(rows["T_semidirect_{1,7}_unit7_axis"]["orbit_count_on_all_supports"], 40)
        self.assertEqual(rows["T_semidirect_{1,5,7,11}_full_affine_Aut"]["orbit_count_on_all_supports"], 25)

        self.assertTrue(rows["T_semidirect_{1}_cyclic_translation"]["branch_chi11_well_defined_on_quotient"])
        self.assertTrue(rows["T_semidirect_{1,11}_D12"]["branch_chi11_well_defined_on_quotient"])
        self.assertFalse(rows["T_semidirect_{1,5}_unit5_axis"]["branch_chi11_well_defined_on_quotient"])
        self.assertFalse(rows["T_semidirect_{1,7}_unit7_axis"]["branch_chi11_well_defined_on_quotient"])
        self.assertFalse(rows["T_semidirect_{1,5,7,11}_full_affine_Aut"]["branch_chi11_well_defined_on_quotient"])

        self.assertEqual(rows["T_semidirect_{1,11}_D12"]["branch_orbit_partition"], {
            "0": ["A1_k1", "A11_k11"],
            "37": ["A5_k5", "A7_k7"],
        })
        self.assertEqual(rows["T_semidirect_{1,5,7,11}_full_affine_Aut"]["branch_orbit_partition"], {
            "0": ["A1_k1", "A5_k5", "A7_k7", "A11_k11"],
        })
        self.assertEqual(len(rows["T_semidirect_{1,5,7,11}_full_affine_Aut"]["mixed_chi11_sign_branch_orbits"]), 1)

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("no quotient orbit contains both signs", proof["well_definedness_test"])
        self.assertIn("cyclic_translation", proof["admitting_subgroups"])
        self.assertIn("full_affine_Aut", proof["killing_subgroups"])
        self.assertIn("unit 5 or unit 7", proof["unit_obstruction"])
        self.assertIn("cannot export chi_11 polarity", proof["full_aut_consequence"])
        self.assertIn("does not discharge QW-2191", proof["strict_limit"])

        interpretation = self.payload["interpretation"]
        self.assertIn("D12 is the largest", interpretation["honest_positive"])
        self.assertIn("unit 5 or unit 7", interpretation["honest_negative"])
        self.assertIn("which subgroup symmetries erase", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives chi_11", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertIn("does not supply an internal strict source", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
