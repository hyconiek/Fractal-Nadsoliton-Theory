import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3131_s2081_epsilon_ot_origin_torsion_twist_invariant_audit.py"
OUT = ROOT / "generated" / "p3131_s2081_epsilon_ot_origin_torsion_twist_invariant_audit.json"
MD = ROOT / "generated" / "p3131_s2081_epsilon_ot_origin_torsion_twist_invariant_audit.md"


class P3131EpsilonOTOriginTorsionTwistAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3131_EPSILON_OT_ORIGIN_TORSION_TWIST_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3130"])

    def test_finite_section_and_character_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3130_accepted_Theta_TO_sources"], 0)
        self.assertEqual(cert["translation_classes"], 351)
        self.assertEqual(cert["equivariant_origin_section_classes"], 1)
        self.assertEqual(cert["orbit_size_distribution"], {"1": 1, "2": 1, "3": 2, "4": 3, "6": 9, "12": 335})
        self.assertEqual(cert["stabilizer_size_distribution"], {"1": 335, "2": 9, "3": 3, "4": 2, "6": 1, "12": 1})
        self.assertEqual(cert["translation_character_twists"], 12)
        self.assertEqual(cert["nontrivial_character_twists"], 11)
        self.assertEqual(cert["origin_selecting_character_twists"], 0)
        self.assertEqual(cert["candidate_Epsilon_OT_sources"], 15)
        self.assertEqual(cert["epsilon_law_rows"], 195)
        self.assertEqual(cert["candidate_gate_rows"], 255)
        self.assertEqual(cert["accepted_Epsilon_OT_sources"], 0)

    def test_candidate_boundaries(self):
        objs = self.payload["constructed_theoretical_objects"]
        section_rows = objs["equivariant_section_obstruction_rows"]
        self.assertEqual(len(section_rows), 351)
        self.assertEqual(sum(row["equivariant_global_origin_section_exists"] for row in section_rows), 1)
        chars = objs["translation_character_twist_rows"]
        self.assertTrue(any(row["character_k"] == 1 and row["nontrivial"] and not row["origin_selecting"] for row in chars))
        candidates = objs["candidate_Epsilon_OT_sources"]
        self.assertTrue(any(row["candidate"] == "c12_character_twist" and row["nonzero_twist_exported"] and not row["absolute_origin_representative_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "equivariant_section_fixed_orbit" and row["equivariant_origin_section_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_origin_torsion" and not row["selector_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Epsilon_OT_source"] for row in objs["candidate_aggregate_certificate"]))

    def test_recommendation_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["finite_equivariant_section_obstruction_proved"])
        self.assertTrue(decision["positive_scoped_flags"]["C12_character_twists_enumerated"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("only the fixed all-support orbit", decision["bounded_result"])
        self.assertIn("Zeta_OS", decision["next_honest_step"])
        self.assertIn("support-local section law", decision["next_honest_step"])
        self.assertIn("P3131/S2081", MD.read_text(encoding="utf-8"))
        self.assertIn("P3131/S2081", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3131/S2081", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3131", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
