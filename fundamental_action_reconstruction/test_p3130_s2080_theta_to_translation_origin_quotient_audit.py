import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3130_s2080_theta_to_translation_origin_quotient_audit.py"
OUT = ROOT / "generated" / "p3130_s2080_theta_to_translation_origin_quotient_audit.json"
MD = ROOT / "generated" / "p3130_s2080_theta_to_translation_origin_quotient_audit.md"


class P3130ThetaTOTranslationOriginQuotientAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3130_THETA_TO_TRANSLATION_ORIGIN_QUOTIENT_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3129"])

    def test_exhaustive_quotient_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3129_accepted_Gamma_SO_generators"], 0)
        self.assertEqual(cert["nonempty_Z12_supports"], 4095)
        self.assertEqual(cert["translation_classes"], 351)
        self.assertEqual(cert["dihedral_classes"], 223)
        self.assertEqual(cert["translation_classes_with_absolute_origin"], 1)
        self.assertGreater(cert["inversion_fixed_translation_classes"], 0)
        self.assertEqual(cert["candidate_Theta_TO_sources"], 15)
        self.assertEqual(cert["quotient_law_rows"], 210)
        self.assertEqual(cert["symmetry_witness_rows"], 105)
        self.assertEqual(cert["candidate_gate_rows"], 270)
        self.assertEqual(cert["accepted_Theta_TO_sources"], 0)
        self.assertGreater(cert["A_phi"], 0)

    def test_candidate_boundaries(self):
        objs = self.payload["constructed_theoretical_objects"]
        rows = objs["full_Z12_translation_quotient_summary"]["translation_class_rows"]
        self.assertEqual(len(rows), 351)
        self.assertTrue(any(row["absolute_origin_fixed_by_quotient"] for row in rows))
        self.assertTrue(any(not row["absolute_origin_fixed_by_quotient"] for row in rows))
        candidates = objs["candidate_Theta_TO_sources"]
        self.assertTrue(any(row["candidate"] == "necklace_translation_quotient" and row["translation_quotient_exported"] and not row["absolute_origin_representative_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "chiral_signed_quotient" and row["nonzero_internal_sign_exported"] and not row["absolute_origin_representative_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_zero_origin" and not row["selector_free"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "lagrangian_coordinate_origin" and not row["not_imported_dynamics"] for row in candidates))
        self.assertTrue(all(not row["accepted_Theta_TO_source"] for row in objs["candidate_aggregate_certificate"]))

    def test_recommendation_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["full_Z12_translation_quotient_exhausted"])
        self.assertTrue(decision["positive_scoped_flags"]["strict_translation_quotient_object_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("4095 supports collapse to 351 translation classes", decision["bounded_result"])
        self.assertIn("Epsilon_OT", decision["next_honest_step"])
        self.assertIn("origin-torsion", decision["next_honest_step"])
        self.assertIn("P3130/S2080", MD.read_text(encoding="utf-8"))
        self.assertIn("P3130/S2080", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3130/S2080", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3130", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
