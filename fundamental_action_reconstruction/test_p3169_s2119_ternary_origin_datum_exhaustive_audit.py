import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3169_s2119_ternary_origin_datum_exhaustive_audit.py"
OUT = ROOT / "generated" / "p3169_s2119_ternary_origin_datum_exhaustive_audit.json"
MD = ROOT / "generated" / "p3169_s2119_ternary_origin_datum_exhaustive_audit.md"


class P3169TernaryOriginDatumExhaustiveAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_exhaustive_counts(self):
        self.assertEqual(self.payload["status"], "P3169_TERNARY_ORIGIN_DATUM_EXHAUSTIVE_AUDIT_BOUNDED_NO_GO")
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["ternary_profiles_nonzero"], 3**12 - 1)
        self.assertEqual(cert["translation_classes"], 44367)
        self.assertEqual(cert["affine_translation_classes"], 12768)
        self.assertEqual(cert["trivial_translation_stabilizer_profiles"], 530640)
        self.assertEqual(cert["nonzero_first_resultant_profiles"], 529416)
        self.assertEqual(cert["unique_max_profiles"], 24588)
        self.assertEqual(cert["unique_min_profiles"], 24588)
        self.assertEqual(cert["accepted_Lambda_origin_sources"], 0)

    def test_receivers_are_not_strict_sources(self):
        gates = self.payload["constructed_theoretical_objects"]["gate_rows"]
        gate_map = {row["gate"]: row["passed"] for row in gates}
        self.assertTrue(gate_map["new_nonbinary_candidate_class_constructed"])
        self.assertTrue(gate_map["translation_breaking_receivers_exist"])
        self.assertTrue(gate_map["phase_resultant_receivers_exist"])
        self.assertFalse(gate_map["strict_nadsoliton_provenance_law_exported"])
        self.assertFalse(gate_map["absolute_origin_representative_without_convention"])
        self.assertFalse(gate_map["Phi_Info_A_phi_coupling_theorem_exported"])
        self.assertFalse(gate_map["accepted_Lambda_origin_source"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_docs_updated_and_recommendation(self):
        self.assertIn("P3169/S2119", MD.read_text(encoding="utf-8"))
        self.assertIn("P3169/S2119", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3169/S2119", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3169", (REPO / "AGENTS.md").read_text(encoding="utf-8"))
        self.assertIn("scale-charged strict S_+", self.payload["decision"]["next_honest_step"])


if __name__ == "__main__":
    unittest.main()
