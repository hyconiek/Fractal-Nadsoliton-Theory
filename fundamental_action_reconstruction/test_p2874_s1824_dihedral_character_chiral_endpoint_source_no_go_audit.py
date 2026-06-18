import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2873_SCRIPT = ROOT / "p2873_s1823_exhaustive_dihedral_z12_endpoint_predicate_no_go_audit.py"
SCRIPT = ROOT / "p2874_s1824_dihedral_character_chiral_endpoint_source_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2874_s1824_dihedral_character_chiral_endpoint_source_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2874_s1824_dihedral_character_chiral_endpoint_source_no_go_audit.md"


class P2874DihedralCharacterChiralEndpointSourceNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2873_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["dihedral_character_chiral_endpoint_source_no_go_audit"]

    def test_status_and_p2873_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2874_DIHEDRAL_CHARACTER_CHIRAL_ENDPOINT_SOURCE_NO_GO_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2873_EXHAUSTIVE_DIHEDRAL_Z12_ENDPOINT_PREDICATE_NO_GO_AUDIT_NO_CLOSURE",
        )

    def test_all_four_d12_characters_are_enumerated(self):
        self.assertEqual(self.audit["character_count"], 4)
        chars = {(record["character"]["chi_r"], record["character"]["chi_s"]) for record in self.audit["character_records"]}
        self.assertEqual(chars, {(1, 1), (1, -1), (-1, 1), (-1, -1)})
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["all_four_real_d12_characters_enumerated"])

    def test_reflection_odd_characters_vanish(self):
        odd_records = [record for record in self.audit["character_records"] if record["reflection_odd"]]
        self.assertEqual(len(odd_records), 2)
        self.assertTrue(all(record["endpoint_field_dimension"] == 0 for record in odd_records))
        self.assertTrue(all(record["support"] == [] for record in odd_records))
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["reflection_odd_characters_have_zero_endpoint_field"])

    def test_reflection_even_fields_are_global_not_singleton(self):
        even_records = [record for record in self.audit["character_records"] if not record["reflection_odd"]]
        self.assertEqual(len(even_records), 2)
        self.assertTrue(all(record["endpoint_field_dimension"] == 1 for record in even_records))
        self.assertTrue(all(record["support"] == list(range(12)) for record in even_records))
        self.assertTrue(all(not record["singleton_11_support"] for record in even_records))
        kinds = {record["field_kind"] for record in even_records}
        self.assertEqual(kinds, {"constant", "alternating_rotation_parity"})

    def test_alternating_field_is_global_parity_not_endpoint_localizer(self):
        alternating = [record for record in self.audit["character_records"] if record["field_kind"] == "alternating_rotation_parity"]
        self.assertEqual(len(alternating), 1)
        self.assertFalse(alternating[0]["separates_11_from_1_by_value"])
        self.assertEqual(alternating[0]["basis_fields"][0][11], -1)
        self.assertEqual(alternating[0]["basis_fields"][0][1], -1)
        self.assertEqual(alternating[0]["support"], list(range(12)))
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["no_character_field_has_singleton_11_support"])
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["alternating_character_is_global_parity_not_endpoint_11_localization"])

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["orientation_source_law_exported"])
        self.assertFalse(flags["chiral_endpoint_11_source_law_exported"])
        self.assertFalse(flags["selector_or_localizer_source_exported"])
        self.assertFalse(flags["unit_bearing_coupling_localization_theorem_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2874/S1824", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2874/S1824", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2874/S1824", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2874", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
