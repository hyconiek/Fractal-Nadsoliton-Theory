import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3057_s2007_source_polarity_carrier_extension_sat_verifier.py"
OUT = ROOT / "generated" / "p3057_s2007_source_polarity_carrier_extension_sat_verifier.json"
MD = ROOT / "generated" / "p3057_s2007_source_polarity_carrier_extension_sat_verifier.md"

class P3057SourcePolarityCarrierExtensionSatVerifierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3057_SOURCE_POLARITY_CARRIER_EXTENSION_SAT_VERIFIER_TEMPLATE_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3056"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["core_rows"], 6)
        self.assertEqual(cert["compatibility_squares"], 6)
        self.assertEqual(cert["extension_atoms"], 12)
        self.assertEqual(cert["minimal_extension_count"], 1)
        self.assertEqual(cert["minimal_extension_size"], 12)
        self.assertEqual(cert["primitive_rows_in_minimal_extension"], 2)
        self.assertEqual(cert["row_import_certificates_in_minimal_extension"], 4)
        self.assertEqual(cert["compatibility_squares_in_minimal_extension"], 6)
        self.assertEqual(cert["single_atom_obstruction_rows"], 12)
        self.assertEqual(cert["single_atom_accepts"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 3)

    def test_template_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        template = obj["theorem_template"]
        self.assertEqual(template["object"], "StrictSourcePolarityCarrierExtensionTheoremTemplate")
        self.assertEqual(template["carrier_symbol"], "G_selector")
        self.assertEqual(template["minimal_atom_count"], 12)
        self.assertEqual(len(template["required_atoms"]), 12)
        for value in obj["current_single_carrier_state"].values():
            self.assertFalse(value)
        for row in obj["obstruction_cut_table"]:
            self.assertFalse(row["accepted_alone"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3057/S2007", MD.read_text(encoding="utf-8"))
        self.assertIn("P3057/S2007", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3057/S2007", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3057", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
