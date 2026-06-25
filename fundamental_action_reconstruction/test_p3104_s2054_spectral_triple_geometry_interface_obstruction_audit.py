import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3104_s2054_spectral_triple_geometry_interface_obstruction_audit.py"
OUT = ROOT / "generated" / "p3104_s2054_spectral_triple_geometry_interface_obstruction_audit.json"
MD = ROOT / "generated" / "p3104_s2054_spectral_triple_geometry_interface_obstruction_audit.md"

class P3104SpectralTripleGeometryInterfaceObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3104_SPECTRAL_TRIPLE_GEOMETRY_INTERFACE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3103"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3103_accepted_nonimported_hilbert_sources"], 0)
        self.assertEqual(cert["algebra_representation_rows"], 12)
        self.assertEqual(cert["dirac_spectrum_rows"], 12)
        self.assertEqual(cert["dirac_rows_with_physical_units"], 0)
        self.assertEqual(cert["distance_formula_rows"], 66)
        self.assertEqual(cert["distance_rows_with_physical_length_units"], 0)
        self.assertEqual(cert["distance_rows_with_alpha_geo_unit_conversion"], 0)
        self.assertEqual(cert["commutator_bound_rows"], 4)
        self.assertEqual(cert["geometry_candidates"], 6)
        self.assertEqual(cert["required_gates"], 9)
        self.assertEqual(cert["candidate_gate_rows"], 54)
        self.assertEqual(cert["accepted_nonimported_geometry_sources"], 0)

    def test_witness_rows_are_formal_and_unsourced(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["represented_on_C12"] for row in objs["algebra_representation_rows"]))
        self.assertTrue(all(not row["physical_coordinate_chart_attached"] for row in objs["algebra_representation_rows"]))
        self.assertTrue(all(row["finite_compact_resolvent_proxy"] for row in objs["dirac_spectrum_rows"]))
        self.assertTrue(all(not row["physical_dirac_unit_attached"] for row in objs["dirac_spectrum_rows"]))
        self.assertTrue(all(row["distance_formula_proxy"] for row in objs["distance_formula_rows"]))
        self.assertTrue(all(not row["physical_length_unit_attached"] for row in objs["distance_formula_rows"]))
        self.assertTrue(all(not row["alpha_geo_converted_to_unit"] for row in objs["distance_formula_rows"]))
        self.assertTrue(all(row["bounded_commutator_proxy"] for row in objs["commutator_bound_rows"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("alpha_geo", decision["next_honest_step"])
        self.assertIn("entropy-to-action/length", decision["next_honest_step"])
        self.assertIn("P3104/S2054", MD.read_text(encoding="utf-8"))
        self.assertIn("P3104/S2054", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3104/S2054", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3104", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
