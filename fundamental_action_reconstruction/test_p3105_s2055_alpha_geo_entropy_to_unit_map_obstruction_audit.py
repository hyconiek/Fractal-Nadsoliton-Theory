import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit.py"
OUT = ROOT / "generated" / "p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit.json"
MD = ROOT / "generated" / "p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit.md"

class P3105AlphaGeoEntropyToUnitMapObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3105_ALPHA_GEO_ENTROPY_TO_UNIT_MAP_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3104"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3104_accepted_nonimported_geometry_sources"], 0)
        self.assertEqual(cert["entropy_identity_rows"], 4)
        self.assertEqual(cert["alpha_geo_reference_rows"], 1)
        self.assertEqual(cert["entropy_rows_with_unit_dimension"], 0)
        self.assertEqual(cert["self_coupled_action_rows"], 4)
        self.assertEqual(cert["action_rows_with_dimensionful_unit"], 0)
        self.assertEqual(cert["scale_orbit_rows"], 4)
        self.assertEqual(cert["canonical_scale_selections"], 0)
        self.assertEqual(cert["target_unit_rows"], 6)
        self.assertEqual(cert["target_rows_with_dimension_assignment"], 0)
        self.assertEqual(cert["unit_candidates"], 6)
        self.assertEqual(cert["required_gates"], 8)
        self.assertEqual(cert["candidate_gate_rows"], 48)
        self.assertEqual(cert["accepted_nonimported_unit_sources"], 0)

    def test_witness_rows_are_internal_but_unitless(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertEqual(objs["alpha_geo_unit_map_audit_object"]["ontology"], "nadsoliton is self-coupled primordial information; standard-physics units are not assumed")
        self.assertTrue(any(row["is_alpha_geo_reference"] for row in objs["entropy_identity_rows"]))
        self.assertTrue(all(not row["unit_dimension_attached"] for row in objs["entropy_identity_rows"]))
        self.assertTrue(all(row["self_coupled_information_proxy"] for row in objs["self_coupled_information_action_rows"]))
        self.assertTrue(all(not row["dimensionful_action_unit_attached"] for row in objs["self_coupled_information_action_rows"]))
        self.assertTrue(all(not row["canonical_unit_selected"] for row in objs["scale_orbit_rows"]))
        self.assertTrue(all(not row["dimension_assignment_exported"] for row in objs["target_unit_rows"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("self-coupled informational action-density normalization theorem", decision["next_honest_step"])
        self.assertIn("P3105/S2055", MD.read_text(encoding="utf-8"))
        self.assertIn("P3105/S2055", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3105/S2055", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3105", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
