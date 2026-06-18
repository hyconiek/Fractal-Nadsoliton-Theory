import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2856_SCRIPT = ROOT / "p2856_s1806_prime5_phase_unit_extension_ambiguity_audit.py"
SCRIPT = ROOT / "p2857_s1807_observer_readout_phase_source_effect_audit.py"
JSON_PATH = ROOT / "generated" / "p2857_s1807_observer_readout_phase_source_effect_audit.json"
MD_PATH = ROOT / "generated" / "p2857_s1807_observer_readout_phase_source_effect_audit.md"


class P2857ObserverReadoutPhaseSourceEffectAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2856_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["observer_readout_phase_source_effect_audit"]

    def test_status_and_ontology_order(self):
        self.assertEqual(self.payload["status"], "P2857_OBSERVER_READOUT_PHASE_SOURCE_EFFECT_AUDIT_NO_CLOSURE")
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2856_PRIME5_PHASE_UNIT_EXTENSION_AMBIGUITY_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(self.audit["ontology_order"], "nadsoliton -> light -> matter -> emergent observer")

    def test_observer_frame_orbit_is_finite_and_nonselecting(self):
        orbit = self.audit["observer_frame_orbit"]
        self.assertEqual(orbit["frame_count"], 48)
        self.assertGreaterEqual(orbit["distinct_profile_count"], 1)
        self.assertGreaterEqual(orbit["stabilizer_count"], 1)
        self.assertIn({"unit": 1, "shift": 0}, orbit["stabilizers"])

    def test_bit_observer_cannot_distinguish_p2856_witnesses(self):
        rows = self.audit["witness_observer_indistinguishability"]
        self.assertGreater(len(rows), 0)
        for row in rows:
            self.assertTrue(row["same_base_bits"])
            self.assertTrue(row["same_observer_orbit"])
            self.assertFalse(row["distinguishable_by_bit_observer"])

    def test_observer_modes_export_no_source(self):
        matrix = {row["candidate"]: row for row in self.audit["candidate_matrix"]}
        self.assertFalse(matrix["observer_bit_profile_readout"]["exports_pre_observer_source_law"])
        self.assertFalse(matrix["observer_origin_or_frame_anchor"]["exports_pre_observer_source_law"])
        self.assertFalse(matrix["observer_full_parameter_meter"]["exports_pre_observer_source_law"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_observer_effect_no_source_audit"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_observer_strict_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_strict_phase_frequency_source_law"])

    def test_no_false_closure_documents(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["observer_source_law_exported"])
        self.assertFalse(flags["strict_phase_frequency_source_law_exported"])
        self.assertFalse(flags["full_kernel_bridge_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2857/S1807", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2857/S1807", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2857/S1807", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2857", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
