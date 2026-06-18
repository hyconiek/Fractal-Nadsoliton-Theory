import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2857_SCRIPT = ROOT / "p2857_s1807_observer_readout_phase_source_effect_audit.py"
SCRIPT = ROOT / "p2858_s1808_phase_bit_cell_continuum_no_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2858_s1808_phase_bit_cell_continuum_no_source_audit.json"
MD_PATH = ROOT / "generated" / "p2858_s1808_phase_bit_cell_continuum_no_source_audit.md"


class P2858PhaseBitCellContinuumNoSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2857_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["phase_bit_cell_continuum_no_source_audit"]

    def test_status_and_p2857_input(self):
        self.assertEqual(self.payload["status"], "P2858_PHASE_BIT_CELL_CONTINUUM_NO_SOURCE_AUDIT_NO_CLOSURE")
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2857_OBSERVER_READOUT_PHASE_SOURCE_EFFECT_AUDIT_NO_CLOSURE",
        )

    def test_every_phase_cell_row_has_positive_margin(self):
        rows = self.audit["cell_rows"]
        self.assertEqual(len(rows), 12)
        for row in rows:
            self.assertTrue(row["sign_constraint_satisfied"])
            self.assertGreater(row["distance_to_nearest_cos_zero"], 0.0)
            self.assertGreater(row["safe_common_epsilon_bound"], 0.0)

    def test_open_box_and_rational_probes_preserve_bits(self):
        box = self.audit["open_box_certificate"]
        target = self.audit["target_phase_bits"]
        self.assertGreater(box["certified_open_box_half_width"], 0.0)
        self.assertTrue(box["rational_probe_delta_inside_box"])
        self.assertGreater(len(box["probe_phase_bits"]), 0)
        for bits in box["probe_phase_bits"]:
            self.assertEqual(bits, target)
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["rational_probes_preserve_phase_bits"])

    def test_no_source_or_closure_export(self):
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_phase_bit_open_cell_no_source_audit"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_strict_phase_frequency_source_law"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["phase_bit_source_law_exported"])
        self.assertFalse(flags["strict_phase_frequency_source_law_exported"])
        self.assertFalse(flags["full_kernel_bridge_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2858/S1808", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2858/S1808", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2858/S1808", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2858", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
