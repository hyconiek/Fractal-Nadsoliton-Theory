from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p2002_s952_strict_cutkosky_missing_channel_vs_structural_classifier_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p2002_s952_strict_cutkosky_missing_channel_vs_structural_classifier_witness.json"


class TestP2002S952(unittest.TestCase):
    def test_classifier(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P2002")
        self.assertEqual(p["result_kind"], "PASS_DELTA_NORM_CLASSIFIER_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["shared_grid_length"])
        self.assertIn(
            p["classifier"],
            {"MISSING_CHANNEL_PRESSURE_SUPPORTED", "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED", "MIXED_OR_INCONCLUSIVE"},
        )


if __name__ == "__main__":
    unittest.main()
