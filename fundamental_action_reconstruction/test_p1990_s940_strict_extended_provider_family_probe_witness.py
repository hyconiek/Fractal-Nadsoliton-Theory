from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1990_s940_strict_extended_provider_family_probe_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1990_s940_strict_extended_provider_family_probe_witness.json"

class TestP1990S940(unittest.TestCase):
    def test_extended_provider_probe(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P1990")
        self.assertEqual(p["result_kind"], "PASS_EXTENDED_PROVIDER_PROBE_WITH_NONSTRICT_LABEL")
        self.assertEqual(p["selector_premise_label"]["status"], "NON_STRICT_AUGMENTED_CLASS")
        self.assertTrue(p["gatekeeper_checks"]["augmented_channel_solution_exists"])
        self.assertTrue(p["gatekeeper_checks"]["mixed_channels_still_persist"])
        self.assertTrue(p["gatekeeper_checks"]["numeric_residual_after_augmented_fit_nonzero"])

if __name__ == "__main__":
    unittest.main()
