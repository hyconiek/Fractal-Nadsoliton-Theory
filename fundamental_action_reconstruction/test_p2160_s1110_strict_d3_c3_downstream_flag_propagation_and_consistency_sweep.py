from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2160_s1110_strict_d3_c3_downstream_flag_propagation_and_consistency_sweep.json"


class TestP2160StrictD3C3DownstreamFlagPropagationAndConsistencySweep(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2160_s1110_strict_d3_c3_downstream_flag_propagation_and_consistency_sweep.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2160_s1110_v1")
        self.assertTrue(payload["gatekeeper_checks"]["sweep_exported"])
        self.assertTrue(payload["downstream_flag_propagation_and_consistency_sweep"]["consistency_ok"])


if __name__ == "__main__":
    unittest.main()
