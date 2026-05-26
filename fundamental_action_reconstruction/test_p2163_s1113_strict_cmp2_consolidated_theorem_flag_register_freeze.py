from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2163_s1113_strict_cmp2_consolidated_theorem_flag_register_freeze.json"


class TestP2163StrictCmp2ConsolidatedTheoremFlagRegisterFreeze(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2163_s1113_strict_cmp2_consolidated_theorem_flag_register_freeze.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2163_s1113_v1")
        self.assertTrue(payload["gatekeeper_checks"]["register_freeze_exported"])
        self.assertGreater(payload["consolidated_theorem_flag_register_freeze"]["n_packets"], 0)


if __name__ == "__main__":
    unittest.main()
