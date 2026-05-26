from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2165_s1115_strict_qw2191_selector_premise_lattice_and_noncyclic_admissibility_checker.json"


class TestP2165StrictQW2191SelectorPremiseLatticeAndNoncyclicAdmissibilityChecker(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2165_s1115_strict_qw2191_selector_premise_lattice_and_noncyclic_admissibility_checker.py")],
            check=True,
        )
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema_version"], "p2165_s1115_v1")
        self.assertTrue(payload["gatekeeper_checks"]["lattice_checker_exported"])
        self.assertTrue(payload["strict_qw2191_selector_premise_lattice_and_noncyclic_admissibility_checker"]["checker_ready"])


if __name__ == "__main__":
    unittest.main()
