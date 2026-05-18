from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1989_s939_strict_minimal_provider_family_no_go_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1989_s939_strict_minimal_provider_family_no_go_witness.json"

class TestP1989S939(unittest.TestCase):
    def test_minimal_provider_no_go(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P1989")
        self.assertEqual(p["result_kind"], "PASS_MINIMAL_PROVIDER_NO_GO_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["minimal_provider_has_no_exact_solution"])
        self.assertTrue(p["gatekeeper_checks"]["numeric_residual_nonzero"])

if __name__ == "__main__":
    unittest.main()
