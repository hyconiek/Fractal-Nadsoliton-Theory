from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2255S1205StrictNuBranchGroupPolicyDominanceEnvelopeScanProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2255_s1205_strict_nu_branch_group_policy_dominance_envelope_scan_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2255_s1205_strict_nu_branch_group_policy_dominance_envelope_scan_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2255_s1205_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["dominance_envelope_exported"])
        self.assertTrue(g["dominance_fraction_bounded"])
        self.assertTrue(g["nonempty_scan_grid"])


if __name__ == "__main__":
    unittest.main()
