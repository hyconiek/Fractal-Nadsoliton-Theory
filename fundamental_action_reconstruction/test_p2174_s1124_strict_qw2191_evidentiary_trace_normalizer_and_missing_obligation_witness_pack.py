from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2174S1124StrictQW2191TraceNormalizer(unittest.TestCase):
    def test_trace_normalizer(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2174_s1124_strict_qw2191_evidentiary_trace_normalizer_and_missing_obligation_witness_pack.py")],
            check=True,
        )
        d = json.loads(
            (G / "p2174_s1124_strict_qw2191_evidentiary_trace_normalizer_and_missing_obligation_witness_pack.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(d["schema_version"], "p2174_s1124_v1")
        self.assertTrue(d["gatekeeper_checks"]["trace_normalization_exported"])
        traces = d["strict_qw2191_evidentiary_trace_normalizer_and_missing_obligation_witness_pack"]["normalized_obligation_traces"]
        self.assertTrue(all("evidence_items" in t and "premise_items" in t for t in traces))


if __name__ == "__main__":
    unittest.main()
