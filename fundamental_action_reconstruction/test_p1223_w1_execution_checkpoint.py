#!/usr/bin/env python3
"""Regression tests for P1223 witness-evaluation semantics."""
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p1223_w1_execution_checkpoint.py"


class P1223CheckpointTests(unittest.TestCase):
    def _run_case(self, strict_theorem: bool, axiom: bool, tag: bool) -> dict:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            p1222 = tmp / "p1222.json"
            p1192 = tmp / "p1192.json"
            inp = tmp / "input.json"
            out = tmp / "out.json"

            p1222.write_text(
                json.dumps(
                    {
                        "packet": "P1222",
                        "selected_open_witness": {"id": "W1_selector_uniqueness_bridge"},
                    }
                ),
                encoding="utf-8",
            )
            p1192.write_text(
                json.dumps(
                    {
                        "packet": "P1192",
                        "witness_obligations": [
                            {
                                "id": "W1_selector_uniqueness_bridge",
                                "target": "QW-2191 strict-core uniqueness obstruction",
                                "pass_condition": "unique strict branch selected with reproducible witness",
                                "required_artifact": "formal selector source theorem OR explicit symmetry-breaking axiom with scoped non-strict tag",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )

            payload = {
                "packet": "P1223_INPUT",
                "as_of": "2026-05-11",
                "witness_id": "W1_selector_uniqueness_bridge",
                "available_inputs": {
                    "strict_selector_source_theorem_exported": strict_theorem,
                    "explicit_symmetry_breaking_axiom_declared": axiom,
                    "scoped_non_strict_tag_if_axiom_used": tag,
                },
                "guardrails": {
                    "qw_2191_is_real_obstruction": True,
                    "strict_closure_claim_allowed": False,
                },
            }
            inp.write_text(json.dumps(payload), encoding="utf-8")

            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--p1222",
                    str(p1222),
                    "--p1192",
                    str(p1192),
                    "--input",
                    str(inp),
                    "--out",
                    str(out),
                ],
                check=True,
                cwd=ROOT.parent,
            )

            return json.loads(out.read_text(encoding="utf-8"))

    def test_open_when_no_path_available(self) -> None:
        result = self._run_case(False, False, False)
        self.assertFalse(result["evaluation"]["pass_condition_met"])
        self.assertEqual(result["witness_status"], "OPEN")
        self.assertEqual(result["evaluation"]["discharge_mode"], "NO_DISCHARGE")
        self.assertFalse(result["evaluation"]["strict_core_closure_eligibility"])

    def test_discharged_when_strict_path_available(self) -> None:
        result = self._run_case(True, False, False)
        self.assertTrue(result["evaluation"]["strict_theorem_path_available"])
        self.assertTrue(result["evaluation"]["pass_condition_met"])
        self.assertEqual(result["witness_status"], "DISCHARGED")
        self.assertEqual(result["evaluation"]["discharge_mode"], "STRICT_PATH_DISCHARGE")
        self.assertTrue(result["evaluation"]["strict_core_closure_eligibility"])

    def test_discharged_when_tagged_axiom_path_available(self) -> None:
        result = self._run_case(False, True, True)
        self.assertTrue(result["evaluation"]["non_strict_axiom_path_available"])
        self.assertTrue(result["evaluation"]["pass_condition_met"])
        self.assertEqual(result["witness_status"], "DISCHARGED")
        self.assertEqual(result["evaluation"]["discharge_mode"], "NON_STRICT_AXIOM_DISCHARGE")
        self.assertFalse(result["evaluation"]["strict_core_closure_eligibility"])

    def test_strict_path_precedence_when_both_paths_available(self) -> None:
        result = self._run_case(True, True, True)
        self.assertEqual(result["evaluation"]["discharge_mode"], "STRICT_PATH_DISCHARGE")
        self.assertTrue(result["evaluation"]["strict_core_closure_eligibility"])


if __name__ == "__main__":
    unittest.main()
