from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1287_qw2191_r6_new_observable_provider_or_axiom_tagged_extension_checkpoint.py"


class TestP1287QW2191R6NewObservableProviderOrAxiomTaggedExtensionCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1286 = td_path / "p1286.json"
            out = td_path / "p1287.json"
            p1286.write_text(
                json.dumps({"next_priority": "R6_NEW_OBSERVABLE_PROVIDER_OR_AXIOM_TAGGED_EXTENSION"}),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1286", str(p1286), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r6_program"]["observable_provider"]["status"], "DECLARED")
            self.assertFalse(payload["r6_program"]["axiom_tagged_extension"]["enabled"])
            self.assertEqual(payload["next_priority"], "R7_PROVIDER_VALIDATION_AND_SELECTOR_SPLIT_TEST")

    def test_requires_priority(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1286 = td_path / "p1286.json"
            out = td_path / "p1287.json"
            p1286.write_text(json.dumps({"next_priority": "OTHER"}), encoding="utf-8")
            proc = subprocess.run(["python3", str(SCRIPT), "--p1286", str(p1286), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("R6_NEW_OBSERVABLE_PROVIDER_OR_AXIOM_TAGGED_EXTENSION", proc.stderr)


if __name__ == "__main__":
    unittest.main()
