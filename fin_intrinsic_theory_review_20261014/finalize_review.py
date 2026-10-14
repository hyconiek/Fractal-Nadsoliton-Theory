#!/usr/bin/env python3
"""Mechanical review packaging; scientific checks are in verify_intake.py."""
from pathlib import Path
import ast
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PLAN = ROOT / "FIN_STOCHASTIC_FOLD_CLOSURE_NEXT_CHAIN_20261014.md"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    required = ["README.md", "HANDOFF.md", "INTAKE_20261014.md", "CLAIM_REGISTER.json",
                "TH_010_DECISION.md", "FAST_MODE_PROOF_ROUTE.md", "NONCONCLUSIONS.md",
                "REPLAY.md", "verification.json", "verify_intake.py", "finalize_review.py"]
    for name in required:
        assert (HERE/name).is_file(), name
    assert PLAN.is_file()
    parsed = []
    for path in sorted(HERE.glob("*.json")):
        if path.name == "PACKAGE_CHECKS.json":
            continue
        json.loads(path.read_text(), parse_constant=lambda val: (_ for _ in ()).throw(ValueError(val)))
        parsed.append(path.name)
    for path in HERE.glob("*.py"):
        ast.parse(path.read_text(), filename=str(path))
    broken, whitespace = [], []
    for path in [*HERE.glob("*.md"), PLAN]:
        for label, target in re.findall(r"\[([^\]]+)\]\(([^)]+)\)", path.read_text()):
            if target.startswith(("https://", "http://", "#")):
                continue
            dest = path.parent/target.split("#")[0]
            if not dest.exists():
                broken.append({"file":str(path.relative_to(ROOT)),"target":target})
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if line != line.rstrip():
                whitespace.append([str(path.relative_to(ROOT)),number])
    assert not broken, broken
    assert not whitespace, whitespace
    verification = json.loads((HERE/"verification.json").read_text())
    assert verification["status"] == "PASS_SCOPED_INTAKE_CHECKS"
    assert verification["replay"]["independent_checks"] == [41,41]
    assert len(verification["checks"]) == 22 and all(r["passed"] for r in verification["checks"])
    text = PLAN.read_text()
    for task in range(1,6):
        assert f"### SF-{task:03d} " in text
    checks = {"status":"PASS_DOCUMENT_AND_PACKAGE_CHECKS", "date":"2026-10-14",
              "json_parsed":parsed, "python_syntax":"PASS", "broken_local_links":broken,
              "trailing_whitespace":whitespace, "required_files":required,
              "next_plan_sha256":digest(PLAN), "next_tasks":[f"SF-{i:03d}" for i in range(1,6)],
              "scope":"packaging/syntax/link checks, not analytical or stochastic proof"}
    (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(checks,indent=2)+"\n")
    files = sorted(p for p in HERE.iterdir() if p.is_file() and p.name != "MANIFEST.sha256")
    (HERE/"MANIFEST.sha256").write_text("".join(f"{digest(p)}  {p.name}\n" for p in files))
    print(json.dumps({"review_manifest_entries":len(files),"json_files":parsed,"status":checks["status"]},indent=2))


if __name__ == "__main__":
    main()
