#!/usr/bin/env python3
"""Mechanical review/plan validation and hashes; no new scientific campaign."""
import ast
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PLAN = ROOT / "FIN_RECORD_INFORMATION_NEXT_CHAIN_20261015.md"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def reject(value):
    raise ValueError(value)


def main():
    required = ["INTAKE_20261015.md", "PROOF_COMPLETIONS_AND_CORRECTIONS.md", "CORRECTION_REGISTER.json",
                "REVERSIBILITY_AND_RECORD_SCOPE.md", "CLAIM_REGISTER.json", "NONCONCLUSIONS.md",
                "RI_DECISION.md", "README.md", "HANDOFF.md", "REPLAY.md", "verification.json",
                "verify_intake.py", "observation_kernels.py", "finalize_review.py"]
    assert all((HERE/name).is_file() for name in required)
    parsed = []
    for p in HERE.glob("*.json"):
        if p.name != "PACKAGE_CHECKS.json":
            json.loads(p.read_text(), parse_constant=reject)
            parsed.append(p.name)
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text(), filename=str(p))
    bad_links, whitespace = [], []
    for p in [*HERE.glob("*.md"), PLAN]:
        for target in re.findall(r"\[[^\]]+\]\(([^)]+)\)", p.read_text()):
            if target.startswith(("https://", "http://", "#")):
                continue
            if not (p.parent/target.split("#")[0]).exists():
                bad_links.append([p.name, target])
        for i, line in enumerate(p.read_text().splitlines(), 1):
            if line != line.rstrip():
                whitespace.append([p.name, i])
    assert not bad_links, bad_links
    assert not whitespace, whitespace
    cards = PLAN.read_text().split("### RI-")[1:]
    fields = ["Title", "Scientific question", "Why it matters", "Exact inputs", "Dependencies", "Method",
              "Deliverables", "Acceptance criterion", "Failure/refutation criterion", "What becomes possible",
              "Must NOT conclude", "Estimated computational class", "Priority"]
    assert len(cards) == 5
    for i, card in enumerate(cards, 1):
        assert card.startswith(f"{i:03d} ") and all(f in card for f in fields), i
    result = json.loads((HERE/"verification.json").read_text())
    assert result["status"] == "PASS_SCOPED_OBS_INTAKE_CHECKS"
    assert len(result["checks"]) == 29 and all(c["passed"] for c in result["checks"])
    assert result["replay"]["source_flags"] == [81, 81]
    for rel, h in result["input_hashes"].items():
        assert sha(ROOT/rel) == h, rel
    out = {"date": "2026-10-15", "status": "PASS_DOCUMENT_AND_PACKAGE_CHECKS", "required_files": required,
           "json_parsed": sorted(parsed), "python_syntax": "PASS", "broken_local_links": bad_links,
           "trailing_whitespace": whitespace, "next_task_contracts": len(cards), "next_plan_sha256": sha(PLAN),
           "scope": "links, syntax, recorded inputs and packaging; not mechanized theorem proof or RI execution"}
    (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(out, indent=2)+"\n")
    files = sorted(p for p in HERE.iterdir() if p.is_file() and p.name != "MANIFEST.sha256")
    (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"status": out["status"], "manifest_entries": len(files), "next_cards": len(cards)}, indent=2))


if __name__ == "__main__":
    main()
