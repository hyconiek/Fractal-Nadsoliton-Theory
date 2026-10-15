#!/usr/bin/env python3
"""Validate/package only this RI review and its unexecuted LI plan."""
import ast
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PLAN = ROOT/"FIN_PARITY_RECORD_SUFFICIENCY_NEXT_CHAIN_20261015.md"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def reject(x):
    raise ValueError(x)


def main():
    required = ["INTAKE_20261015.md","PROOF_COMPLETIONS.md","TAG_FEEDBACK_REMOVAL_LEMMA.md",
                "LI_DECISION.md","CLAIM_REGISTER.json","CORRECTIONS_AND_GATES.json",
                "NONCONCLUSIONS.md","REPLAY.md","HANDOFF.md","README.md",
                "verification.json","verify_intake.py","finalize_review.py"]
    assert all((HERE/p).is_file() for p in required)
    parsed = []
    for p in HERE.glob("*.json"):
        if p.name != "PACKAGE_CHECKS.json":
            json.loads(p.read_text(),parse_constant=reject)
            parsed.append(p.name)
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text(),filename=str(p))
    broken,spaces = [],[]
    for p in [*HERE.glob("*.md"),PLAN]:
        for target in re.findall(r"\[[^\]]+\]\(([^)]+)\)",p.read_text()):
            if target.startswith(("https://","http://","#")):
                continue
            if not (p.parent/target.split("#")[0]).exists():
                broken.append([p.name,target])
        for i,line in enumerate(p.read_text().splitlines(),1):
            if line != line.rstrip():
                spaces.append([p.name,i])
    assert not broken,broken
    assert not spaces,spaces
    cards = PLAN.read_text().split("### LI-")[1:]
    fields = ["Title","Scientific question","Why it matters","Exact inputs","Dependencies","Method",
              "Deliverables","Acceptance criterion","Failure/refutation criterion","What becomes possible",
              "Must NOT conclude","Estimated computational class","Priority"]
    assert len(cards)==5
    for i,card in enumerate(cards,1):
        assert card.startswith(f"{i:03d} ") and all(f in card for f in fields),i
    result = json.loads((HERE/"verification.json").read_text())
    assert result["status"]=="PASS_SCOPED_RI_INTAKE_CHECKS"
    assert len(result["checks"])==24 and all(c["passed"] for c in result["checks"])
    assert result["replay"]["source_flags"]==[5626,5626]
    for rel,h in result["input_hashes"].items():
        assert sha(ROOT/rel)==h,rel
    package = {"date":"2026-10-15","status":"PASS_DOCUMENT_AND_PACKAGE_CHECKS","required_files":required,
               "json_parsed":sorted(parsed),"python_syntax":"PASS","broken_links":broken,"trailing_whitespace":spaces,
               "next_cards":len(cards),"next_plan_sha256":sha(PLAN),
               "scope":"syntax, links, saved input hashes and packaging; not formal theorem verification or LI execution"}
    (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(package,indent=2)+"\n")
    files = sorted(p for p in HERE.iterdir() if p.is_file() and p.name != "MANIFEST.sha256")
    (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"status":package["status"],"manifest_entries":len(files),"next_cards":len(cards)},indent=2))


if __name__ == "__main__":
    main()
