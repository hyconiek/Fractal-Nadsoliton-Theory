#!/usr/bin/env python3
"""Mechanical review packaging; not an analytical proof checker."""
from pathlib import Path
import ast
import hashlib
import json
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
PLAN=ROOT/"FIN_KINETIC_INVARIANTS_NEXT_CHAIN_20261014.md"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    required=["README.md","INTAKE_20261014.md","STOPPED_LIMIT_COMPLETION.md","KI_DECISION.md",
              "CLAIM_REGISTER.json","NONCONCLUSIONS.md","REPLAY.md","HANDOFF.md",
              "verify_intake.py","verification.json","finalize_review.py"]
    assert all((HERE/p).is_file() for p in required)
    parsed=[]
    def invalid_constant(value):
        raise ValueError(value)
    for path in HERE.glob("*.json"):
        if path.name=="PACKAGE_CHECKS.json":continue
        json.loads(path.read_text(),parse_constant=invalid_constant);parsed.append(path.name)
    for path in HERE.glob("*.py"):
        ast.parse(path.read_text(),filename=str(path))
    broken=[];whitespace=[]
    for path in [*HERE.glob("*.md"),PLAN]:
        for _,target in re.findall(r"\[([^\]]+)\]\(([^)]+)\)",path.read_text()):
            if target.startswith(("http://","https://","#")):continue
            if not (path.parent/target.split("#")[0]).exists():broken.append([path.name,target])
        for number,line in enumerate(path.read_text().splitlines(),1):
            if line!=line.rstrip():whitespace.append([path.name,number])
    assert not broken,broken
    assert not whitespace,whitespace
    verification=json.loads((HERE/"verification.json").read_text())
    assert verification["status"]=="PASS_SCOPED_ARITHMETIC_AND_REPLAY"
    assert len(verification["checks"])==26 and all(x["passed"] for x in verification["checks"])
    assert verification["replay"]["source_checks"]==[41,41]
    cards=PLAN.read_text().split("### KI-")[1:]
    assert len(cards)==5
    fields=["Scientific question","Why it matters","Exact inputs","Dependencies","Method",
            "Deliverables","Acceptance criterion","Failure/refutation criterion","What becomes possible",
            "Must NOT conclude","Estimated computational class","Priority"]
    for index,card in enumerate(cards,1):
        assert card.startswith(f"{index:03d} ")
        assert all(field in card for field in fields),(index,fields)
    checks={"status":"PASS_DOCUMENT_AND_PACKAGE_CHECKS","date":"2026-10-14",
            "required_files":required,"json_parsed":sorted(parsed),"python_syntax":"PASS",
            "broken_local_links":broken,"trailing_whitespace":whitespace,
            "next_plan_sha256":sha(PLAN),"task_card_contracts":5,
            "scope":"packaging and syntax, not stochastic or cross-kinetic proof"}
    (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(checks,indent=2)+"\n")
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256")
    (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"review_manifest_entries":len(files),"task_cards":len(cards),"status":checks["status"]},indent=2))

if __name__=="__main__":main()
