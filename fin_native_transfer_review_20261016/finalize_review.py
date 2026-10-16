#!/usr/bin/env python3
"""Mechanical NT intake/HM plan validation, not a proof assistant."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
PLAN=ROOT/"FIN_AUXILIARY_MEDIATOR_NEXT_CHAIN_20261016.md"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reject(x):raise ValueError(x)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--check",action="store_true"); args=parser.parse_args()
    required=["INTAKE_20261016.md","PROOF_COMPLETIONS_AND_SCOPE.md","MEDIATOR_BRIDGE_SEED.md",
              "HM_DECISION.md","CLAIM_REGISTER.json","NONCONCLUSIONS.md","README.md","HANDOFF.md",
              "REPLAY.md","verify_intake.py","verification.json","finalize_review.py"]
    assert all((HERE/name).is_file() for name in required)
    for p in HERE.glob("*.json"):json.loads(p.read_text(),parse_constant=reject)
    for p in HERE.glob("*.py"):ast.parse(p.read_text(),filename=str(p))
    for p in [*HERE.glob("*.md"),PLAN]:
        assert all(line==line.rstrip() for line in p.read_text().splitlines()),p
        for target in re.findall(r"\[[^\]]+\]\(([^)]+)\)",p.read_text()):
            if not target.startswith(("http://","https://","#")):
                assert (p.parent/target.split("#",1)[0]).exists(),(p.name,target)
    saved=json.loads((HERE/"verification.json").read_text())
    assert saved["status"]=="PASS_SCOPED_NT_INTAKE_CHECKS"
    assert len(saved["checks"])==18 and all(row["passed"] for row in saved["checks"])
    assert saved["source_manifest"]["matched"]==67 and saved["replay"]["bounded_source_flags"]==14
    for rel,h in saved["input_hashes"].items():assert sha(ROOT/rel)==h,rel
    cards=PLAN.read_text().split("### HM-")[1:]
    fields=["Title","Scientific question","Why it matters","Exact inputs","Dependencies","Method",
            "Deliverables","Acceptance criterion","Failure/refutation criterion","What becomes possible",
            "Must NOT conclude","Estimated computational class","Priority"]
    assert len(cards)==5
    for i,card in enumerate(cards,1):assert card.startswith(f"{i:03d} ") and all(field in card for field in fields),i
    result={"date":"2026-10-16","status":"PASS_REVIEW_PACKAGE_CHECKS","required_files":required,
            "HM_cards":5,"plan_sha256":sha(PLAN),"syntax_json_links_whitespace":"PASS",
            "scope":"mechanical intake/plan/hash checks, not formal NT proof or HM execution"}
    if args.check:
        assert result==json.loads((HERE/"PACKAGE_CHECKS.json").read_text())
        names=[]
        for row in (HERE/"MANIFEST.sha256").read_text().splitlines():
            h,name=row.split("  ",1); assert name not in names and sha(HERE/name)==h,name; names.append(name)
        assert set(names)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256"}
    else:
        (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(result,indent=2)+"\n")
        files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256")
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"status":result["status"],"read_only":args.check,"manifest_entries":len(required)+1,"HM_cards":5},indent=2))

if __name__=="__main__":main()
