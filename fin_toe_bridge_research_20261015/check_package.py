#!/usr/bin/env python3
"""Research package integrity, links, input/capture hashes and task cards.

--write records metadata/manifest for this new package only. Default is read-only.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
from check_locality_certificate import verify

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
PLAN=ROOT/"FIN_TOE_THEORETICAL_NEXT_CHAIN_20261015.md"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reject(value):raise ValueError(value)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--write",action="store_true"); args=parser.parse_args()
    required=["README.md","RESEARCH_SYNTHESIS_PL.md","DERIVED_RESULTS.md","TRANSFER_PROOF_REVIEW.md",
              "REPOSITORY_FRONTIER_AUDIT.md","LITERATURE_DATABASE.md","LITERATURE_REGISTER.json",
              "CLAIM_REGISTER.json","PRIORITY_TASKS.json","AGENTS_PROPOSED_PATCH.md","NONCONCLUSIONS.md",
              "REPLAY.md","HANDOFF.md","verify_bridges.py","check_locality_certificate.py","check_package.py","verification.json"]
    assert all((HERE/name).is_file() for name in required)
    for p in HERE.glob("*.json"):json.loads(p.read_text(),parse_constant=reject)
    for p in HERE.glob("*.py"):ast.parse(p.read_text(),filename=str(p))
    for p in [*HERE.glob("*.md"),PLAN]:
        assert all(line==line.rstrip() for line in p.read_text().splitlines()),p
        for target in re.findall(r"\[[^\]]+\]\(([^)]+)\)",p.read_text()):
            if not target.startswith(("https://","http://","#")):
                assert (p.parent/target.split("#",1)[0]).exists(),(p.name,target)
    evidence=json.loads((HERE/"verification.json").read_text()); verify(evidence["certificate"])
    assert len(evidence["checks"])==10 and all(c["passed"] for c in evidence["checks"])
    for rel,h in evidence["input_hashes"].items():assert sha(ROOT/rel)==h,rel
    sources=json.loads((HERE/"LITERATURE_REGISTER.json").read_text())
    assert len(sources["sources"])==sources["source_count"]==15
    missing=[]; captures=0
    for record in [*sources["sources"],*sources["unread_originals"]]:
        p=ROOT/record["cache_path"]
        if p.is_file():assert sha(p)==record["cache_sha256"],p; captures+=1
        else:missing.append(record["cache_path"])
    cards=PLAN.read_text().split("### NT-")[1:]
    fields=["Title","Scientific question","Why it matters","Exact inputs","Dependencies","Method",
            "Deliverables","Acceptance criterion","Failure/refutation criterion","What becomes possible",
            "Must NOT conclude","Estimated computational class","Priority"]
    assert len(cards)==6
    for i,card in enumerate(cards,1):assert card.startswith(f"{i:03d} ") and all(f in card for f in fields),i
    tasklist=json.loads((HERE/"PRIORITY_TASKS.json").read_text())
    assert [t["id"] for t in tasklist["tasks"]]==[f"NT-{i:03d}" for i in range(1,7)]
    status={"date":"2026-10-15","status":"PASS_RESEARCH_PACKAGE_CHECKS","required_files":required,
            "tasks":6,"source_records":15,"unread_scan_records":1,"plan_sha256":sha(PLAN),
            "checks":"JSON/syntax/links/whitespace/scientific input hashes/exact rational witnesses/task-card fields",
            "scope":"research submission, not formal proof, physical evidence, historical custody or NT campaign execution"}
    if args.write:
        (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(status,indent=2)+"\n")
        files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256")
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    else:
        assert status==json.loads((HERE/"PACKAGE_CHECKS.json").read_text())
        names=[]
        for row in (HERE/"MANIFEST.sha256").read_text().splitlines():
            h,name=row.split("  ",1); assert name not in names and sha(HERE/name)==h,name; names.append(name)
        assert set(names)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256"}
    print(json.dumps({"status":status["status"],"read_only":not args.write,"captures_matched":captures,
                      "optional_captures_absent":missing,"tasks":6,"research_checks":10},indent=2))

if __name__=="__main__":main()
