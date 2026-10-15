#!/usr/bin/env python3
"""Mechanical review/activation-plan checks. --check never writes."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
from certify_angular_gap import build

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
PLAN=ROOT/"FIN_POST_SRN_BRIDGE_ACTIVATION_CHAIN_20261015.md"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def reject(value):
    raise ValueError(value)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    required=["INTAKE_20261015.md","PROOF_COMPLETIONS.md","READOUT_ERROR_SCOPE.md",
              "ANGULAR_GAP_CERTIFICATE.json","certify_angular_gap.py","verify_intake.py",
              "verification.json","BRIDGE_DECISION.md","CLAIM_REGISTER.json","NONCONCLUSIONS.md",
              "README.md","HANDOFF.md","REPLAY.md","finalize_review.py"]
    assert all((HERE/p).is_file() for p in required)
    parsed=[]
    for p in HERE.glob("*.json"):
        if p.name!="PACKAGE_CHECKS.json":
            json.loads(p.read_text(),parse_constant=reject)
            parsed.append(p.name)
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text(),filename=str(p))
    broken=[]; whitespace=[]
    for p in [*HERE.glob("*.md"),PLAN]:
        for target in re.findall(r"\[[^\]]+\]\(([^)]+)\)",p.read_text()):
            if target.startswith(("http://","https://","#")):
                continue
            if not (p.parent/target.split("#")[0]).exists():
                broken.append([p.name,target])
        for n,line in enumerate(p.read_text().splitlines(),1):
            if line!=line.rstrip():
                whitespace.append([p.name,n])
    assert not broken,broken
    assert not whitespace,whitespace
    fields=["Title","Scientific question","Why it matters","Exact inputs","Dependencies","Method",
            "Deliverables","Acceptance criterion","Failure/refutation criterion","What becomes possible",
            "Must NOT conclude","Estimated computational class","Priority"]
    cards=PLAN.read_text().split("### BR-")[1:]
    assert len(cards)==3
    for i,card in enumerate(cards,1):
        assert card.startswith(f"{i:03d} ") and all(f in card for f in fields),i
    evidence=json.loads((HERE/"verification.json").read_text())
    assert evidence["status"]=="PASS_SCOPED_SRN_INTAKE_CHECKS"
    assert len(evidence["checks"])==28 and all(c["passed"] for c in evidence["checks"])
    assert evidence["source_manifest"]["matched"]==60 and evidence["replay"]["source_checks"]==535
    for rel,h in evidence["input_hashes"].items():
        assert sha(ROOT/rel)==h,rel
    assert build()==json.loads((HERE/"ANGULAR_GAP_CERTIFICATE.json").read_text())==evidence["ideal_certificate"]
    result={"date":"2026-10-15","status":"PASS_DOCUMENT_AND_PACKAGE_CHECKS","required_files":required,
            "json_parsed":sorted(parsed),"python_syntax":"PASS","broken_links":broken,"trailing_whitespace":whitespace,
            "activation_cards":3,"next_plan_sha256":sha(PLAN),
            "scope":"syntax, links, saved input hashes, ideal arithmetic, packaging; not mechanized native proof or BR execution"}
    if args.check:
        assert result==json.loads((HERE/"PACKAGE_CHECKS.json").read_text())
        names=[]
        for row in (HERE/"MANIFEST.sha256").read_text().splitlines():
            h,name=row.split("  ",1)
            assert name not in names and sha(HERE/name)==h,name
            names.append(name)
        actual={p.name for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256"}
        assert set(names)==actual
    else:
        (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(result,indent=2)+"\n")
        files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256")
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"status":result["status"],"read_only":args.check,"manifest_entries":len(required)+1,"activation_cards":len(cards)},indent=2))

if __name__=="__main__":
    main()
