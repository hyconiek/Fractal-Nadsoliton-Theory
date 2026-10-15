#!/usr/bin/env python3
"""Validate/package this DKO review and unexecuted NOC plan only."""
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
PLAN=ROOT/"FIN_NATIVE_OBSERVATION_COST_NEXT_CHAIN_20261015.md"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reject(x):raise ValueError(x)


def main():
    required=["INTAKE_20261015.md","PROOF_COMPLETIONS.md","SIGN_EVENT_CERTIFICATE.json","certify_sign_event.py",
              "verify_intake.py","verification.json","CLAIM_REGISTER.json","NOC_DECISION.md","NONCONCLUSIONS.md",
              "REPLAY.md","HANDOFF.md","README.md","finalize_review.py"]
    assert all((HERE/p).is_file() for p in required)
    parsed=[]
    for p in HERE.glob("*.json"):
        if p.name!="PACKAGE_CHECKS.json":json.loads(p.read_text(),parse_constant=reject); parsed.append(p.name)
    for p in HERE.glob("*.py"):ast.parse(p.read_text(),filename=str(p))
    broken=[]; spaces=[]
    for p in [*HERE.glob("*.md"),PLAN]:
        for target in re.findall(r"\[[^\]]+\]\(([^)]+)\)",p.read_text()):
            if target.startswith(("http://","https://","#")):continue
            if not (p.parent/target.split("#")[0]).exists():broken.append([p.name,target])
        for i,line in enumerate(p.read_text().splitlines(),1):
            if line!=line.rstrip():spaces.append([p.name,i])
    assert not broken,broken
    assert not spaces,spaces
    cards=PLAN.read_text().split("### NOC-")[1:]
    fields=["Title","Scientific question","Why it matters","Exact inputs","Dependencies","Method","Deliverables",
            "Acceptance criterion","Failure/refutation criterion","What becomes possible","Must NOT conclude","Estimated computational class","Priority"]
    assert len(cards)==5
    for i,card in enumerate(cards,1):assert card.startswith(f"{i:03d} ") and all(f in card for f in fields),i
    evidence=json.loads((HERE/"verification.json").read_text())
    assert evidence["status"]=="PASS_SCOPED_DKO_INTAKE_CHECKS" and len(evidence["checks"])==22 and all(c["passed"] for c in evidence["checks"])
    assert evidence["replay"]["source_checks"]==32 and evidence["replay"]["historical_cache_executed"] is False
    for rel,h in evidence["input_hashes"].items():assert sha(ROOT/rel)==h,rel
    from certify_sign_event import build
    assert build()==json.loads((HERE/"SIGN_EVENT_CERTIFICATE.json").read_text())==evidence["certificate"]
    result={"date":"2026-10-15","status":"PASS_DOCUMENT_AND_PACKAGE_CHECKS","required_files":required,"json_parsed":sorted(parsed),
            "python_syntax":"PASS","broken_links":broken,"trailing_whitespace":spaces,"NOC_cards":5,"next_plan_sha256":sha(PLAN),
            "scope":"syntax/links/saved inputs/event certificate/packaging; not mechanized DKO proof or NOC execution"}
    (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(result,indent=2)+"\n")
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256")
    (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"status":result["status"],"manifest_entries":len(files),"next_cards":5},indent=2))


if __name__=="__main__":main()
