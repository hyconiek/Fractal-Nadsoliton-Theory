#!/usr/bin/env python3
"""Mechanical LI review/KO-plan validation; not new scientific execution."""
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
PLAN=ROOT/"FIN_KERNEL_OBSERVABILITY_NEXT_CHAIN_20261015.md"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reject(x):raise ValueError(x)


def main():
    names=["INTAKE_20261015.md","PROOF_COMPLETIONS.md","CLAIM_REGISTER.json","NONCONCLUSIONS.md",
           "KO_DECISION.md","KERNEL_FAMILY_SEED.md","KERNEL_FAMILY_SEED.json","kernel_family_certificate.py",
           "verify_intake.py","verification.json","README.md","HANDOFF.md","REPLAY.md","finalize_review.py"]
    assert all((HERE/n).is_file() for n in names)
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
    cards=PLAN.read_text().split("### KO-")[1:]
    fields=["Title","Scientific question","Why it matters","Exact inputs","Dependencies","Method","Deliverables",
            "Acceptance criterion","Failure/refutation criterion","What becomes possible","Must NOT conclude","Estimated computational class","Priority"]
    assert len(cards)==5
    for i,card in enumerate(cards,1):assert card.startswith(f"{i:03d} ") and all(f in card for f in fields),i
    evidence=json.loads((HERE/"verification.json").read_text())
    assert evidence["status"]=="PASS_SCOPED_LI_INTAKE_CHECKS" and len(evidence["checks"])==23 and all(c["passed"] for c in evidence["checks"])
    assert evidence["replay"]["status"]=="PASS" and evidence["replay"]["source_assertions"]==20527
    for rel,h in evidence["input_hashes"].items():assert sha(ROOT/rel)==h,rel
    from kernel_family_certificate import build
    assert json.loads((HERE/"KERNEL_FAMILY_SEED.json").read_text())==build()==evidence["kernel_family_seed"]
    result={"date":"2026-10-15","status":"PASS_DOCUMENT_AND_PACKAGE_CHECKS","required_files":names,
            "json_parsed":sorted(parsed),"python_syntax":"PASS","broken_links":broken,"trailing_whitespace":spaces,
            "KO_cards":5,"next_plan_sha256":sha(PLAN),"scope":"syntax/links/inputs/exact seed/packaging, not mechanized LI proof or KO execution"}
    (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(result,indent=2)+"\n")
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256")
    (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"status":result["status"],"manifest_entries":len(files),"next_cards":5},indent=2))


if __name__=="__main__":main()
