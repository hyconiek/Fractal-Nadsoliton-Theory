#!/usr/bin/env python3
"""Mechanical MT intake packaging. Does not execute new observation research."""
from pathlib import Path
import ast,hashlib,json,re
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
PLAN=ROOT/"FIN_INTRINSIC_OBSERVATION_NEXT_CHAIN_20261015.md"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    required=["README.md","INTAKE_20261015.md","PROOF_COMPLETIONS.md","LOCAL_GLOBAL_SEPARATION.md",
              "OBS_DECISION.md","CLAIM_REGISTER.json","NONCONCLUSIONS.md","REPLAY.md","HANDOFF.md",
              "GEOMETRY_CERTIFICATE.json","check_geometry.py","verify_intake.py","verification.json","finalize_review.py"]
    assert all((HERE/p).is_file() for p in required)
    def reject(x):raise ValueError(x)
    parsed=[]
    for p in HERE.glob("*.json"):
        if p.name=="PACKAGE_CHECKS.json":continue
        json.loads(p.read_text(),parse_constant=reject);parsed.append(p.name)
    for p in HERE.glob("*.py"):ast.parse(p.read_text(),filename=str(p))
    broken=[];spaces=[]
    for p in [*HERE.glob("*.md"),PLAN]:
        for _,target in re.findall(r"\[([^\]]+)\]\(([^)]+)\)",p.read_text()):
            if target.startswith(("http://","https://","#")):continue
            if not (p.parent/target.split("#")[0]).exists():broken.append([p.name,target])
        for num,line in enumerate(p.read_text().splitlines(),1):
            if line!=line.rstrip():spaces.append([p.name,num])
    assert not broken,broken
    assert not spaces,spaces
    cards=PLAN.read_text().split("### OBS-C")[1:];assert len(cards)==5
    fields=["Scientific question","Why it matters","Exact inputs","Dependencies","Method","Deliverables",
            "Acceptance criterion","Failure/refutation criterion","What becomes possible","Must NOT conclude","Estimated computational class","Priority"]
    for i,card in enumerate(cards,1):assert card.startswith(f"{i:02d} ") and all(f in card for f in fields),i
    checked=json.loads((HERE/"verification.json").read_text())
    assert len(checked["checks"])==24 and all(c["passed"] for c in checked["checks"])
    assert checked["status"]=="PASS_SCOPED_MT_INTAKE_CHECKS"
    assert checked["replay"]["producer_flags"]==[57,57] and checked["replay"]["independent_flags"]==[58,58]
    from check_geometry import build
    assert json.loads((HERE/"GEOMETRY_CERTIFICATE.json").read_text())==build()
    package={"date":"2026-10-15","status":"PASS_DOCUMENT_AND_PACKAGE_CHECKS","required_files":required,
             "json_parsed":sorted(parsed),"python_syntax":"PASS","broken_local_links":broken,
             "trailing_whitespace":spaces,"next_task_contracts":5,"next_plan_sha256":sha(PLAN),
             "scope":"syntax, payments, links and packaging; not formal MT proof or executed OBS theorem"}
    (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(package,indent=2)+"\n")
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256")
    (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"review_entries":len(files),"next_cards":5,"status":package["status"]},indent=2))

if __name__=="__main__":main()
