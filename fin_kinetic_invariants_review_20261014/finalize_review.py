#!/usr/bin/env python3
"""Mechanical review packaging; does not verify analytical weak-limit proofs."""
from pathlib import Path
import ast,hashlib,json,re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
PLAN=ROOT/"FIN_UNIFORM_CRITICAL_NEXT_CHAIN_20261014.md"
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    required=["README.md","INTAKE_20261014.md","UF_DECISION.md","CLAIM_REGISTER.json",
              "BENCHMARK_SCOPE.md","NONCONCLUSIONS.md","REPLAY.md","HANDOFF.md",
              "verify_intake.py","verification.json","finalize_review.py"]
    assert all((HERE/name).is_file() for name in required)
    parsed=[]
    def invalid(value):raise ValueError(value)
    for p in HERE.glob("*.json"):
        if p.name=="PACKAGE_CHECKS.json":continue
        json.loads(p.read_text(),parse_constant=invalid);parsed.append(p.name)
    for p in HERE.glob("*.py"):ast.parse(p.read_text(),filename=str(p))
    broken=[];whitespace=[]
    for p in [*HERE.glob("*.md"),PLAN]:
        for _,target in re.findall(r"\[([^\]]+)\]\(([^)]+)\)",p.read_text()):
            if target.startswith(("http://","https://","#")):continue
            if not (p.parent/target.split("#")[0]).exists():broken.append([p.name,target])
        for number,line in enumerate(p.read_text().splitlines(),1):
            if line!=line.rstrip():whitespace.append([p.name,number])
    assert not broken,broken
    assert not whitespace,whitespace
    cards=PLAN.read_text().split("### UF-")[1:];assert len(cards)==5
    fields=["Scientific question","Why it matters","Exact inputs","Dependencies","Method",
            "Deliverables","Acceptance criterion","Failure/refutation criterion","What becomes possible",
            "Must NOT conclude","Estimated computational class","Priority"]
    for i,card in enumerate(cards,1):
        assert card.startswith(f"{i:03d} ") and all(f in card for f in fields),i
    verification=json.loads((HERE/"verification.json").read_text())
    assert verification["status"]=="PASS_SCOPED_KI_INTAKE_CHECKS"
    assert len(verification["checks"])==46 and all(c["passed"] for c in verification["checks"])
    assert verification["replay"]["source_flags"]==[41,41]
    result={"status":"PASS_DOCUMENT_AND_PACKAGE_CHECKS","date":"2026-10-14",
            "required_files":required,"json_parsed":sorted(parsed),"python_syntax":"PASS",
            "broken_links":broken,"trailing_whitespace":whitespace,"task_card_contracts":5,
            "next_plan_sha256":sha(PLAN),"scope":"syntax and packaging, not KI/UF analytical proof"}
    (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(result,indent=2)+"\n")
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256")
    (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"manifest_entries":len(files),"task_cards":5,"status":result["status"]},indent=2))

if __name__=="__main__":main()
