#!/usr/bin/env python3
"""Mechanical UF review packaging; analytical MT tasks are not run."""
from pathlib import Path
import ast,hashlib,json,re
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
PLAN=ROOT/"FIN_CRITICAL_METASTABILITY_NEXT_CHAIN_20261014.md"
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    required=["README.md","INTAKE_20261014.md","CLAIM_REGISTER.json","GIBBS_EXCLUSION_PROOF.md",
              "GIBBS_EXCLUSION_CERTIFICATE.json","check_gibbs_certificate.py","verification.json",
              "verify_intake.py","MT_DECISION.md","TRAPPING_PROOF_ROUTE.md","NONCONCLUSIONS.md",
              "REPLAY.md","HANDOFF.md","finalize_review.py"]
    assert all((HERE/p).is_file() for p in required)
    parsed=[]
    def reject(value):raise ValueError(value)
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
    cards=PLAN.read_text().split("### MT-")[1:];assert len(cards)==5
    fields=["Scientific question","Why it matters","Exact inputs","Dependencies","Method","Deliverables",
            "Acceptance criterion","Failure/refutation criterion","What becomes possible","Must NOT conclude",
            "Estimated computational class","Priority"]
    for i,card in enumerate(cards,1):
        assert card.startswith(f"{i:03d} ") and all(f in card for f in fields),i
    result=json.loads((HERE/"verification.json").read_text())
    assert result["status"]=="PASS_SCOPED_UF_ARITHMETIC_AND_REPLAY"
    assert len(result["checks"])==40 and all(c["passed"] for c in result["checks"])
    assert result["replay"]["source_flags"]==[44,44]
    from check_gibbs_certificate import build,check_public
    check_public(json.loads((HERE/"GIBBS_EXCLUSION_CERTIFICATE.json").read_text()),build())
    package={"status":"PASS_DOCUMENT_AND_PACKAGE_CHECKS","date":"2026-10-14","required_files":required,
             "json_parsed":sorted(parsed),"python_syntax":"PASS","broken_local_links":broken,
             "trailing_whitespace":spaces,"task_card_contracts":5,"next_plan_sha256":sha(PLAN),
             "scope":"syntax, exact certificate payments and packaging; not an MT theorem or formal UF proof"}
    (HERE/"PACKAGE_CHECKS.json").write_text(json.dumps(package,indent=2)+"\n")
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!="MANIFEST.sha256")
    (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"review_entries":len(files),"task_cards":5,"status":package["status"]},indent=2))

if __name__=="__main__":main()
