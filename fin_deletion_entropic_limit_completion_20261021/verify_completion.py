#!/usr/bin/env python3
"""Read-only completion of the now-committed DEL packet; never reseal it."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import tempfile

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/"fin_deletion_entropic_limit_research_20261020"
COMMIT="a3397ccad56c4d74c7c676a066876ac271783b62"
SOURCE_MANIFEST="607d32da2774029780b36699ebdc55e6c94abead40429af3bd08f09355dad49f"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha(p):
    return digest(p.read_bytes())


def snapshot(folder):
    return {p.relative_to(folder).as_posix():sha(p) for p in sorted(folder.rglob("*"))
            if p.is_file() and "__pycache__" not in p.parts and p.suffix!=".pyc"}


def manifest(folder):
    seen=set()
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        h,name=line.split("  ",1)
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name)==h,(folder,name)
        seen.add(name)
    assert set(snapshot(folder))==seen|{"MANIFEST.sha256"},folder
    return len(seen)


def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(55,55))
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))


def replay_source():
    with tempfile.TemporaryDirectory(prefix="fin-del-completion-") as temp:
        output=Path(temp)/"read_only_result.json"
        child=subprocess.run([sys.executable,str(SOURCE/"verify_research.py"),"--check","--output",str(output)],
            cwd=ROOT,env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1"),
            capture_output=True,text=True,timeout=55,preexec_fn=limits)
        assert child.returncode==0,(child.stdout,child.stderr)
        result=json.loads(output.read_text())
    assert result==json.loads((SOURCE/"RESULTS.json").read_text())
    assert result["bounded_checks"]==len(result["checks"])==41
    assert len(result["input_hashes"])==17 and len(result["protected_manifests"])==9
    assert sum(v["entries"] for v in result["protected_manifests"].values())==182
    return result,{"returncode":child.returncode,"stdout":child.stdout,"stderr":child.stderr,
                  "mode":"existing source --check; isolated temporary output; not --record"}


def linkcheck():
    for p in HERE.glob("*.md"):
        lines=[]
        fence=False
        for line in p.read_text().splitlines():
            assert line==line.rstrip(),p.name
            if line.lstrip().startswith("```"):
                fence=not fence
            elif not fence and not line.startswith(("    ","\t")):
                lines.append(line)
        txt=re.sub(r"`[^`\n]*`","","\n".join(lines))
        for target in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)",txt):
            target=target.strip("<>").split("#")[0]
            if target and not target.startswith(("http:","https:")):
                assert (p.parent/target).exists(),(p.name,target)


def main():
    ap=argparse.ArgumentParser()
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--record",action="store_true")
    group.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    if args.check:
        manifest(HERE)
    assert sha(SOURCE/"MANIFEST.sha256")==SOURCE_MANIFEST
    entries=manifest(SOURCE)
    before=snapshot(SOURCE)
    assert entries==18 and len(before)==19
    for name,h in before.items():
        relative=f"{SOURCE.name}/{name}"
        blob=subprocess.check_output(["git","show",f"{COMMIT}:{relative}"],cwd=ROOT,timeout=10)
        assert digest(blob)==h,(name,"source differs from committed bytes")
    source,log=replay_source()
    assert before==snapshot(SOURCE)
    input_hashes=dict(source["input_hashes"])
    input_hashes.update({f"{SOURCE.name}/{n}":h for n,h in before.items()})
    result={"verification_date":"2026-10-21","source_start_date":"2026-10-20",
            "status":"SCOPED_DEL_COMPLETION_SOURCE_UNCHANGED_IBA_CLARIFIED",
            "entry_and_source_commit":COMMIT,"source_manifest_sha256":SOURCE_MANIFEST,
            "source_manifest_entries":entries,"source_present_files":len(before),
            "source_files_match_local_commit":True,"source_before_and_after":before,
            "source_replay":source,"replay_log":log,"bound_inputs":input_hashes,
            "specific_deletion_comparison_IBA":"HOLDS_BY_STATIONARITY_AND_SOURCE_NC",
            "general_native_KR_IBA":"OPEN","F01":"OPEN_SOURCE",
            "new_theorem_count":0,"same_author":True,"independent_peer_review":False,
            "AGENTS_or_originals_modified":False,"commit_created_by_this_assistant":False,
            "broad_goal":"ACTIVE_NOT_COMPLETE"}
    intext="".join(f"{h}  {n}\n" for n,h in sorted(input_hashes.items()))
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result,indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(intext)
        linkcheck()
        files=sorted(p for p in HERE.rglob("*") if p.is_file() and p.name!="MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        assert json.loads((HERE/"RESULTS.json").read_text())==result
        assert (HERE/"INPUTS.sha256").read_text()==intext
        linkcheck()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"verification_date":result["verification_date"],
                      "source_files_preserved":len(before),"source_checks":41,"bound_inputs":len(input_hashes),
                      "specific_comparison_IBA":"HOLDS","general_KR_IBA":"OPEN","F01":"OPEN_SOURCE",
                      "source_or_AGENTS_modified":False,"independent_review":False},indent=2))


if __name__=="__main__":
    main()
