#!/usr/bin/env python3
"""Exact scope controls and dependency replay; imported math is not mechanized."""
import argparse
import ast
from collections import defaultdict
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BI = ROOT/"fin_critical_bit_information_budget_20261021"
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "fin_rank7_followup/src/model.py",
    "FIN_POST_AFTER_CONTINUATION_FULL_HANDOFF_20260926/01_REPORTS/FISHER_KINETIC_METRIC_44.md",
    "FIN_ST2532_ST2621_Configuration_Generator_and_Dynamic_RG_Report_EN.md",
    "fin_critical_bit_information_budget_20261021/PROOF_OR_OPEN_LEMMA.md",
    "fin_critical_bit_information_budget_20261021/FINITE_RECORD_AND_INITIAL_BUDGET.md",
    "fin_critical_bit_information_budget_20261021/RESULTS.json",
    "fin_critical_bit_information_budget_20261021/verify_research.py",
    "fin_unified_strict_evolution_audit_20261017/STRICT_EVOLUTION_SOURCE_AUDIT.md",
    "fin_coupled_evolution_source_review_20261017/PROOF_AUDIT_AND_COMPLETIONS.md",
    "fin_nonlinear_stability_transfer_review_20261017/PROOF_COMPLETIONS_AND_SCOPE.md",
]
PROTECTED = ["fin_critical_bit_information_budget_20261021",
             "fin_coupled_evolution_source_review_20261017",
             "fin_nonlinear_stability_transfer_review_20261017"]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder, coverage=False):
    seen = set()
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        h, name = line.split("  ", 1)
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name) == h, (folder, name)
        seen.add(name)
    if coverage:
        actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts}
        assert actual == seen | {"MANIFEST.sha256"}
    return {"entries": len(seen), "sha256": sha(folder/"MANIFEST.sha256"),
            "coverage": "complete_noncache" if coverage else "listed_members_only"}


def word(x, queries):
    return tuple(int(a*x+b > 0) for a, b in queries)


def all_affine_words(queries):
    roots = sorted({-b/a for a, b in queries if a})
    if roots:
        samples = roots+[roots[0]-1, roots[-1]+1]
        samples += [(a+b)/2 for a, b in zip(roots, roots[1:])]
    else:
        samples = [F(0)]
    return {word(x, queries) for x in samples}, samples


def parity(weights):
    assert len(weights) == 12 and min(weights) > 0
    normalizer = sum(weights)
    mean = sum(((-1)**i)*weights[i] for i in range(12))/normalizer
    return mean, int(mean > 0)


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (55, 55))


def replay_bi():
    before = {p.relative_to(BI).as_posix(): sha(p) for p in BI.rglob("*")
              if p.is_file() and "__pycache__" not in p.parts}
    with tempfile.TemporaryDirectory(prefix="fin-ls-bi-") as tmp:
        result = Path(tmp)/"check.json"
        run = subprocess.run([sys.executable, str(BI/"verify_research.py"), "--check", "--output", str(result)],
                             cwd=ROOT, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
                             capture_output=True, text=True, timeout=55, preexec_fn=limits)
        assert run.returncode == 0, (run.stdout, run.stderr)
        out = json.loads(result.read_text())
    assert out == json.loads((BI/"RESULTS.json").read_text())
    assert before == {p.relative_to(BI).as_posix(): sha(p) for p in BI.rglob("*")
                      if p.is_file() and "__pycache__" not in p.parts}
    return {"status": out["status"], "bounded_checks": out["bounded_checks"],
            "source_inputs": len(out["input_hashes"]), "source_files_preserved": len(before),
            "mode": "existing --check; isolated temporary output; not source --record or independent peer review"}


def links():
    for p in HERE.glob("*.md"):
        assert all(line == line.rstrip() for line in p.read_text().splitlines()), p
        text = re.sub(r"```.*?```", "", p.read_text(), flags=re.S)
        text = re.sub(r"`[^`\n]*`", "", text)
        for target in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", text):
            if not target.startswith(("http:", "https:")):
                assert (p.parent/target.split("#")[0]).is_file(), (p, target)


def main():
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    if args.check:
        manifest(HERE, True)
    inputs = {n: sha(ROOT/n) for n in INPUTS}
    protected = {n: manifest(ROOT/n) for n in PROTECTED}
    checks = []
    def ck(name, condition, **data):
        assert bool(condition), name
        checks.append({"name": name, "passed": True, **data})

    # Read actual source AST, do not import/run its numerical code.
    tree = ast.parse((ROOT/INPUTS[2]).read_text())
    dual = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "dual7")
    assignments = {n.targets[0].id: n.value for n in dual.body if isinstance(n, ast.Assign)
                   and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name)}
    h, p = assignments["h"], assignments["p"]
    ck("actual canonical map is typed as affine logits then softmax",
       isinstance(h, ast.BinOp) and isinstance(h.op, ast.MatMult)
       and isinstance(h.left, ast.Name) and h.left.id == "X7"
       and isinstance(h.right, ast.Name) and h.right.id == "theta"
       and isinstance(p, ast.Call) and isinstance(p.func, ast.Name)
       and p.func.id == "softmax" and isinstance(p.args[0], ast.Name) and p.args[0].id == "h",
       scope="structural source binding, not derivation of a physical state map or dynamics")

    for exponent in (-3, -1, 0, 1, 3):
        even_weight = F(2)**exponent
        weights = [even_weight if i % 2 == 0 else F(1) for i in range(12)]
        mean, bit = parity(weights)
        ck(f"exact softmax parity reduction exponent{exponent}", mean == (even_weight-1)/(even_weight+1)
           and bit == int(exponent > 0))
        shifted = [F(3, 7)*w for w in weights]
        ck(f"common logit gauge cannot alter the bit exponent{exponent}", parity(shifted) == (mean, bit))
    general = [F((i+2)**2, i+1) for i in range(12)]
    mean, bit = parity(general)
    ck("general normalizer is positive and cancels in the sign", bit == int(sum(((-1)**i)*general[i] for i in range(12)) > 0)
       and -1 < mean < 1)

    tie_queries = [(F(1), F(0)), (F(-1), F(0))]
    tie_words, _ = all_affine_words(tie_queries)
    ck("atoms on zero faces are not discarded", tie_words == {(1, 0), (0, 1), (0, 0)},
       scope="three equal-mass seed atoms would yield three equiprobable words; open chambers alone miss one")
    ck("constant and identically zero query edge cases", all_affine_words([(F(0), F(1)), (F(0), F(0))])[0] == {(1, 0)})
    first = [(F(1), F(1)), (F(1), F(0)), (F(1), F(-1))]
    second = [(-a, -b) for a, b in first]
    w0, samples0 = all_affine_words(first)
    w1, _ = all_affine_words(second)
    ck("initial discrete-branch union needs its resource", len(w0) == len(w1) == 4
       and len(w0 | w1) == 6 <= 2*4)
    cantor = [sum((F(2*b, 3**(j+1)) for j, b in enumerate(digits)), F(0))
              for digits in product((0, 1), repeat=4)]
    ck("a singular-support subset does not create new words", {word(x, first) for x in cantor} <= w0,
       initial_points=len(cantor), scope="finite Cantor-type control; the arbitrary-measure theorem is analytical")
    ck("univariate full-cell enumeration includes all boundary samples", len(samples0) == 7
       and len(w0) <= 2*3+1)

    # Hypothesis countercontrols, not proposed FIN mechanisms.
    bits = 8
    seeds = [F(2*k+1, 2**(bits+1)) for k in range(2**bits)]
    adaptive, nonlinear = set(), set()
    same = True
    for x in seeds:
        low, high = F(0), F(1)
        adaptive_word = []
        for _ in range(bits):
            mid = (low+high)/2
            value = int(x > mid)
            adaptive_word.append(value)
            if value:
                low = mid
            else:
                high = mid
        z = x
        nonlinear_word = []
        for _ in range(bits):
            value = int(z >= F(1, 2))
            nonlinear_word.append(value)
            z = 2*z-value
            assert 0 <= z < 1
        adaptive.add(tuple(adaptive_word))
        nonlinear.add(tuple(nonlinear_word))
        same &= adaptive_word == nonlinear_word
    ck("adaptive affine queries cannot be treated as fixed queries", len(adaptive) == 2**bits > 2*bits+1,
       scope="exact binary-search countercontrol; not a computation of C for the general softmax family")
    ck("nonlinear definable step does not imply a uniform family of iterates", len(nonlinear) == 2**bits and same,
       scope="doubling/folding map countercontrol, not a FIN provider or an excluded nonlinear action")

    # Fixed-d source capacity is sublinear in the frozen native record length.
    n, d = 144, 7
    transitions = n*n
    observed_bits = transitions+1
    ck("record indexing and canonical source dimension bookkeeping", observed_bits == 20737
       and 12*(d+1) == 96 and d == 7)
    ck("discrete initial-resource normalization", F(n, transitions) == F(1, n),
       scope="coefficient of log12 after dividing by M=N^2; no finite-N error certified")
    ck("actual model-clock horizon retains sqrtN", 12*transitions == 248832 and 12**2 == n,
       scope="sqrt(N)Mt at t=1 only; no numerical mixing threshold")
    dependency = replay_bi()
    ck("read-only native information dependency replay", dependency["bounded_checks"] == 36
       and dependency["status"] == "CRITICAL_ONE_BIT_INFORMATION_BUDGET_WITH_INITIAL_RESOURCE_PAID")

    literature = json.loads((HERE/"LITERATURE_REGISTER.json").read_text())
    cache_states = []
    for source in literature["sources"]:
        cache = Path(source["optional_local_cache"])
        if cache.exists():
            assert cache.stat().st_size == source["pdf_bytes"] and sha(cache) == source["pdf_sha256"]
        cache_states.append("PRESENT_AND_HASH_MATCHED" if cache.exists() else "ABSENT_OPTIONAL_CACHE")
    named = [n for n in INPUTS if n.endswith(".md") and n != "AGENTS.md"]
    pattern = r"o[- ]minimal|VC[ -](dimension|density)|Vapnik|Sauer|shatter|continuous.{0,25}information|softmax"
    search = subprocess.run(["rg", "-n", "-i", "--max-columns", "300", "--", pattern, *named],
                            cwd=ROOT, capture_output=True, text=True, timeout=10)
    assert search.returncode in (0, 1), search.stderr
    assert inputs == {n: sha(ROOT/n) for n in INPUTS}
    assert protected == {n: manifest(ROOT/n) for n in PROTECTED}
    for pth in HERE.glob("*.py"):
        ast.parse(pth.read_text())
    for pth in HERE.glob("*.json"):
        json.loads(pth.read_text())
    result = {"date": "2026-10-21", "status": "CANONICAL_LINEAR_SOURCE_CAPACITY_BOUND_WITH_CONTINUOUS_SEEDS",
              "bounded_checks": len(checks), "checks": checks, "input_hashes": inputs,
              "protected_manifests": protected, "BI_dependency_replay": dependency,
              "published_math": "Wilkie Theorem5.2; ADHMS Theorem1.1 and duality, imported analytically not proved by fixtures",
              "word_bound": "R_N C(B+1)^d, d and formula fixed; C not computed",
              "source_measure": "arbitrary joint discrete/continuous law, no finite-entropy or definable-measure assumption",
              "native_record": "M_N=N^2 transitions, B_N=N^2+1 bits at tau_N+sqrt(N)jt",
              "source_error_floor": "liminf e_N >= b2^-1(c(t))>0 for R_N<=12^N",
              "duplicate_search": {"files": named, "pattern": pattern,
                                   "matching_lines": sorted(search.stdout.splitlines()),
                                   "scope": "bounded named-source search, not global novelty"},
              "imported_mathematical_theorems_mechanized": False,
              "full_nonlinear_FIN_source_excluded": False,
              "F01": "OPEN_SOURCE", "independent_review": False,
              "AGENTS_or_originals_modified": False}
    intext = "".join(f"{h}  {n}\n" for n, h in sorted(inputs.items()))
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result, indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(intext)
        files = sorted(pth for pth in HERE.rglob("*") if pth.is_file() and pth.name != "MANIFEST.sha256"
                       and "__pycache__" not in pth.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(pth)}  {pth.relative_to(HERE).as_posix()}\n" for pth in files))
    else:
        saved = json.loads((HERE/"RESULTS.json").read_text())
        assert result == saved, [key for key in result if result[key] != saved.get(key)]
        assert intext == (HERE/"INPUTS.sha256").read_text()
    links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "bounded_LS_checks": len(checks),
                      "BI_replayed_checks_not_new_LS_assertions": dependency["bounded_checks"],
                      "inputs": len(inputs), "protected_manifest_entries": sum(v["entries"] for v in protected.values()),
                      "optional_PDF_caches": cache_states, "C_or_usable_N0_computed": False,
                      "F01": "OPEN_SOURCE", "AGENTS_or_originals_modified": False}, indent=2))


if __name__ == "__main__":
    main()
