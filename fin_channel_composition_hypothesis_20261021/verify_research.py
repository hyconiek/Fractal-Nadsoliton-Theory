#!/usr/bin/env python3
"""Bounded exact controls; not an asymptotic or physical proof."""
import argparse
import ast
from collections import Counter
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "FIN_CAUSAL_SOURCE_TRANSFER_NEXT_CHAIN_20261020.md",
    "FIN_POST_DEL_SOURCE_DECISION_20261021.md",
    "fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md",
    "fin_kinetic_composition_research_20261016/LITERATURE_COMPARISON.md",
    "fin_spectator_composition_research_20261020/REUSE_AND_PROGRESS.md",
    "fin_rate_coupling_review_20261016/ACTIVITY_DECISION.md",
    "fin_rate_coupling_derivation_20261016/CERTIFICATE.json",
    "fin_sqrt_critical_process_research_20261020/CONTRACT.md",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
    "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py",
    "FIN_ST2532_ST2621_Configuration_Generator_and_Dynamic_RG_Report_EN.md",
    "FIN_POST_AFTER_CONTINUATION_FULL_HANDOFF_20260926/01_REPORTS/KINETIC_NATURALITY_43.md",
    "FIN_POST_AFTER_CONTINUATION_FULL_HANDOFF_20260926/01_REPORTS/FISHER_KINETIC_METRIC_44.md",
]
PROTECTED = [
    "fin_hyp_history_entropy_20261020",
    "fin_rate_coupling_review_20261016",
    "fin_sqrt_critical_process_research_20261020",
    "fin_critical_metastability_review_20261015",
    "fin_kinetic_composition_research_20261016",
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder, coverage=False):
    seen = set()
    for line in (folder / "MANIFEST.sha256").read_text().splitlines():
        h, name = line.split("  ", 1)
        path = Path(name)
        assert name not in seen and not path.is_absolute() and ".." not in path.parts
        assert sha(folder / name) == h, (folder, name)
        seen.add(name)
    if coverage:
        actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
        assert actual == seen | {"MANIFEST.sha256"}, folder
    return {"entries": len(seen), "sha256": sha(folder / "MANIFEST.sha256"),
            "coverage": "complete_noncache" if coverage else "listed_members_only"}


def rank(rows):
    a = [list(map(F, row)) for row in rows]
    pivot = 0
    for j in range(len(a[0])):
        k = next((i for i in range(pivot, len(a)) if a[i][j]), None)
        if k is None:
            continue
        a[pivot], a[k] = a[k], a[pivot]
        d = a[pivot][j]
        a[pivot] = [x / d for x in a[pivot]]
        for i in range(len(a)):
            if i != pivot:
                d = a[i][j]
                a[i] = [x - d*y for x, y in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


def coefficients(q, a, b, s):
    return [s*((i == a)+(i == b)-F(2, q))
            + F((i == b)-(i == a), 2) for i in range(q)]


def dot(x, y):
    return sum((a*b for a, b in zip(x, y)), F(0))


def generator(off):
    out = [row[:] for row in off]
    for i, row in enumerate(out):
        assert row[i] == 0
        row[i] = -sum(row)
    return out


def balance(q, weights):
    return all(weights[i]*q[i][j] == weights[j]*q[j][i]
               for i, j in it.product(range(len(weights)), repeat=2))


def equal_pairs(x):
    return sum(x[i] == x[j] for i in range(len(x)) for j in range(i))


def conditional_weights(x, slot, q):
    bg = Counter(x[:slot] + x[slot+1:])
    return [F(2**bg[j]) for j in range(q)]


def check_links():
    for p in HERE.glob("*.md"):
        assert all(line == line.rstrip() for line in p.read_text().splitlines()), p
        txt = re.sub(r"```.*?```", "", p.read_text(), flags=re.S)
        txt = re.sub(r"`[^`\n]*`", "", txt)
        for link in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", txt):
            if not link.startswith(("https:", "http:")):
                assert (p.parent / link.split("#")[0]).is_file(), (p, link)


def main():
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    if args.check:
        manifest(HERE, True)
    inputs = {n: sha(ROOT / n) for n in INPUTS}
    protected = {n: manifest(ROOT / n) for n in PROTECTED}
    checks = []

    def ck(name, value, **data):
        assert bool(value), name
        checks.append({"name": name, "passed": True, **data})

    # Completeness checks the linear constraints, rather than only substituting examples.
    for q in (2, 3, 4, 12):
        rows = [[F(1)]*q, [F((i == 1)-(i == 0)) for i in range(q)]]
        rows += [[F((i == j)-(i == 2)) for i in range(q)] for j in range(3, q)]
        ck(f"coefficient constraint rank q{q}", rank(rows) == (2 if q == 2 else q-1))
        for s in (F(-1, 2), F(0), F(1, 2)):
            c = coefficients(q, 0, 1, s)
            reverse = coefficients(q, 1, 0, s)
            u = [F(i*i-2*i+1, 7) for i in range(q)]
            v = [F((-1)**i*(i+1), 5) for i in range(q)]
            ck(f"response identities q{q} s{s}",
               sum(c) == 0 and c[1]-c[0] == 1
               and [a-b for a, b in zip(c, reverse)] == [F((i == 1)-(i == 0)) for i in range(q)]
               and dot(c, [a+b for a, b in zip(u, v)]) == dot(c, u)+dot(c, v))
        if q == 2:
            ck("q2 symmetric direction identically zero", coefficients(q, 0, 1, F(17)) == coefficients(q, 0, 1, F(0)))
        else:
            ck(f"nonzero remaining direction q{q}", coefficients(q, 0, 1, F(1)) != coefficients(q, 0, 1, F(0)))

    q = 3
    c = coefficients(q, 0, 1, F(2, 7))
    covariance = []
    for p in it.permutations(range(q)):
        cp = coefficients(q, p[0], p[1], F(2, 7))
        covariance.append(all(cp[p[i]] == c[i] for i in range(q)))
    ck("all q3 label permutations", all(covariance))

    # h=2 log(w), product(w)=1 gives exact rational rates at the three fixed fixtures.
    w = [F(1, 2), F(1), F(2)]
    weights = [a*a for a in w]
    off0 = [[F(0) if i == j else w[j]/w[i]/q for j in range(q)] for i in range(q)]
    offp = [[F(0) if i == j else w[j]**2/q for j in range(q)] for i in range(q)]
    offm = [[F(0) if i == j else F(1)/(q*w[i]**2) for j in range(q)] for i in range(q)]
    matrices = [generator(x) for x in (off0, offp, offm)]
    ck("three different positive same-Gibbs generators", all(balance(x, weights) for x in matrices)
       and len({tuple(map(tuple, x)) for x in matrices}) == 3)
    ck("remaining activity is not just common clock", offp[0][1]/off0[0][1] != offp[0][2]/off0[0][2])
    hb = [[F(0) if i == j else weights[j]/sum(weights) for j in range(q)] for i in range(q)]
    ck("s half matches isolated HB modulo one clock", len({offp[i][j]/hb[i][j]
       for i, j in it.product(range(q), repeat=2) if i != j}) == 1)
    z, z2 = sum(weights), sum(a*a for a in weights)
    ck("HB doubling composition ratio strictly below one", z*z/(q*z2) == F(7, 13))
    constant_z, constant_z2 = sum([F(4)]*q), sum([F(16)]*q)
    ck("uniform field equality retained", constant_z**2/(q*constant_z2) == 1)
    ck("distance-three thinning violates multiplicativity", F(1, 5) != F(1, 2)**2)
    ck("endpoint-only dependence was not silently imposed", coefficients(3, 0, 1, F(1, 2))[2] != 0)

    # A separate integer-potential model, not A7. Check exact whole Gibbs edge balance.
    for q in (3, 12):
        states = list(it.product(range(q), repeat=3)) if q == 3 else [(0, 0, 1)]
        edges = 0
        okay = True
        for x in states:
            wx = F(2**equal_pairs(x))
            for slot in range(3):
                old = x[slot]
                p = conditional_weights(x, slot, q)
                for new in range(q):
                    if new == old:
                        continue
                    y = x[:slot] + (new,) + x[slot+1:]
                    wy = F(2**equal_pairs(y))
                    rev = conditional_weights(y, slot, q)
                    okay &= wx*p[new] == wy*rev[old]
                    okay &= wx*p[new]/sum(p) == wy*rev[old]/sum(rev)
                    edges += 1
        ck(f"full-field Gibbs balance exact control q{q}", okay, directed_edges=edges,
           scope="all q3 N3 states; one q12 N3 state and reverse edges; different matrix")
        x = (0, 0, 1)
        hom, mixed = conditional_weights(x, 2, q), conditional_weights(x, 0, q)
        ck(f"strict clock split q{q}", sum(hom) == q+3 and sum(mixed) == q+2
           and hom[0]*hom[1] == mixed[0]*mixed[1])
        roots = [F(2) if j == 0 else F(1) for j in range(q)]
        shifted_roots = [roots[(j-1) % q] for j in range(q)]
        ck(f"square identity for clock split q{q}", sum(hom)-sum(mixed)
           == sum((a-b)**2 for a, b in zip(roots, shifted_roots))/2 == 1)
        a = (sum(hom)-hom[1])/sum(hom)
        b = 2*(sum(mixed)-mixed[0])/sum(mixed)
        p_hb = a/(a+b)
        p_comp = (sum(hom)-hom[1])/((sum(hom)-hom[1])+2*(sum(mixed)-mixed[0]))
        ck(f"bounded next-count-event gap q{q}", p_comp > p_hb,
           HB=str(p_hb), composition=str(p_comp), difference=str(p_comp-p_hb))
        if q == 12:
            ck("published rational control endpoints", (p_hb, p_comp, p_comp-p_hb)
               == (F(49, 139), F(7, 19), F(42, 2641)))

    old, new = (0, 0, 0), (1, 0, 0)
    u, v = conditional_weights(old, 0, 3), conditional_weights(new, 0, 3)
    for name, normalizer in (("common-clock composition", F(3)), ("HB", sum(u))):
        forward, reverse = 3*u[1]/normalizer, v[0]/normalizer
        ck(f"count multinomial and repeated-label factor: {name}",
           8*forward == (3*2)*reverse and 8*forward != 2*reverse)

    # Independent direct differentiation of fluxes at uniformity.
    for q in (3, 12):
        w = [F(i, q*q) for i in range(q)]
        w = [a-sum(w)/q for a in w]
        h = [F((-1)**i*(i+1), 17) for i in range(q)]
        h = [a-sum(h)/q for a in h]
        linear = []
        for s in (F(-1, 2), F(0), F(1, 2)):
            row = []
            for j in range(q):
                inflow = sum(w[i]/q + dot(coefficients(q, i, j, s), h)/q**2
                             for i in range(q) if i != j)
                outflow = (q-1)*w[j]/q + sum(dot(coefficients(q, j, b, s), h)/q**2
                                                          for b in range(q) if b != j)
                row.append(inflow-outflow)
            linear.append(row)
        ck(f"direct full-flux linear response q{q}", all(row == [-w[j]+h[j]/q for j in range(q)] for row in linear))

    q, m, ex = 12, F(1, 5), F(2)
    sigma = [(-1)**j for j in range(q)]
    p = [(1+m*a)/q for a in sigma]
    for two_s in (-1, 0, 1):
        drift = F(0)
        for a, b in it.product(range(q), repeat=2):
            if a == b:
                continue
            exponent = F(sigma[b]-sigma[a], 2)+F(two_s, 2)*(sigma[a]+sigma[b])
            assert exponent.denominator == 1
            rate = ex**int(exponent)/q
            drift += p[a]*rate*(sigma[b]-sigma[a])
        ck(f"exact parity-line drift 2s={two_s}", drift == ((1-m)*ex-(1+m)/ex)/2)
    second = sum(F(1, q*q)*(sigma[b]-sigma[a])**2
                 for a, b in it.product(range(q), repeat=2) if a != b)
    ck("neutral parity second jump moment", second == 2)
    wcap = [F(1, 4), F(1), F(4)]
    ck("unit-attempt countercontrol is an exit not a single marked jump", sum(wcap[j]/wcap[0]/3 for j in (1, 2)) == F(20, 3) > 1)

    # Inherited actual inputs, not a synthetic reconstruction called A7.
    geometry = json.loads((ROOT/INPUTS[10]).read_text())
    ck("actual MT spectral input hash", sha(ROOT/geometry["spectral_path"]) == geometry["spectral_sha256"])
    certificate = json.loads((ROOT/INPUTS[8]).read_text())
    lo, hi = map(F, certificate["laplacian_active_spectra"]["6"])
    ck("actual nonzero parity premise bound by saved certificate", 0 < lo <= hi)
    ck("actual canonical-gap premise is not flat matrix control", F(certificate["gap_A00_minus_A03"][0]) > 0)

    files = [n for n in INPUTS if n.endswith(".md") and n != "AGENTS.md"]
    pattern = r"multiplicat|channel|square.root|kinetic.{0,30}(ratio|metric)|activity"
    search = subprocess.run(["rg", "-n", "-i", "--max-columns", "300", "--", pattern, *files],
                            cwd=ROOT, capture_output=True, text=True, timeout=10)
    assert search.returncode in (0, 1), search.stderr
    assert inputs == {n: sha(ROOT/n) for n in INPUTS}
    assert protected == {n: manifest(ROOT/n) for n in PROTECTED}
    for pth in HERE.glob("*.py"):
        ast.parse(pth.read_text())
    for pth in HERE.glob("*.json"):
        json.loads(pth.read_text())
    result = {
        "date": "2026-10-21", "model_id": "FIN-HYP-CHANNEL-COMPOSITION-002",
        "status": "CONDITIONAL_ACTIVITY_RESTRICTION_WITH_FREE_PARAMETER_AND_N3_HB_OBSTRUCTION",
        "entry_head": "a3397ccad56c4d74c7c676a066876ac271783b62",
        "bounded_checks": len(checks), "checks": checks,
        "input_hashes": inputs, "protected_manifests": protected,
        "duplicate_search": {"files": files, "pattern": pattern,
                             "matching_lines": sorted(search.stdout.splitlines()),
                             "scope": "bounded named-file comparison, not global novelty"},
        "actual_A7_N3_result": "analytical from supplied symmetry and nonzero parity; not the rational-control number",
        "F01": "OPEN_SOURCE", "independent_peer_review": False,
        "physical_model_adoption": False, "source_or_AGENTS_rewritten": False,
    }
    intext = "".join(f"{h}  {n}\n" for n, h in sorted(inputs.items()))
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result, indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(intext)
        paths = sorted(p for p in HERE.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256"
                       and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in paths))
    else:
        saved = json.loads((HERE/"RESULTS.json").read_text())
        if result != saved:
            changed = [key for key in result if result[key] != saved.get(key)]
            raise AssertionError("replay differs in: " + ", ".join(changed))
        assert intext == (HERE/"INPUTS.sha256").read_text()
    check_links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "bounded_checks": len(checks),
                      "inputs": len(inputs), "protected_manifest_entries": sum(v["entries"] for v in protected.values()),
                      "free_parameter": "s", "N3_native_HB_transfer": "REJECTED_FOR_ALL_s",
                      "F01": result["F01"], "source_or_AGENTS_rewritten": False}, indent=2))


if __name__ == "__main__":
    main()
