#!/usr/bin/env python3
"""Small exact controls and own-packet integrity, not an asymptotic proof."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools as it
import json
import math
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = [
    "AGENTS.md",
    "research_tasks/00_EXECUTION_RULES.md",
    "FIN_RESEARCH_EXECUTION_PLAN.md",
    "fin_causal_source_contract_20261021/AUTHOR_REQUIREMENT.md",
    "fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md",
    "fin_hyp_history_entropy_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_hyp_history_entropy_review_20261020/PROOF_AUDIT_AND_DUAL.md",
    "fin_hyp_history_entropy_review_20261020/INTAKE_20261020.md",
    "fin_reference_transport_research_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_channel_composition_hypothesis_20261021/AUTHORIZATION_AND_MODEL.md",
    "fin_channel_composition_hypothesis_20261021/PROOF_OR_OPEN_LEMMA.md",
    "fin_spectator_composition_research_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_global_activity_defect_research_20261020/GLOBAL_APPROXIMATE_SCALING.md",
    "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_hyp_history_entropy_20261020/MANIFEST.sha256",
    "fin_hyp_history_entropy_review_20261020/MANIFEST.sha256",
    "fin_reference_transport_research_20261020/MANIFEST.sha256",
    "fin_channel_composition_hypothesis_20261021/MANIFEST.sha256",
    "fin_spectator_composition_research_20261020/MANIFEST.sha256",
    "fin_global_activity_defect_research_20261020/MANIFEST.sha256",
]
PROTECTED = [p for p in INPUTS if p.endswith("/MANIFEST.sha256")]
SEARCH_FILES = [p for p in INPUTS if p.endswith(".md") and p != "AGENTS.md"]
SEARCH_PATTERN = r"Hamming|suscept|interaction.order|composition|next.event|clock|causal"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest(path, own=False):
    folder = path.parent
    seen = set()
    for line in path.read_text().splitlines():
        h, name = line.split("  ", 1)
        p = Path(name)
        assert not p.is_absolute() and ".." not in p.parts and name not in seen
        assert digest(folder / name) == h, (str(folder), name)
        seen.add(name)
    actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*")
              if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
    complete = actual == seen | {"MANIFEST.sha256"}
    if own:
        assert complete, (seen, actual)
    return {"entries": len(seen), "manifest_sha256": digest(path),
            "coverage": "complete_noncache" if complete else "listed_members_only"}


def solve(a, b):
    """Independent exact dense elimination; at most 27 states."""
    n = len(b)
    rows = [[F(v) for v in a[i]] + [F(b[i])] for i in range(n)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        value = rows[j][j]
        rows[j] = [x/value for x in rows[j]]
        for i in range(n):
            if i != j and rows[i][j]:
                value = rows[i][j]
                rows[i] = [x-value*y for x, y in zip(rows[i], rows[j])]
    return [row[-1] for row in rows]


def carrier(q, n, self_weight=1):
    states = list(it.product(range(q), repeat=n))
    k = [[self_weight if x == y else int(sum(a != b for a, b in zip(x, y)) == 1)
          for y in states] for x in states]
    return states, k


def matvec(k, v):
    return [sum((a*b for a, b in zip(row, v)), F(0)) for row in k]


def pot(x, a):
    return sum((a[x[i]][x[j]] for i in range(len(x)) for j in range(i)), F(0))/len(x)


def field(x, i, a):
    return [sum((a[b][x[j]] for j in range(len(x)) if j != i), F(0))/len(x)
            for b in range(len(a))]


def log_rate_coeff(x, i, b, a):
    h = field(x, i, a)
    old = x[i]
    return ((h[b]-h[old])/2, h[b]+h[old]-2*sum(h)/len(h))


def count_key(x):
    return tuple(sorted(x))


def direct_jets(states, k, u, v, x):
    ix = states.index(x)
    legal = [j for j, y in enumerate(states) if j != ix and k[ix][j]]
    slopes = {j: u[ix]+u[j]-v[ix] for j in legal}
    avg = sum(slopes.values())/len(legal)
    mass, derivative = {}, {}
    for j in legal:
        key = count_key(states[j])
        mass[key] = mass.get(key, F(0))+F(1, len(legal))
        derivative[key] = derivative.get(key, F(0))+(slopes[j]-avg)/len(legal)
    return mass, derivative


def numerical_scaling(states, k, a, g):
    """One q3/N3 centred-difference fixture; NOT an outward enclosure."""
    weights = [math.exp(g*float(pot(x, a))) for x in states]
    pi = [p/sum(weights) for p in weights]
    d = sum(k[0])
    z = [math.sqrt(p/d) for p in pi]
    residual = math.inf
    for iteration in range(2000):
        kz = [sum(c*t for c, t in zip(row, z)) for row in k]
        residual = max(abs(s*t/p-1) for s, t, p in zip(z, kz, pi))
        if residual < 4e-15:
            break
        z = [math.sqrt(s*p/t) for s, p, t in zip(z, pi, kz)]
    assert residual < 4e-15, residual
    x = (0, 0, 1)
    i = states.index(x)
    row = {}
    for j, y in enumerate(states):
        if j != i and k[i][j]:
            key = count_key(y)
            row[key] = row.get(key, 0.0)+z[j]
    total = sum(row.values())
    row = {key: value/total for key, value in row.items()}
    ja, jb, jc = [row[key] for key in ((0, 0, 0), (0, 0, 2), (0, 1, 2))]
    delta = float(a[0][0]-a[0][1])
    fvalue = math.log(jc/(2*math.sqrt(ja*jb)))+g*delta/3
    dvalue = jc*jc-4*math.exp(-2*g*delta/3)*ja*jb
    return {"F": fvalue, "D": dvalue, "residual": residual,
            "iterations": iteration, "kind": "ordinary_binary64_diagnostic"}


def diagnostics():
    checks = []

    def ck(name, value, kind="exact_fraction", **data):
        assert bool(value), (name, data)
        checks.append({"name": name, "passed": True, "kind": kind, **data})

    # Direct tensor-product function checks, not a precomputed graph spectrum.
    q, n = 3, 3
    states, k = carrier(q, n)
    d = sum(k[0])
    centred = [F(-1), F(-2), F(3)]
    for r in range(4):
        vector = [math.prod(centred[x[i]] for i in range(r)) for x in states]
        ck(f"direct Hamming interaction-order action r{r}",
           matvec(k, vector) == [(d-r*q)*t for t in vector])
    ck("positive all-order margin denominators q3 N3",
       min(2*d-r*q for r in range(n+1)) == n*(q-2)+2 > 0)

    a = [[F(q*int(i == j)-1) for j in range(q)] for i in range(q)]
    v = [pot(x, a) for x in states]
    ck("partition derivative vanishes for the row-centred pair potential", sum(v) == 0)
    matrix = [[F(k[i][j])+d*int(i == j) for j in range(len(states))]
              for i in range(len(states))]
    u = solve(matrix, [d*t for t in v])
    ck("independent 27-state rational susceptibility solve", u == [F(7, 8)*t for t in v])
    ck("full differentiated margin equations", matvec(matrix, u) == [d*t for t in v])
    ck("normalization derivative from susceptibility", sum(u) == 0)

    x = (0, 0, 1)
    keys = ((0, 0, 0), (0, 0, 2), (0, 1, 2))
    j, jp = direct_jets(states, k, u, v, x)
    ja, jb, jc = [j[t] for t in keys]
    ja1, jb1, jc1 = [jp[t] for t in keys]
    delta = a[0][0]-a[0][1]
    fp = jc1/jc-(ja1/ja+jb1/jb)/2+delta/3
    dp = 2*jc*jc1-4*((-2*delta/3)*ja*jb+ja1*jb+ja*jb1)
    ck("count multiplicities from direct legal labelled neighbours",
       (ja, jb, jc) == (F(1, 6), F(1, 6), F(1, 3)))
    ck("direct normalized-branch derivatives retain denominator term", sum(jp.values()) == 0)
    ck("direct logarithmic contrast derivative q3", fp == F(-3, 4))
    ck("direct bounded contrast derivative q3", dp == F(-1, 6))

    # All s coefficients cancel: these are exact affine identities, not sample s values.
    for q in (3, 4, 12):
        for name in ("potts", "cycle"):
            a_test = [[F(q*int(i == j)-1) if name == "potts" else
                       F(2*int(i == j)-int((i-j) % q in (1, q-1)))
                       for j in range(q)] for i in range(q)]
            h_a = log_rate_coeff(x, 2, 0, a_test)
            h_b = log_rate_coeff(x, 2, 2, a_test)
            h_c = log_rate_coeff(x, 0, 2, a_test)
            delta_test = a_test[0][0]-a_test[0][1]
            ck(f"all-parameter exact CC contrast q{q} {name}",
               h_c[0]-(h_a[0]+h_b[0])/2+delta_test/3 == 0
               and h_c[1]-(h_a[1]+h_b[1])/2 == 0)
            vv = [pot(t, a_test) for t in keys]
            ck(f"destination-potential identity q{q} {name}",
               vv[2]-(vv[0]+vv[1])/2 == -2*delta_test/3)

    q, n = 12, 3
    degree = 1+n*(q-1)
    c = F(degree, 2*(degree-q))
    fp_delta = (1-2*c)/3
    dp_delta = F(4, (3*(q-1))**2)*2*fp_delta
    tv_delta = -dp_delta/6
    cert = json.loads((HERE/"EXACT_CERTIFICATE.json").read_text())
    ck("q12 exact coefficient certificate", degree == cert["degree"]
       and c == F(cert["pair_scaling_susceptibility"])
       and fp_delta == F(cert["F_derivative_per_Delta"])
       and dp_delta == F(cert["D_derivative_per_Delta"])
       and tv_delta == F(cert["TV_liminf_per_Delta"]))
    ck("q12 direct neutral count factors", [F(cert["neutral_event_probabilities"][t])
       for t in ("A", "B", "C")] == [F(1, 33), F(1, 33), F(2, 33)])
    ck("self multiplicity mutation changes susceptibility",
       F(35, 2*(35-12)) != c and F(8, 2*(8-3)) == F(4, 5))
    ck("omitting the C multiplicity fails already at zero", jc*jc-ja*jb != 0)
    rates = [F(2, 17), F(3, 17), F(4, 17)]
    scale = F(3, 7)
    ck("common state clock cancels exactly from branch law",
       [r/sum(rates) for r in rates] == [scale*r/sum(scale*t for t in rates) for r in rates])
    ck("bounded readout payoff range for g nonnegative",
       max(F(1), F(0), F(-2))-min(F(1), F(0), F(-2)) == 3)

    # Exactly one non-A7 nonlinear fixture, with a centred finite difference.
    a = [[F(3*int(i == j)-1) for j in range(3)] for i in range(3)]
    minus = numerical_scaling(states, k, a, -1e-4)
    plus = numerical_scaling(states, k, a, 1e-4)
    f_est = (plus["F"]-minus["F"])/2e-4
    d_est = (plus["D"]-minus["D"])/2e-4
    ck("nonlinear q3 control agrees with exact first derivatives",
       abs(f_est+0.75) < 1e-7 and abs(d_est+1/6) < 1e-7,
       "ordinary_binary64", F_central=f_est, D_central=d_est)
    ck("same positive-g control has the adverse sign", plus["F"] < 0 and plus["D"] < 0,
       "ordinary_binary64")

    # Bind and use versioned A7 only as a nominal finite arithmetic control.
    module = ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
    spec = importlib.util.spec_from_file_location("ch_safe_seed", module)
    seed = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(seed)
    _, _, lam, _, a7 = seed.build_rank7()
    delta7 = float(a7[0, 0]-a7[0, 1])
    ck("nominal versioned A7 parity lower bound",
       delta7 > 0 and delta7 >= float(lam[6])/6-1e-12,
       "ordinary_binary64", Delta=delta7, parity_lower=float(lam[6])/6)
    nominal = {"Delta": delta7, "F_prime": -2*delta7/11,
               "D_prime": -16*delta7/11979,
               "TV_liminf_lower_coefficient": 8*delta7/35937,
               "scope": "nominal non-outward diagnostics, not certificate endpoints"}

    return {"status": "CONDITIONAL_NEAR_ZERO_INCOMPATIBILITY",
            "date": "2026-10-21", "checks": checks, "check_count": len(checks),
            "exact_check_count": sum(c["kind"] == "exact_fraction" for c in checks),
            "nonlinear_control": {"minus": minus, "plus": plus},
            "nominal_A7": nominal}


def lint():
    ast.parse((HERE/"verify_research.py").read_text())
    for path in HERE.glob("*.json"):
        json.loads(path.read_text())
    for path in HERE.glob("*.md"):
        txt = path.read_text()
        assert all(line == line.rstrip() for line in txt.splitlines()), path
        prose = re.sub(r"```.*?```", "", txt, flags=re.S)
        prose = re.sub(r"(?m)^(?: {4}|\t).*(?:\n|$)", "", prose)
        prose = re.sub(r"`[^`\n]*`", "", prose)
        for target in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", prose):
            if not target.startswith(("https:", "http:")):
                assert (path.parent/target.split("#")[0]).is_file(), (path, target)


def main():
    ap = argparse.ArgumentParser()
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--record", action="store_true")
    group.add_argument("--check", action="store_true")
    group.add_argument("--diagnose", action="store_true")
    args = ap.parse_args()
    before = {p: digest(ROOT/p) for p in INPUTS}
    protected = {p: manifest(ROOT/p) for p in PROTECTED}
    if args.check:
        manifest(HERE/"MANIFEST.sha256", True)
    result = diagnostics()
    result["inputs"] = before
    result["protected_manifests"] = protected
    result["protected_entries"] = sum(v["entries"] for v in protected.values())
    search = subprocess.run(["rg", "-n", "-i", "-e", SEARCH_PATTERN, *SEARCH_FILES],
                            cwd=ROOT, text=True, capture_output=True, timeout=15)
    assert search.returncode in (0, 1), search.stderr
    search_data = {"scope": SEARCH_FILES, "pattern": SEARCH_PATTERN,
                   "status": "bounded_discovery_aid_not_exhaustive",
                   "matches": search.stdout.splitlines()}
    assert before == {p: digest(ROOT/p) for p in INPUTS}, "source mutation during checks"
    for p in PROTECTED:
        assert protected[p] == manifest(ROOT/p)
    result_text = json.dumps(result, indent=2, sort_keys=True)+"\n"
    input_text = "".join(f"{h}  {p}\n" for p, h in before.items())
    search_text = json.dumps(search_data, indent=2, sort_keys=True)+"\n"
    if args.record:
        (HERE/"RESULTS.json").write_text(result_text)
        (HERE/"INPUTS.sha256").write_text(input_text)
        (HERE/"DUPLICATE_SEARCH.json").write_text(search_text)
        lint()
        files = sorted(p for p in HERE.rglob("*") if p.is_file()
                       and p.name != "MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(
            f"{digest(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    elif args.check:
        assert (HERE/"RESULTS.json").read_text() == result_text
        assert (HERE/"INPUTS.sha256").read_text() == input_text
        assert (HERE/"DUPLICATE_SEARCH.json").read_text() == search_text
        lint()
        manifest(HERE/"MANIFEST.sha256", True)
    print(json.dumps({"status": "PASS_BOUNDED_CHECKS", "mode":
          "record_own_metadata" if args.record else "read_only",
          "checks": result["check_count"], "exact_checks": result["exact_check_count"],
          "inputs": len(INPUTS), "protected_entries": result["protected_entries"],
          "nominal_A7": result["nominal_A7"]}, sort_keys=True))


if __name__ == "__main__":
    main()
