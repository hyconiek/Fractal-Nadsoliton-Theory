#!/usr/bin/env python3
"""Exact finite information controls, not a mechanical native-limit proof."""
import argparse
import ast
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import hashlib
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
    "fin_native_transfer_review_20261016/PROOF_COMPLETIONS_AND_SCOPE.md",
    "fin_native_transfer_review_20261016/INTAKE_20261016.md",
    "fin_native_transfer_campaign_20261015/NT-002/STRONG_LIFTED_TRANSFER_PROOF.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "fin_critical_metastability_campaign_20261014/MT-003_CONDITIONAL_GIBBS.md",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md",
    "fin_native_time_reversal_entropy_20261021/CONTRACT.md",
    "fin_native_time_reversal_entropy_20261021/EXPONENTIAL_PREPARATION_AND_NATIVE_WINDOW.md",
    "fin_native_time_reversal_entropy_20261021/RESULTS.json",
    "fin_intrinsic_entropy_source_review_20261019/PROOF_AUDIT_AND_SCOPE.md",
    "fin_intrinsic_entropy_source_review_20261019/SUBSAMPLING_ERROR_FLOOR.md",
    "fin_finite_history_source_research_20261020/NATIVE_LIMIT_ORDER.md",
    "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
]
PROTECTED = [
    "fin_native_transfer_review_20261016", "fin_critical_metastability_review_20261015",
    "fin_native_time_reversal_entropy_20261021", "fin_intrinsic_entropy_source_review_20261019",
    "fin_finite_history_source_research_20261020",
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


@dataclass(frozen=True)
class Box:
    lo: F
    hi: F

    def __post_init__(self):
        assert self.lo <= self.hi

    def __add__(self, b):
        return Box(self.lo+b.lo, self.hi+b.hi)

    def __neg__(self):
        return Box(-self.hi, -self.lo)

    def __sub__(self, b):
        return self+-b

    def __mul__(self, c):
        values = (self.lo*c, self.hi*c)
        return Box(min(values), max(values))

    __rmul__ = __mul__


def series_log(x):
    assert 1 <= x <= 2
    z = (x-1)/(x+1)
    m = 65
    partial = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(m)), F(0))
    tail = 2*z**(2*m+1)/(F(2*m+1)*(1-z*z))
    return Box(partial, partial+tail)


@lru_cache(None)
def logarithm(x):
    x = F(x)
    assert x > 0
    power = 0
    while x < 1:
        x *= 2
        power -= 1
    while x >= 2:
        x /= 2
        power += 1
    return series_log(x)+series_log(F(2))*power


@lru_cache(None)
def binary_entropy(p):
    p = F(p)
    assert 0 <= p <= 1
    if p in (0, 1):
        return Box(F(0), F(0))
    return logarithm(p)*(-p)+logarithm(1-p)*(-(1-p))


def entropy(p):
    out = Box(F(0), F(0))
    assert sum(p.values()) == 1 and min(p.values()) >= 0
    for mass in p.values():
        if mass:
            out = out+logarithm(mass)*(-mass)
    return out


def tv(p, q):
    return sum((abs(p.get(x, F(0))-q.get(x, F(0))) for x in set(p) | set(q)), F(0))/2


def paths(mu, kernel, transitions):
    out = {}
    for word in product(range(len(mu)), repeat=transitions+1):
        mass = mu[word[0]]
        for a, b in zip(word, word[1:]):
            mass *= kernel[a][b]
        if mass:
            out[word] = mass
    return out


def bit_records(p):
    out = defaultdict(F)
    for word, mass in p.items():
        out[tuple(int(k == 2) for k in word)] += mass
    return dict(out)


def de_bruijn(order):
    work = [0]*(2*order+1)
    sequence = []
    def build(t, period):
        if t > order:
            if order % period == 0:
                sequence.extend(work[1:period+1])
        else:
            work[t] = work[t-period]
            build(t+1, period)
            for value in range(work[t-period]+1, 2):
                work[t] = value
                build(t+1, t)
    build(1, 1)
    return sequence


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


def links():
    for p in HERE.glob("*.md"):
        assert all(line == line.rstrip() for line in p.read_text().splitlines()), p
        text = re.sub(r"```.*?```", "", p.read_text(), flags=re.S)
        text = re.sub(r"`[^`\n]*`", "", text)
        for target in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", text):
            if not target.startswith(("https:", "http:")):
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

    log2 = logarithm(F(2))
    ck("natural-log normalization", F(2, 3) < log2.lo < log2.hi < F(7, 10))
    ck("binary endpoints and symmetry", binary_entropy(F(0)) == binary_entropy(F(1)) == Box(F(0), F(0))
       and binary_entropy(F(1, 8)) == binary_entropy(F(7, 8)))
    for p in (F(1, 8), F(1, 2), F(7, 8)):
        hb = binary_entropy(p)
        ck(f"entropy analytic bounds fixture p{p}", hb.lo >= 2*p*(1-p)
           and hb.hi**2 <= 4*min(p, 1-p))
    ck("Gaussian threshold positivity inequality", F(103, 120)**2/6 < F(3, 8)**2
       and 2*(2*F(1, 8)*F(7, 8)) == F(7, 16))
    ck("local cubic stopping constants", F(1, 3)*F(3, 4) == F(1, 4)
       and F(1, 3)**2*F(3, 4)/2 < F(1, 4))

    # Actual g0 labelled HB: at model time log2, P_copy=(I+uniform)/2.
    kernel = [[F(9, 16), F(6, 16), F(1, 16)],
              [F(3, 16), F(10, 16), F(3, 16)],
              [F(1, 16), F(6, 16), F(9, 16)]]
    pi = [F(1, 4), F(1, 2), F(1, 4)]
    labels = list(product(range(12), repeat=2))
    count_rows_ok = True
    for x in labels:
        even_x = sum(a % 2 == 0 for a in x)
        row = [F(0)]*3
        for y in labels:
            mass = (F(x[0] == y[0], 2)+F(1, 24))*(F(x[1] == y[1], 2)+F(1, 24))
            row[sum(a % 2 == 0 for a in y)] += mass
        count_rows_ok &= row == kernel[even_x]
    ck("actual q12 N2 g0 collective-count skeleton", count_rows_ok,
       labelled_starts=len(labels), target_states_per_start=len(labels),
       scope="exact neutral finite control at model time log2; not an MT-domain or critical-lag approximation")
    ck("actual g0 stationary orbit masses and detailed balance", all(sum(row) == 1 for row in kernel)
       and all(pi[i]*kernel[i][j] == pi[j]*kernel[j][i] for i, j in product(range(3), repeat=2)))
    c = Box(F(0), F(0))
    for i in range(3):
        c = c+binary_entropy(kernel[i][2])*pi[i]
    ck("full-state conditional bit entropy is strictly positive", c.lo > 0)
    twice = [[sum(kernel[i][k]*kernel[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    c_twice = Box(F(0), F(0))
    for i in range(3):
        c_twice = c_twice+binary_entropy(twice[i][2])*pi[i]
    ck("conditional entropy increases with a longer native g0 lag", c_twice.lo > c.hi)
    memory = bit_records(paths(pi, kernel, 2))
    p_one = memory[(1, 0, 1)]/sum(memory[(1, 0, y)] for y in (0, 1))
    p_zero = memory[(0, 0, 1)]/sum(memory[(0, 0, y)] for y in (0, 1))
    ck("observed bit is genuinely non-Markov", p_one == F(19, 112)
       and p_zero == F(93, 656) and p_one != p_zero,
       previous_one=str(p_one), previous_zero=str(p_zero))

    epsilon = F(1, 100)
    mu = [pi[0]+epsilon, pi[1], pi[2]-epsilon]
    for transitions in (1, 2, 3):
        observed = bit_records(paths(pi, kernel, transitions))
        h_observed = entropy(observed)
        ck(f"full-state conditional lower bound for a finite bit block M{transitions}",
           h_observed.lo >= transitions*c.hi)
        fano = binary_entropy(F(1, 4))*(transitions+1)
        ck(f"binary Fano control with constant source M{transitions}",
           fano.lo >= h_observed.hi,
           scope="source constant-zero record, mismatch 1/4; no iid error assumption")
        actual = bit_records(paths(mu, kernel, transitions))
        h_actual = entropy(actual)
        error = tv(actual, observed)
        continuity = binary_entropy(epsilon)+log2*(epsilon*(transitions+1))
        ck(f"finite-alphabet entropy continuity M{transitions}", error <= epsilon
           and max(h_actual.hi-h_observed.lo, h_observed.hi-h_actual.lo) <= continuity.lo)
    pair = entropy(bit_records(paths(pi, kernel, 1)))-binary_entropy(F(1, 4))
    ck("pair entropy is not silently called the full-state bound", pair.lo > c.hi)

    # A finite deterministic periodic source can spend initial information on short fair blocks.
    order = 5
    sequence = de_bruijn(order)
    ck("one fixed de Bruijn cycle", len(sequence) == 2**order and set(sequence) == {0, 1})
    for bits in (1, 3, 5, 8):
        words = Counter(tuple(sequence[(phase+j) % len(sequence)] for j in range(bits))
                        for phase in range(len(sequence)))
        masses = {word: F(count, len(sequence)) for word, count in words.items()}
        coefficient = min(bits, order)
        ck(f"finite seed pays exact periodic block entropy B{bits}", len(masses) == 2**coefficient
           and set(masses.values()) == {F(1, 2**coefficient)},
           entropy_in_units_log2=coefficient, phase_entropy_in_units_log2=order,
           scope="finite deterministic permutation, entropy rate zero; not a natural FIN provider")
    ck("finite KS alone cannot bound finite-block entropy", order > 0 and len(sequence) == 32,
       explanation="one 32-state deterministic cycle has positive five-bit block entropy and zero entropy rate")
    ck("finite codebook capacity rejects a too-small distortion", binary_entropy(F(1, 10)).hi < (log2*F(4, 5)).lo,
       scope="generic 25 fair bits versus at most32 initial seeds: (25-5)/25 log2; no native-N performance claim")

    # Exact arithmetic behind the coarse-lag and source-bound optimizations.
    t = F(3, 7)
    for m in (5, 11, 101):
        k = max(1, (m*t.numerator)//t.denominator)
        lag = F(k, m)
        ck(f"positive coarse-lag rounding m{m}", 0 < lag <= t and t-lag < F(1, m))
    coefficient, budget = F(3, 4), F(2)
    center = coefficient/(2*budget)
    ck("deficit quadratic completed square", coefficient**2/(4*budget)-budget*center**2 == 0
       and 2*budget*center == coefficient,
       scope="fixed algebraic control, not an optimization scan or a measured source budget")
    ck("critical/model clock factor retained", F(12)*7*F(2, 3) == 56 and 12**2 == 144,
       scope="clock arithmetic only: sqrt(N)Mt, not a finite-N entropy certificate")
    geometry = json.loads((ROOT/INPUTS[7]).read_text())
    ck("actual inherited spectral bytes", sha(ROOT/geometry["spectral_path"]) == geometry["spectral_sha256"])
    tr = json.loads((ROOT/INPUTS[11]).read_text())
    ck("actual late-iid coupling dependency is the paid TR result",
       tr["status"] == "NATIVE_POST_PREPARATION_REVERSAL_KL_WINDOW_PAID"
       and F(tr["iid_decay_rate"]) > 0)

    named = [n for n in INPUTS if n.endswith(".md") and n != "AGENTS.md"]
    pattern = r"entropy.rate|sign.{0,35}entropy|binary.entropy|conditional.entropy|initial.{0,30}(budget|information)|decimation"
    search = subprocess.run(["rg", "-n", "-i", "--max-columns", "300", "--", pattern, *named],
                            cwd=ROOT, capture_output=True, text=True, timeout=10)
    assert search.returncode in (0, 1), search.stderr
    assert inputs == {n: sha(ROOT/n) for n in INPUTS}
    assert protected == {n: manifest(ROOT/n) for n in PROTECTED}
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    result = {"date": "2026-10-21", "status": "CRITICAL_ONE_BIT_INFORMATION_BUDGET_WITH_INITIAL_RESOURCE_PAID",
              "bounded_checks": len(checks), "checks": checks, "input_hashes": inputs,
              "protected_manifests": protected,
              "conditional_entropy_asymptotic": "c(t)=sqrt(2)rho(0)[integral b2(Phi(u))du]sqrt(t)+o(sqrt(t))",
              "rate_scope": "c(t) is a lower bound, not the exact observed bit entropy rate",
              "source_floor_bound_asymptotic": "e_bit(H)~C_bit^2/[4H log H], a necessary bound not an achievable optimum",
              "initial_resource": "I_N is required in finite-record source budgets; finite KS alone is insufficient",
              "finite_carrier_corollary": "At most12^N complete deterministic initial states cannot yield vanishing mean-bit error on N^2-transition actual-parent records, conditionally on the stated native premises",
              "clock": "critical time; actual model horizon sqrt(N)Mt is explicit",
              "joint_small_lag_native_CLT": "NOT_USED_OR_CLAIMED",
              "duplicate_search": {"files": named, "pattern": pattern,
                                   "matching_lines": sorted(search.stdout.splitlines()),
                                   "scope": "bounded named-predecessor search, not exhaustive novelty"},
              "F01": "OPEN_SOURCE", "independent_review": False,
              "new_microscopic_law": False, "AGENTS_or_originals_modified": False}
    input_text = "".join(f"{h}  {n}\n" for n, h in sorted(inputs.items()))
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result, indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(input_text)
        files = sorted(p for p in HERE.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        saved = json.loads((HERE/"RESULTS.json").read_text())
        assert result == saved, [key for key in result if result[key] != saved.get(key)]
        assert input_text == (HERE/"INPUTS.sha256").read_text()
    links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "bounded_checks": len(checks), "inputs": len(inputs),
                      "protected_manifest_entries": sum(v["entries"] for v in protected.values()),
                      "non_Markov_native_g0_bit_control": True, "initial_resource_required": True,
                      "Fano_decimation": "REUSED", "F01": "OPEN_SOURCE",
                      "AGENTS_or_originals_modified": False}, indent=2))


if __name__ == "__main__":
    main()
