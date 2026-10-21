#!/usr/bin/env python3
"""Exact finite entropy controls, not a mechanized MT/asymptotic proof."""
import argparse
import ast
from collections import defaultdict
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
    "fin_critical_metastability_review_20261015/INTAKE_20261015.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_critical_metastability_campaign_20261014/MT-001_COERCIVITY_CONNECTIVITY.md",
    "fin_critical_metastability_campaign_20261014/MT-002_IID_TRAPPING.md",
    "fin_critical_metastability_campaign_20261014/MT-004_REFLECTED_MIXING_WINDOW.md",
    "fin_intrinsic_observation_review_20261015/REVERSIBILITY_AND_RECORD_SCOPE.md",
    "fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md",
    "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py",
    "fin_intrinsic_entropy_source_review_20261019/PROOF_AUDIT_AND_SCOPE.md",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
]
PROTECTED = ["fin_critical_metastability_review_20261015", "fin_intrinsic_observation_review_20261015"]


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


@dataclass(frozen=True)
class Box:
    lo: F
    hi: F

    def __post_init__(self):
        assert self.lo <= self.hi

    def __add__(self, b):
        return Box(self.lo+b.lo, self.hi+b.hi)

    def __mul__(self, c):
        p = (self.lo*c, self.hi*c)
        return Box(min(p), max(p))

    __rmul__ = __mul__


@lru_cache(None)
def log_prime_interval(p):
    assert p >= 2
    z = F(p-1, p+1)
    m = 300
    partial = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(m)), F(0))
    tail = 2*z**(2*m+1)/(F(2*m+1)*(1-z*z))
    return Box(partial, partial+tail)


def factor(n):
    assert isinstance(n, int) and n > 0
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def log_terms(x):
    x = F(x)
    assert x > 0
    out = {p: F(k) for p, k in factor(x.numerator).items()}
    for p, k in factor(x.denominator).items():
        out[p] = out.get(p, F(0))-k
    return {p: k for p, k in out.items() if k}


def combine(a, b, scale=F(1)):
    out = dict(a)
    for p, k in b.items():
        out[p] = out.get(p, F(0))+scale*k
    return {p: k for p, k in out.items() if k}


def evaluate(terms):
    out = Box(F(0), F(0))
    for p, k in terms.items():
        out = out+log_prime_interval(p)*k
    return out


def kl(p, q):
    out = {}
    for x, mass in p.items():
        if mass:
            if q.get(x, F(0)) <= 0:
                raise ValueError("infinite relative entropy: support mismatch")
            out = combine(out, log_terms(mass/q[x]), mass)
    return out


def tv(p, q):
    return sum((abs(p.get(x, F(0))-q.get(x, F(0))) for x in set(p) | set(q)), F(0))/2


def paths(mu, kernel, steps):
    out = {}
    for word in product(range(len(mu)), repeat=steps+1):
        mass = mu[word[0]]
        for a, b in zip(word, word[1:]):
            mass *= kernel[a][b]
        if mass:
            out[word] = mass
    return out


def reverse(p):
    return {word[::-1]: v for word, v in p.items()}


def aggregate(p, function):
    out = defaultdict(F)
    for x, v in p.items():
        out[function(x)] += v
    return dict(out)


def vec(v):
    return dict(enumerate(v))


def reversible(kernel, pi):
    return all(pi[i]*kernel[i][j] == pi[j]*kernel[j][i]
               for i, j in product(range(len(pi)), repeat=2))


def integral(states, holds, values, time):
    out = F(0)
    remaining = time
    for x, duration in zip(states, holds):
        used = min(remaining, duration)
        out += values[x]*used
        remaining -= used
        if remaining <= 0:
            break
    return out


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
    input_hashes = {n: sha(ROOT/n) for n in INPUTS}
    protected = {n: manifest(ROOT/n) for n in PROTECTED}
    checks = []
    def ck(name, condition, **data):
        assert bool(condition), name
        checks.append({"name": name, "passed": True, **data})

    # A tiny reset CTMC has exp(t(Pi-I))=exp(-t)I+(1-exp(-t))Pi.
    pi = [F(1, 3)]*3
    mu = [F(1, 7), F(2, 7), F(4, 7)]
    projection = [pi[:] for _ in pi]
    kernel = [[F(i == j, 2)+pi[j]/2 for j in range(3)] for i in range(3)]
    ck("exact reset projection and log2-time embedding",
       all(sum(projection[i][k]*projection[k][j] for k in range(3)) == projection[i][j]
           for i, j in product(range(3), repeat=2)) and reversible(kernel, pi))
    log2 = log_prime_interval(2)
    ck("rational log-series sanity enclosure", F(2, 3) < log2.lo < log2.hi < F(7, 10))
    for steps in (1, 2, 3):
        p = paths(mu, kernel, steps)
        rp = reverse(p)
        end = aggregate(p, lambda x: x[-1])
        mu_t = [end[i] for i in range(3)]
        ck(f"endpoint likelihood for all sampled words M{steps}", sum(p.values()) == 1
           and all(p[x]/rp[x] == (mu[x[0]]/pi[x[0]])/(mu[x[-1]]/pi[x[-1]]) for x in p))
        a = kl(p, rp)
        endpoints = aggregate(p, lambda x: (x[0], x[-1]))
        ck(f"endpoint entropy sufficiency M{steps}", a == kl(endpoints, reverse(endpoints)))
        e = kl(p, reverse(paths(mu_t, kernel, steps)))
        drop = combine(kl(vec(mu), vec(pi)), kl(vec(mu_t), vec(pi)), F(-1))
        ck(f"correct final-preparation entropy identity M{steps}", e == drop
           and a == combine(e, kl(vec(mu_t), vec(mu))))
        ck(f"production nonnegative and at most same-preparation KL M{steps}",
           evaluate(e).lo >= 0 and evaluate(combine(a, e, F(-1))).lo >= 0)
        ratios = [mu[i]/pi[i] for i in range(3)]
        bound = {p0: tv(p, rp)*c for p0, c in log_terms(max(ratios)/min(ratios)).items()}
        ck(f"TV times actual likelihood-span bound M{steps}", evaluate(combine(bound, a, F(-1))).lo >= 0)
        observed = aggregate(p, lambda x: tuple(int(i == 2) for i in x))
        obs_kl = kl(observed, reverse(observed))
        ck(f"equivariant binary record entropy contraction M{steps}", evaluate(combine(a, obs_kl, F(-1))).lo >= 0)
    stationary = paths(pi, kernel, 3)
    ck("stationary reversible path is exactly symmetric", stationary == reverse(stationary)
       and kl(stationary, reverse(stationary)) == {})

    # Separate unmarked CTMC words: real holding factors, not only a skeleton identity.
    pi_c = [F(1, 2), F(1, 3), F(1, 6)]
    conductance = {(0, 1): F(1, 7), (0, 2): F(1, 11), (1, 2): F(1, 13)}
    rates = [[F(0) if i == j else conductance[tuple(sorted((i, j)))]/pi_c[i]
              for j in range(3)] for i in range(3)]
    exits = list(map(sum, rates))
    for word in ((0, 1, 2), (2, 0, 1, 2), (1, 0)):
        length = len(word)
        holds = [F(2*(i+1), length*(length+1)) for i in range(length)]
        forward = mu[word[0]]
        backward = mu[word[-1]]
        for i, j in zip(word, word[1:]):
            assert i != j
            forward *= rates[i][j]
            backward *= rates[j][i]
        ck(f"CTMC jump and holding cancellation {word}", sum(holds) == 1
           and forward/backward == (mu[word[0]]/pi_c[word[0]])/(mu[word[-1]]/pi_c[word[-1]])
           and sum(exits[x]*t for x, t in zip(word, holds))
           == sum(exits[x]*t for x, t in zip(word[::-1], holds[::-1])))
    word, holds = (0, 1, 2), [F(1, 6), F(2, 6), F(3, 6)]
    ck("incorrect unreversed holding intervals detected", sum(exits[x]*t for x, t in zip(word, holds))
       != sum(exits[x]*t for x, t in zip(word[::-1], holds)))

    for n in (2, 5, 11):
        p = {0: 1-F(1, n)-F(1, n*2**n), 1: F(1, n), 2: F(1, n*2**n)}
        rp = {0: p[0], 1: p[2], 2: p[1]}
        ck(f"positive-support TV-to-KL countercontrol n{n}", min(p.values()) > 0
           and tv(p, rp) == (1-F(1, 2**n))/n
           and kl(p, rp) == {2: 1-F(1, 2**n)})
        ck(f"sharp reversal entropy bound n{n}", kl(p, rp) == {2: n*tv(p, rp)})
    try:
        kl({0: F(127, 128), 1: F(1, 128)}, {0: F(1)})
    except ValueError:
        ck("small TV to reflected-support law can have infinite KL", True)
    else:
        raise AssertionError("support mismatch lost")
    cycle = [[F(1, 2) if j in (i, (i+1) % 3) else F(0) for j in range(3)] for i in range(3)]
    cp = paths(pi, cycle, 1)
    ck("stationarity alone does not supply reversibility", not reversible(cycle, pi)
       and aggregate(cp, lambda x: x[-1]) == vec(pi))
    try:
        kl(cp, reverse(cp))
    except ValueError:
        ck("nonreversible endpoint-zero shortcut rejected", True)
    else:
        raise AssertionError("directed cycle unexpectedly reversible")

    # Actual fixed-ball geometry, using the true saved MT constants.
    geometry = json.loads((ROOT/INPUTS[4]).read_text())
    r, q = F(1, 10**6), 12
    h, eps, delta, b = map(F, (geometry["h"], geometry["epsilon"], geometry["delta"], geometry["b"]))
    kappa = 2*r*r/q
    ck("actual source domain and fixed-ball margin", F(1, q)-r > F(1, 24)
       and 4*r < eps and r < delta and 12*r*r < h and b == h/2)
    ck("exponential iid rate is positive and explicitly bound", kappa == F(1, 6*10**12) and 0 < kappa < b)
    ck("source spectral inputs remain actual bound bytes", sha(ROOT/geometry["spectral_path"]) == geometry["spectral_sha256"])
    n = 12*10**6
    counts = [n//q]*q
    counts[0] += 1
    counts[1] -= 1
    v = [F(x, n)-F(1, q) for x in counts]
    norm = sum(x*x for x in v)
    decrement = 2*(v[1]-v[0])/n+F(2, n*n)
    ck("large integer count fixture in the fixed ball, no enumeration", 0 < norm <= r*r
       and decrement == -F(2, n*n) and norm+decrement == 0,
       N=n, scope="one explicit geometric type; no MT mixing error at this N certified")
    ck("Bernoulli tilted-variance algebra", all(p*(1-p) == F(1, 4)-(p-F(1, 2))**2 <= F(1, 4)
       for p in (F(1, 7), F(1, 2), F(5, 6))))
    ck("long horizon decay direction", b-b/2 > 0 and kappa > 0,
       scope="rational sign check only; the all-a polynomial-exponential limit is analytical")

    # Exact actual q12/N2 g0 rate normalization, not a tiny critical-domain claim.
    states = list(product(range(12), repeat=2))
    off = {}
    totals = defaultdict(F)
    for x in states:
        for slot in range(2):
            for label in range(12):
                if label != x[slot]:
                    y = x[:slot]+(label,)+x[slot+1:]
                    off[x, y] = F(1, 12)
                    totals[x] += F(1, 12)
    ck("native g0 unit-attempt and reversible normalization", len(off) == 3168
       and all(off[x, y] == off[y, x] for x, y in off)
       and max(totals.values()) == F(11, 6) <= 2,
       scope="g0 algebra only, not the critical gain or MT domain")

    # Exact exchangeable-orbit TV and endpoint entropy isometry.
    labelled = list(product(range(3), repeat=2))
    groups = defaultdict(list)
    for x in labelled:
        groups[tuple(x.count(i) for i in range(3))].append(x)
    names = sorted(groups)
    cmu = {x: F(i+1, 21) for i, x in enumerate(names)}
    cnu = {x: F(6-i, 21) for i, x in enumerate(names)}
    lmu = {x: cmu[c]/len(group) for c, group in groups.items() for x in group}
    lnu = {x: cnu[c]/len(group) for c, group in groups.items() for x in group}
    ck("exact orbit-TV isometry", tv(lmu, lnu) == tv(cmu, cnu))
    lpi = {x: F(1, 9) for x in labelled}
    joint = {(x, y): lmu[x]*(F(x == y, 2)+lpi[y]/2) for x, y in product(labelled, repeat=2)}
    coarse = aggregate(joint, lambda xy: tuple(tuple(x.count(i) for i in range(3)) for x in xy))
    ck("count endpoints retain the full reversal KL", kl(joint, reverse(joint)) == kl(coarse, reverse(coarse)),
       scope="generic reversible exchangeable reset control, not native single-copy dynamics")

    holds = [F(1, 4), F(1, 2), F(1, 4)]
    word, values = [0, 1, 2], [F(-1), F(0), F(2)]
    total = integral(word, holds, values, F(1))
    times = [F(0), F(1, 8), F(1, 4), F(1, 2), F(3, 4), F(7, 8), F(1)]
    ck("integrated-record reversal preserves the origin", all(
        integral(word[::-1], holds[::-1], values, t) == total-integral(word, holds, values, 1-t) for t in times))
    ck("wrong unrepaired record reversal rejected", total != 0)

    named = [x for x in INPUTS if x.endswith(".md") and x != "AGENTS.md"]
    pattern = r"entropy production|reversal.{0,50}(KL|entropy)|entropy.{0,40}(unpromot|transfer)|path.TV|iid initial"
    search = subprocess.run(["rg", "-n", "-i", "--max-columns", "300", "--", pattern, *named],
                            cwd=ROOT, capture_output=True, text=True, timeout=10)
    assert search.returncode in (0, 1), search.stderr
    assert input_hashes == {n: sha(ROOT/n) for n in INPUTS}
    assert protected == {n: manifest(ROOT/n) for n in PROTECTED}
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    result = {"date": "2026-10-21", "status": "NATIVE_POST_PREPARATION_REVERSAL_KL_WINDOW_PAID",
              "bounded_checks": len(checks), "checks": checks, "input_hashes": input_hashes,
              "protected_manifests": protected, "fixed_ball_radius": str(r), "iid_decay_rate": str(kappa),
              "native_bound": "A_N(T)<=2 N[tau_N+2log(q)+G M_A] epsilon_N(T); E_N(T)<=A_N(T)",
              "preparation": "tau_N=(N+1)^(K+2), inherited existential gap exponent",
              "horizon": "uniform 0<=T<=exp(aN), every fixed 0<a<b; no numerical N0",
              "actual_to_reflected_KL": "INFINITE; not used",
              "finite_control_scope": "exact small diagnostics; inherited MT theorem and new native limit are analytical",
              "duplicate_search": {"files": named, "pattern": pattern, "matching_lines": sorted(search.stdout.splitlines()),
                                   "scope": "bounded named-predecessor search; not exhaustive novelty"},
              "F01": "OPEN_SOURCE", "independent_review": False,
              "new_microscopic_law": False, "originals_or_AGENTS_modified": False}
    input_text = "".join(f"{h}  {n}\n" for n, h in sorted(input_hashes.items()))
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
    literature = json.loads((HERE/"LITERATURE_REGISTER.json").read_text())["sources"][0]
    cache = Path(literature["optional_local_cache"])
    if cache.exists():
        assert cache.stat().st_size == literature["pdf_bytes"] and sha(cache) == literature["pdf_sha256"]
    print(json.dumps({"status": result["status"], "bounded_checks": len(checks), "inputs": len(input_hashes),
                      "protected_manifest_entries": sum(x["entries"] for x in protected.values()),
                      "iid_decay_rate": str(kappa), "same_and_final_preparations_distinguished": True,
                      "actual_to_reflected_KL": "INFINITE_NOT_USED", "F01": "OPEN_SOURCE",
                      "optional_pdf_cache": "PRESENT_AND_HASH_MATCHED" if cache.exists() else "ABSENT_OPTIONAL_CACHE",
                      "AGENTS_or_originals_modified": False}, indent=2))


if __name__ == "__main__":
    main()
