#!/usr/bin/env python3
"""Exact bounded certificate for one supplied finite operator, not all locality."""
import argparse
import ast
from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations
import hashlib
import json
from math import factorial
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PORTS = (0, 1, 2, 3, 6)
CORE = [((0, 1, 2, 3), ((0, 2), (1, 3))),
        ((0, 1, 3, 6), ((0, 6), (1, 3))),
        ((0, 2, 3, 6), ((0, 3), (2, 6)))]
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "FIN_RESEARCH_STRATEGY.md", "FIN_POST_PHY_TASK_CATALOG.md",
    "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py",
    "fin_toe_bridge_research_20261015/DERIVED_RESULTS.md",
    "fin_toe_bridge_research_20261015/verify_bridges.py",
    "fin_toe_bridge_research_20261015/check_locality_certificate.py",
    "fin_toe_bridge_research_20261015/verification.json",
    "fin_toe_bridge_research_20261015/LITERATURE_DATABASE.md",
    "fin_toe_bridge_research_20261015/PRIORITY_TASKS.json",
    "fin_channel_composition_hypothesis_20261021/HANDOFF.md",
]
PROTECTED = ["fin_toe_bridge_research_20261015", "fin_channel_composition_hypothesis_20261021"]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


@dataclass(frozen=True)
class Box:
    lo: F
    hi: F

    def __post_init__(self):
        assert self.lo <= self.hi

    @staticmethod
    def of(x):
        return x if isinstance(x, Box) else Box(F(x), F(x))

    def __add__(self, other):
        b = self.of(other)
        return Box(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return Box(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.of(other)

    def __mul__(self, other):
        b = self.of(other)
        terms = [a*c for a in (self.lo, self.hi) for c in (b.lo, b.hi)]
        return Box(min(terms), max(terms))

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.of(other)
        assert b.lo > 0
        return self * Box(1/b.hi, 1/b.lo)


def interval(values):
    return Box(*map(F, values))


def public(b, digits=12):
    scale = 10**digits
    low = b.lo.numerator*scale//b.lo.denominator
    high = -((-b.hi.numerator*scale)//b.hi.denominator)
    def dec(n):
        sign = "-" if n < 0 else ""
        n = abs(n)
        return sign+str(n//scale)+"."+str(n % scale).zfill(digits)
    out = [dec(low), dec(high)]
    assert F(out[0]) <= b.lo <= b.hi <= F(out[1])
    return out


def root5(n):
    scale = 10**24
    target = n*scale**5
    low, high = 0, (n+1)*scale
    while high-low > 1:
        mid = (low+high)//2
        if mid**5 <= target:
            low = mid
        else:
            high = mid
    assert low**5 <= target < (low+1)**5
    top = low if low**5 == target else low+1
    return Box(F(low, scale), F(top, scale))


def weight(d):
    x = F(743, 4000)*d+F(13, 80)
    assert 0 <= x <= F(13, 10) and x*x < 2
    upper = sum((F((-1)**j, factorial(2*j))*x**(2*j) for j in range(11)), F(0))
    lower = upper-x**22/factorial(22)
    assert lower > 0
    r = root5(d**4)
    return Box(lower, upper)/(1+d*r)


def distance(a, b):
    assert a != b and 0 <= a < 12 and 0 <= b < 12
    return min((a-b) % 12, (b-a) % 12)


def matching(pairs):
    return tuple(sorted(tuple(sorted(p)) for p in pairs))


def possibilities(q):
    a, b, c, d = q
    return [matching(((a, b), (c, d))), matching(((a, c), (b, d))), matching(((a, d), (b, c)))]


def product(w, pairs):
    (a, b), (c, d) = pairs
    return w[distance(a, b)]*w[distance(c, d)]


def orders():
    return [(0,)+p for p in permutations(PORTS[1:]) if p[0] < p[-1]]


def crossing(order, quartet):
    a, b, c, d = [p for p in order if p in quartet]
    return matching(((a, c), (b, d)))


def circular(order, p, q):
    assert len(p) == len(q) == 2 and len(set(p+q)) == 4
    actual = [v for v in order if v in set(p+q)]
    target = p+q[::-1]
    return any(actual == target[k:]+target[:k] for k in range(4))


def determinant(w, p, q):
    return w[distance(p[0], q[0])]*w[distance(p[1], q[1])] - w[distance(p[0], q[1])]*w[distance(p[1], q[0])]


def build_certificate():
    raw = {d: weight(d) for d in range(1, 7)}
    public_w = {str(d): public(v) for d, v in raw.items()}
    w = {int(d): interval(v) for d, v in public_w.items()}
    constraints = []
    for quartet, forced in CORE:
        forced = matching(forced)
        gaps = []
        for alternative in possibilities(quartet):
            if alternative != forced:
                gap = product(w, alternative)-product(w, forced)
                assert gap.lo > F(1, 20000)
                gaps.append({"alternative": alternative, "gap": public(gap)})
        constraints.append({"quartet": quartet, "forced_crossing": forced, "gaps": gaps})
    witnesses = []
    for order in orders():
        for quartet, forced in CORE:
            forced = matching(forced)
            if crossing(order, quartet) == forced:
                continue
            a, b, c, d = [x for x in order if x in quartet]
            if forced == matching(((a, b), (c, d))):
                p, q = [b, c], [a, d]
            else:
                assert forced == matching(((a, d), (b, c)))
                p, q = [a, b], [d, c]
            assert circular(order, p, q)
            det = determinant(w, p, q)
            assert det.hi < -F(1, 20000)
            witnesses.append({"order": order, "quartet": quartet, "P": p, "Q": q, "minor": public(det)})
            break
        else:
            raise AssertionError("compatible order unexpectedly exists")
    epsilon = F(1, 100000)
    delta = F(1, 20000)-2*epsilon-2*epsilon**2
    assert delta > 0 and max(v.hi for v in w.values()) < F(1, 2)
    return {
        "date": "2026-10-21", "operator": "M=-L12, full unprojected supplied finite operator",
        "parameters": {"q": 12, "omega": "743/4000", "phi": "13/80", "eta": "9/5"},
        "method": "Exact fifth-root brackets at denominator 10^24; alternating cosine Taylor degree20/22; outward rational decimal rounding; no interval backend dependency",
        "kernel_weights_raw": {str(d): [str(b.lo), str(b.hi)] for d, b in raw.items()},
        "kernel_weights": public_w,
        "ports": PORTS, "constraints": constraints, "order_witnesses": witnesses,
        "matching_gap_floor": "1/20000", "perturbation_epsilon": str(epsilon),
        "perturbed_negative_minor_margin": str(delta),
        "scope": "Any port order, any finite positive one-disk network; not arbitrary locality or A7",
    }


def check_certificate(c):
    assert c["parameters"] == {"q": 12, "omega": "743/4000", "phi": "13/80", "eta": "9/5"}
    w = {int(k): interval(v) for k, v in c["kernel_weights"].items()}
    for d in range(1, 7):
        rebuilt = weight(d)
        assert w[d].lo <= rebuilt.lo <= rebuilt.hi <= w[d].hi
    assert list(c["ports"]) == list(PORTS)
    assert len(c["constraints"]) == 3
    for saved, (quartet, forced) in zip(c["constraints"], CORE):
        assert tuple(saved["quartet"]) == quartet and matching(saved["forced_crossing"]) == matching(forced)
        assert len(saved["gaps"]) == 2
        used = set()
        for item in saved["gaps"]:
            alternative = matching(item["alternative"])
            assert alternative in possibilities(quartet) and alternative != matching(forced)
            used.add(alternative)
            value = product(w, alternative)-product(w, matching(forced))
            shown = interval(item["gap"])
            assert shown.lo <= value.lo <= value.hi <= shown.hi
            assert value.lo > F(1, 20000)
        assert len(used) == 2
    seen = set()
    for item in c["order_witnesses"]:
        order = tuple(item["order"])
        assert order in orders() and order not in seen
        seen.add(order)
        assert circular(order, item["P"], item["Q"])
        actual = determinant(w, item["P"], item["Q"])
        shown = interval(item["minor"])
        assert shown.lo <= actual.lo <= actual.hi <= shown.hi < -F(1, 20000)
    assert seen == set(orders())
    epsilon = F(c["perturbation_epsilon"])
    assert epsilon == F(1, 100000)
    assert F(c["matching_gap_floor"]) == F(1, 20000)
    assert F(c["perturbed_negative_minor_margin"]) == F(1, 20000)-2*epsilon-2*epsilon**2 > 0


def manifest(folder, coverage=False):
    seen = set()
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        h, name = line.split("  ", 1)
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name) == h, (folder, name)
        seen.add(name)
    if coverage:
        present = {p.relative_to(folder).as_posix() for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
        assert present == seen | {"MANIFEST.sha256"}
    return {"entries": len(seen), "sha256": sha(folder/"MANIFEST.sha256"),
            "coverage": "complete_noncache" if coverage else "listed_members_only"}


def links():
    for p in HERE.glob("*.md"):
        assert all(line == line.rstrip() for line in p.read_text().splitlines()), p
        text = re.sub(r"```.*?```", "", p.read_text(), flags=re.S)
        text = re.sub(r"`[^`\n]*`", "", text)
        for link in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", text):
            if not link.startswith(("http:", "https:")):
                assert (p.parent/link.split("#")[0]).is_file(), (p, link)


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

    c = build_certificate()
    check_certificate(c)
    ck("complete certificate including ordered minor signs", len(c["order_witnesses"]) == 12)
    # A second combinatorial formulation ignores saved witnesses and intersects constraints.
    for subset in combinations(range(3), 2):
        possible = [o for o in orders() if all(crossing(o, CORE[i][0]) == matching(CORE[i][1]) for i in subset)]
        ck(f"each pair of core constraints is compatible {subset}", bool(possible), compatible_orders=len(possible))
    all_three = [o for o in orders() if all(crossing(o, q) == matching(m) for q, m in CORE)]
    ck("independent crossing-set intersection is empty", not all_three)
    ck("root perfect-power endpoint", root5(1) == Box(F(1), F(1)))
    for d in range(1, 7):
        b = root5(d**4)
        ck(f"exact algebraic root endpoints d{d}", b.lo**5 <= d**4 <= b.hi**5)
    old = json.loads((ROOT/"fin_toe_bridge_research_20261015/verification.json").read_text())["certificate"]
    for d in range(1, 7):
        older = interval(old["kernel_weight_inputs"][str(d)])
        current = interval(c["kernel_weights"][str(d)])
        ck(f"old exact-decimal weight enclosure agrees d{d}", current.lo <= older.lo <= older.hi <= current.hi)
    # Bind displayed safe_seed decimal parameters without executing numpy or confusing float with real arithmetic.
    p = ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
    text = p.read_text()
    values = {node.targets[0].id: ast.get_source_segment(text, node.value)
              for node in ast.parse(text).body if isinstance(node, ast.Assign)
              and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)}
    ck("supplied decimal kernel parameters bound to versioned source",
       F(values["OMEGA"]) == F(743, 4000) and F(values["PHI"]) == F(13, 80)
       and F(values["ETA"]) == F(9, 5) and int(values["Q"]) == 12)

    bad = deepcopy(c)
    bad["order_witnesses"][0]["minor"] = ["0", "1"]
    try:
        check_certificate(bad)
    except AssertionError:
        ck("altered positive minor certificate rejected", True)
    else:
        raise AssertionError("bad positive witness accepted")
    bad = deepcopy(c)
    bad["constraints"][0]["forced_crossing"] = ((0, 1), (2, 3))
    try:
        check_certificate(bad)
    except AssertionError:
        ck("wrong forced crossing rejected", True)
    else:
        raise AssertionError("bad matching accepted")
    witness = c["order_witnesses"][0]
    ck("reversing columns is not silently kept circular",
       not circular(witness["order"], witness["P"], witness["Q"][::-1]))

    # Actual passive-star response, independently obtained by eliminating its central node.
    conductances = dict(zip(PORTS, map(F, (1, 2, 3, 4, 5))))
    total = sum(conductances.values())
    def star(a, b):
        return conductances[a]*conductances[b]/total
    star_checks = []
    for order in orders():
        for q in combinations(PORTS, 4):
            a, b, d, e = [v for v in order if v in q]
            star_checks.append(star(a, e)*star(b, d)-star(a, d)*star(b, e) == 0)
            star_checks.append(star(a, b)*star(d, e)-star(a, d)*star(b, e) == 0)
    ck("positive planar-star control has no forbidden quartet signs", all(star_checks), minors_checked=len(star_checks))
    epsilon = F(1, 100000)
    ck("uniform perturbation bound leaves strict margin",
       F(c["perturbed_negative_minor_margin"]) > F(2, 100000))
    ck("diagonal corrections cannot alter disjoint-port minors",
       all(not set(w["P"]) & set(w["Q"]) for w in c["order_witnesses"]))

    named = [n for n in INPUTS if n.endswith(".md") and n != "AGENTS.md"]
    pattern = r"REF-004|other.{0,20}order|natural.{0,15}order|larger minor|planar"
    search = subprocess.run(["rg", "-n", "-i", "--max-columns", "300", "--", pattern, *named],
                            cwd=ROOT, capture_output=True, text=True, timeout=10)
    assert search.returncode in (0, 1), search.stderr
    assert inputs == {n: sha(ROOT/n) for n in INPUTS}
    assert protected == {n: manifest(ROOT/n) for n in PROTECTED}
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    result = {"date": "2026-10-21", "status": "ORDER_FREE_PASSIVE_DISK_OBSTRUCTION_WITH_RATIONAL_ROBUSTNESS",
              "bounded_checks": len(checks), "checks": checks, "input_hashes": inputs,
              "protected_manifests": protected, "order_witnesses": 12, "quartet_constraints": 3,
              "ports": 5, "perturbation_epsilon": str(epsilon),
              "duplicate_search": {"files": named, "pattern": pattern,
                                   "matching_lines": sorted(search.stdout.splitlines()),
                                   "scope": "bounded named-artifact search; not exhaustive novelty"},
              "external_theorem": "Alman-Lian-Tran 1309.2697v1 Theorem 2.2.6(a), necessity only",
              "actual_input_arithmetic": "integer/Fraction only; no libm or interval backend for kernel signs",
              "F01": "OPEN_SOURCE", "unrestricted_locality_no_go": False,
              "independent_review": False, "originals_or_AGENTS_rewritten": False}
    input_text = "".join(f"{h}  {n}\n" for n, h in sorted(inputs.items()))
    if args.record:
        (HERE/"CERTIFICATE.json").write_text(json.dumps(c, indent=2)+"\n")
        (HERE/"RESULTS.json").write_text(json.dumps(result, indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(input_text)
        paths = sorted(p for p in HERE.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in paths))
    else:
        assert json.loads(json.dumps(c)) == json.loads((HERE/"CERTIFICATE.json").read_text())
        assert result == json.loads((HERE/"RESULTS.json").read_text())
        assert input_text == (HERE/"INPUTS.sha256").read_text()
    links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    literature = json.loads((HERE/"LITERATURE_REGISTER.json").read_text())["sources"][0]
    cache = Path(literature["optional_local_cache"])
    if cache.exists():
        assert cache.stat().st_size == literature["pdf_bytes"] and sha(cache) == literature["pdf_sha256"]
    print(json.dumps({"status": result["status"], "bounded_checks": len(checks),
                      "inputs": len(inputs), "protected_manifest_entries": sum(x["entries"] for x in protected.values()),
                      "certified_orders_of_five_ports": 12, "required_crossing_constraints": 3,
                      "all_twelve_port_orders": "EXCLUDED_BY_SUBORDER_CONTRADICTION",
                      "optional_pdf_cache": "PRESENT_AND_HASH_MATCHED" if cache.exists() else "ABSENT_OPTIONAL_CACHE",
                      "F01": "OPEN_SOURCE", "AGENTS_or_originals_rewritten": False}, indent=2))


if __name__ == "__main__":
    main()
