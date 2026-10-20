#!/usr/bin/env python3
"""Native fixed-population AD row stress; exact outward arithmetic only."""
import argparse
import ast
import hashlib
import json
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BITS = 80
SCALE = 1 << BITS


def ceildiv(a, b):
    return -((-a)//b)


class I:
    __slots__ = ("lo", "hi")

    def __init__(self, x=0):
        if isinstance(x, I):
            self.lo, self.hi = x.lo, x.hi
        else:
            f = F(x)
            self.lo = (f.numerator*SCALE)//f.denominator
            self.hi = ceildiv(f.numerator*SCALE, f.denominator)

    @classmethod
    def raw(cls, lo, hi):
        assert lo <= hi
        x = object.__new__(cls)
        x.lo, x.hi = lo, hi
        return x

    @staticmethod
    def cast(x):
        return x if isinstance(x, I) else I(x)

    def __add__(self, other):
        if isinstance(other, J):
            return NotImplemented
        y = I.cast(other)
        return I.raw(self.lo+y.lo, self.hi+y.hi)

    __radd__ = __add__

    def __neg__(self):
        return I.raw(-self.hi, -self.lo)

    def __sub__(self, other):
        if isinstance(other, J):
            return NotImplemented
        return self+-I.cast(other)

    def __rsub__(self, other):
        return I.cast(other)+-self

    def __mul__(self, other):
        if isinstance(other, J):
            return NotImplemented
        y = I.cast(other)
        if self.lo == self.hi == 0 or y.lo == y.hi == 0:
            return I(0)
        vals = (self.lo*y.lo, self.lo*y.hi, self.hi*y.lo, self.hi*y.hi)
        return I.raw(min(vals)//SCALE, ceildiv(max(vals), SCALE))

    __rmul__ = __mul__

    def inv(self):
        assert self.lo > 0 or self.hi < 0, "division interval crosses zero"
        return I.raw((SCALE*SCALE)//self.hi, ceildiv(SCALE*SCALE, self.lo))

    def __truediv__(self, other):
        if isinstance(other, J):
            return NotImplemented
        return self*I.cast(other).inv()

    def __rtruediv__(self, other):
        return I.cast(other)*self.inv()

    def sqrt(self):
        assert self.lo >= 0
        low = isqrt(self.lo*SCALE)
        upper = isqrt(self.hi*SCALE)
        if upper*upper < self.hi*SCALE:
            upper += 1
        return I.raw(low, upper)

    def exp(self):
        lower = exp_endpoint(self.lo)
        upper = exp_endpoint(self.hi)
        return I.raw(lower.lo, upper.hi)

    def intersect(self, other):
        y = I.cast(other)
        return I.raw(max(self.lo, y.lo), min(self.hi, y.hi))

    def contains(self, x):
        y = F(x)*SCALE
        return self.lo <= y <= self.hi

    def width(self):
        return F(self.hi-self.lo, SCALE)

    def box(self, digits=14):
        d = 10**digits
        lower = (self.lo*d)//SCALE
        upper = ceildiv(self.hi*d, SCALE)
        def fmt(n):
            sign = "-" if n < 0 else ""
            n = abs(n)
            return f"{sign}{n//d}.{n % d:0{digits}d}"
        return {"lower": fmt(lower), "upper": fmt(upper), "outward_decimal_digits": digits, "exact_dyadic": {"lo": self.lo, "hi": self.hi, "bits": BITS}}


@lru_cache(maxsize=50000)
def exp_endpoint(endpoint):
    if endpoint < 0:
        return exp_endpoint(-endpoint).inv()
    y = F(endpoint, SCALE)
    shifts = 0
    while y > F(1, 8):
        y /= 2
        shifts += 1
    terms = 20
    lower = sum((y**k/F(factorial(k)) for k in range(terms+1)), F(0))
    tail = y**(terms+1)/F(factorial(terms+1))/(1-y/F(terms+2))
    lo, hi = I(lower), I(lower+tail)
    result = I.raw(lo.lo, hi.hi)
    for _ in range(shifts):
        result = result*result
    return result


class J:
    """Outward intervals for value, first derivative and half second derivative."""
    __slots__ = ("v", "d", "s")

    def __init__(self, value=0, derivative=0, second=0):
        if isinstance(value, J):
            self.v, self.d, self.s = value.v, value.d, value.s
        else:
            self.v, self.d, self.s = I(value), I(derivative), I(second)

    @staticmethod
    def cast(x):
        return x if isinstance(x, J) else J(x)

    def __add__(self, other):
        y = J.cast(other)
        return J(self.v+y.v, self.d+y.d, self.s+y.s)

    __radd__ = __add__

    def __neg__(self):
        return J(-self.v, -self.d, -self.s)

    def __sub__(self, other):
        return self+-J.cast(other)

    def __rsub__(self, other):
        return J.cast(other)+-self

    def __mul__(self, other):
        y = J.cast(other)
        return J(self.v*y.v, self.v*y.d+self.d*y.v, self.v*y.s+self.d*y.d+self.s*y.v)

    __rmul__ = __mul__

    def inv(self):
        inv = self.v.inv()
        return J(inv, -self.d*inv*inv, self.d*self.d*inv*inv*inv-self.s*inv*inv)

    def __truediv__(self, other):
        return self*J.cast(other).inv()

    def __rtruediv__(self, other):
        return J.cast(other)*self.inv()

    def sqrt(self):
        v = self.v.sqrt()
        return J(v, self.d/(2*v), self.s/(2*v)-self.d*self.d/(8*v*v*v))

    def exp(self):
        v = self.v.exp()
        return J(v, v*self.d, v*(self.s+self.d*self.d/2))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def native_context():
    path = ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"
    geom = json.loads((ROOT/"fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json").read_text())
    assert sha(path) == geom["spectral_sha256"]
    raw = json.loads(path.read_text())["exact"]["laplacian_intervals"]
    seed_name = "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
    admitted = json.loads((ROOT/"fin_critical_metastability_review_20261015/verification.json").read_text())
    assert sha(ROOT/seed_name) == admitted["input_hashes"][seed_name]
    boxes = {k: I.raw(I(F(raw[k][0])).lo, I(F(raw[k][1])).hi) for k in (3, 4, 5, 6)}
    root3 = I(3).sqrt()
    cos3 = [I(1), I(0), I(-1), I(0)]
    cos4 = [I(1), I(F(-1, 2)), I(F(-1, 2))]
    cos5 = [I(1), -root3/2, I(F(1, 2)), I(0), I(F(-1, 2)), root3/2, I(-1), root3/2, I(F(-1, 2)), I(0), I(F(1, 2)), -root3/2]
    assert all((cos5[(j+1) % 12]+cos5[(j-1) % 12]+root3*cos5[j]).contains(0) for j in range(12))
    row = [(boxes[3]*cos3[j % 4]+boxes[4]*cos4[j % 3]+boxes[5]*cos5[j])/6+boxes[6]*((-1)**j)/12 for j in range(12)]
    matrix = [[row[(j-i) % 12] for j in range(12)] for i in range(12)]
    g = 12/boxes[6]
    assert sum(row, I(0)).contains(0)
    assert (sum((x*((-1)**j) for j, x in enumerate(row)), I(0))-boxes[6]).contains(0)
    assert (g*boxes[6]).contains(12)
    assert (12*row[0]-(2*boxes[3]+2*boxes[4]+2*boxes[5]+boxes[6])).contains(0)
    assert all(matrix[i][j].lo == matrix[j][i].lo and matrix[i][j].hi == matrix[j][i].hi for i in range(12) for j in range(12))
    for k, feature in ((3, [cos3[j % 4] for j in range(12)]), (4, [cos4[j % 3] for j in range(12)]), (5, cos5)):
        assert all((sum((matrix[i][j]*feature[j] for j in range(12)), I(0))-boxes[k]*feature[i]).contains(0) for i in range(12))
    ell = [[g*(row[0]-matrix[j][i])/2 for j in range(12)] for i in range(12)]
    ma = max(max(abs(x.lo), abs(x.hi)) for x in row)
    hmax = I.raw(0, ma)*g
    tau_lower = 11*(-I.raw(hmax.hi, hmax.hi)).exp()
    tau_upper = 11*I.raw(hmax.hi, hmax.hi).exp()
    return {"A": matrix, "diag": row[0], "g": g, "ell": ell, "tau_bound": I.raw(tau_lower.lo, tau_upper.hi), "spectral_sha256": sha(path)}


def mat(matrix, vector):
    return [sum((x*a for x, a in zip(vector, row)), 0) for row in matrix]


def intersect_value(x, interval):
    if isinstance(x, J):
        return J(x.v.intersect(interval), x.d, x.s)
    return x.intersect(interval)


def model(p, context):
    ap = mat(context["A"], p)
    hh = [x*context["g"]/2 for x in ap]
    ee = [x.exp() for x in hh]
    vv = [x/e for x, e in zip(p, ee)]
    plus, minus = sum(ee, 0), sum(vv, 0)
    tau = intersect_value(plus*minus-1, context["tau_bound"])
    drift = [e*minus-v*plus for e, v in zip(ee, vv)]
    au = mat(context["A"], ee)
    kappa = (plus*minus*context["diag"]-sum((v*x for v, x in zip(vv, au)), 0))*context["g"]/2
    hdot = [x*context["g"]/2 for x in mat(context["A"], drift)]
    tau_dot = minus*sum((e*x for e, x in zip(ee, hdot)), 0)+plus*(sum((b/e for b, e in zip(drift, ee)), 0)-sum((v*x for v, x in zip(vv, hdot)), 0))
    a0 = (11/tau).sqrt()
    a1 = (1-(1+kappa)*a0*a0+a0*a0*tau_dot/(2*tau))/(2*tau*a0)
    return {"a0": a0, "a1": a1, "tau": tau, "kappa": kappa, "drift": drift, "tau_dot": tau_dot, "ee": ee}


def phi(population, eps, context, corrected=True):
    dtype = J if isinstance(eps, J) else I
    p = [dtype(x) for x in population]
    base = model(p, context)
    av = base["a0"]+(eps*base["a1"] if corrected else 0)
    accum = eps*av
    amin = av.v.lo if isinstance(av, J) else av.lo
    for old in range(12):
        if not population[old]:
            continue
        for new in range(12):
            if new == old:
                continue
            neighbor = [x+eps*((j == new)-(j == old)) for j, x in enumerate(p)]
            nxt = model(neighbor, context)
            an = nxt["a0"]+(eps*nxt["a1"] if corrected else 0)
            amin = min(amin, an.v.lo if isinstance(an, J) else an.lo)
            c = p[old]*base["ee"][new]/base["ee"][old]
            accum += c*(eps*context["ell"][old][new]).exp()*an
    return av*accum/(11+eps), amin


def direct_field_check(population, context):
    p = [I(x) for x in population]
    data = model(p, context)
    tau, kappa, drift = I(0), I(0), [I(0) for _ in p]
    for old in range(12):
        for new in range(12):
            if old == new:
                continue
            c = p[old]*data["ee"][new]/data["ee"][old]
            tau += c
            kappa += c*context["ell"][old][new]
            drift[old] -= c
            drift[new] += c
    ab = [x*context["g"]/2 for x in mat(context["A"], drift)]
    derivative = I(0)
    for old in range(12):
        for new in range(12):
            if old != new:
                derivative += data["ee"][new]/data["ee"][old]*(drift[old]+p[old]*(ab[new]-ab[old]))
    assert (tau-data["tau"]).contains(0)
    assert (kappa-data["kappa"]).contains(0)
    assert (derivative-data["tau_dot"]).contains(0)
    assert all((x-y).contains(0) for x, y in zip(drift, data["drift"]))
    return True


POPULATIONS = {
    "uniform": [F(1, 12)]*12,
    "interior": [F(1, 8), F(1, 24)]+[F(1, 12)]*10,
    "face": [F(0)]+[F(1, 11)]*11,
    "vertex": [F(1)]+[F(0)]*11,
}


def run_case(name, context):
    p = POPULATIONS[name]
    assert sum(p) == 1 and min(p) >= 0
    checks = []
    def ck(label, value):
        assert value, (name, label)
        checks.append(label)
    ck("factorized_model_matches_direct_edge_sums", direct_field_check(p, context))
    at_zero, _ = phi(p, J(0, 1, 0), context)
    ck("analytic_origin_value_consistent", at_zero.v.contains(1))
    ck("analytic_zero_first_derivative_consistent_not_proved_by_containment", at_zero.d.contains(0))
    eps_range = I.raw(0, I(F(1, 264)).hi)
    span, positive = phi(p, J(eps_range, 1, 0), context)
    remainder = span.s
    points = []
    for n in (264, 528):
        counts = [x*n for x in p]
        ck(f"N{n}_integer_counts_and_support", all(x.denominator == 1 for x in counts) and sum(counts) == n)
        exact, amin = phi(p, I(F(1, n)), context, corrected=True)
        residual = exact-1
        scaled = residual*(n*n)
        intersected = scaled.intersect(remainder)
        ck(f"N{n}_direct_row_and_Taylor_remainder_consistent", intersected.lo <= intersected.hi)
        leading, _ = phi(p, I(F(1, n)), context, corrected=False)
        leading_residual = leading-1
        aw = model([I(x) for x in p], context)
        an = aw["a0"]+aw["a1"]/n
        diag = an*an/(11*n+1)-residual
        points.append({"N": n, "counts": [int(x) for x in counts], "corrected_row_residual": residual.box(), "N2_residual_intersection": intersected.box(), "leading_row_residual": leading_residual.box(), "trial_positive_at_all_legal_neighbors_certified": amin > 0, "analysis_self_probability_box": diag.box(), "analysis_self_probability_positive_certified": diag.lo > 0})
    return {"point": name, "population": list(map(str, p)), "gain": "exact critical lambda=0, g=12/Lambda6", "epsilon_interval": ["0", "1/264"], "checks": checks, "origin_value_box": at_zero.v.box(), "origin_derivative_box": at_zero.d.box(), "origin_half_second_derivative_box": at_zero.s.box(), "uniform_half_second_derivative_box": remainder.box(), "uniform_trial_positivity_on_epsilon_interval_certified": positive > 0, "minimum_trial_enclosure_lower": I.raw(positive, positive).box(), "direct_points": points, "scope": "actual-A7 conditional interval diagnostic, uniform in epsilon at this fixed p only; not full simplex, not independent review"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=list(POPULATIONS), required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    assert Path("/tmp") in args.output.resolve().parents
    context = native_context()
    # Small arithmetic guards are independent of the scientific row trial.
    assert (I(F(1, 3))+I(F(2, 7))).contains(F(13, 21))
    assert (I(3).sqrt()*I(3).sqrt()).contains(3)
    assert I(0).exp().contains(1)
    x = J(I(2), 1, 0)
    assert (x*x).d.contains(4) and (x*x).s.contains(1)
    data = run_case(args.case, context)
    output = {"date": "2026-10-20", "status": "PASS_FIXED_POPULATION_INTERVAL_STRESS_NOT_INDEPENDENT_AUDIT", "bits": BITS, "spectral_sha256": context["spectral_sha256"], "gain_box": context["g"].box(), "case": data}
    args.output.write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps({"status": output["status"], "point": args.case, "checks": len(data["checks"]), "remainder_box": data["uniform_half_second_derivative_box"], "trial_positive": data["uniform_trial_positivity_on_epsilon_interval_certified"]}, indent=2))


if __name__ == "__main__":
    main()
