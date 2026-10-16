#!/usr/bin/env python3
"""Exact kernel/rate-contrast certificate plus separately labelled float checks.

Only the standard library is used. General identities are proved in PROOF.md.
No simulation, fitting, external data, or modified FIN generator is adopted.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
import hashlib
import itertools
import json
import math

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    def __post_init__(self):
        assert self.lo <= self.hi

    @staticmethod
    def of(x):
        return x if isinstance(x, I) else I(F(x), F(x))

    def __add__(self, other):
        o = I.of(other)
        return I(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-I.of(other))

    def __rsub__(self, other):
        return I.of(other) + (-self)

    def __mul__(self, other):
        o = I.of(other)
        p = [x*y for x in (self.lo, self.hi) for y in (o.lo, o.hi)]
        return I(min(p), max(p))

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = I.of(other)
        assert o.lo > 0 or o.hi < 0
        return self * I(1/o.hi, 1/o.lo)


def public(box, digits=24):
    scale = 10**digits
    low = box.lo.numerator * scale // box.lo.denominator
    high = -((-box.hi.numerator * scale) // box.hi.denominator)

    def dec(n):
        sign = '-' if n < 0 else ''
        n = abs(n)
        return sign + str(n//scale) + '.' + str(n % scale).zfill(digits)

    out = [dec(low), dec(high)]
    assert F(out[0]) <= box.lo <= box.hi <= F(out[1])
    return out


def root_interval(integer, degree, digits=55):
    scale = 10**digits
    target = integer * scale**degree
    low, high = 0, 1 << ((target.bit_length() + degree - 1)//degree)
    while high - low > 1:
        mid = (low + high)//2
        if mid**degree <= target:
            low = mid
        else:
            high = mid
    assert low**degree <= target < (low+1)**degree
    return I(F(low, scale), F(low if low**degree == target else low+1, scale))


def cos_rational(x):
    # All declared arguments satisfy 0<x<1.3. Alternating terms decrease
    # starting at n=0. Consecutive even/odd partial sums enclose cos(x).
    assert 0 <= x <= F(13, 10)
    total = F(0)
    previous = F(2)
    upper = lower = None
    for n in range(34):
        term = x**(2*n)/math.factorial(2*n)
        assert term <= previous
        previous = term
        total += (-1)**n * term
        if n == 32:
            upper = total
        if n == 33:
            lower = total
    return I(lower, upper)


def build_certificate():
    sqrt3 = root_interval(3, 2)
    c = [I.of(1), sqrt3/2, I.of(F(1, 2)), I.of(0), I.of(F(-1, 2)),
         -sqrt3/2, I.of(-1), -sqrt3/2, I.of(F(-1, 2)), I.of(0),
         I.of(F(1, 2)), sqrt3/2]
    weights = {}
    for d in range(1, 7):
        angle = F('0.18575')*d + F('0.1625')
        power = root_interval(d**9, 5)
        weights[d] = cos_rational(angle)/(1 + power)
        assert weights[d].lo > 0
    lam = {k: sum((2 if d < 6 else 1)*weights[d]*(1-c[(k*d) % 12])
                  for d in range(1, 7)) for k in (3, 4, 5, 6)}
    prior_path = ROOT/'fin_rank7_followup/inputs/fin_handoff_audit/results.json'
    prior = json.loads(prior_path.read_text())['exact']['laplacian_intervals']
    for k, box in lam.items():
        assert F(prior[k][0]) <= box.lo <= box.hi <= F(prior[k][1])
    row = [
        (2*sum(lam[k]*c[(k*d) % 12] for k in (3, 4, 5))
         + lam[6]*((-1)**d))/12
        for d in range(12)
    ]
    gap3 = (lam[3] + lam[5] + lam[6])/6
    gap6 = (lam[3] + lam[5])/3
    ratio = (1 + lam[6]/(lam[3] + lam[5]))/2
    assert gap3.lo > 0 and gap6.lo > 0
    for d, gap in [(3, gap3), (6, gap6)]:
        direct = row[0]-row[d]
        assert max(direct.lo, gap.lo) <= min(direct.hi, gap.hi)
    # Conservative public interval, independently checked against raw endpoints.
    promised = ["0.77490314788261", "0.77490314788262"]
    assert F(promised[0]) <= ratio.lo <= ratio.hi <= F(promised[1])
    assert ratio.lo > F(3, 4)
    certificate = {
        'status': 'EXACT_RATIONAL_INPUT_TO_OUTPUT_ENCLOSURE',
        'kernel': {'omega': '0.18575', 'phi': '0.1625', 'eta': '9/5',
                   'weight': 'cos(omega*d+phi)/(1+d^(9/5))',
                   'retained_modes': [3, 4, 5, 6, 7, 8, 9]},
        'method': 'rational alternating cosine series, integer-certified fifth/square roots, exact interval arithmetic; no mpmath/libm decisions',
        'laplacian_active_spectra': {str(k): public(v) for k, v in lam.items()},
        'gap_A00_minus_A03': public(gap3),
        'gap_A00_minus_A06': public(gap6),
        'log_R3_over_log_R6': public(ratio),
        'conservative_public_ratio': promised,
        'flat_rank7_ratio': '3/4',
        'separation_from_flat_rank7': public(ratio-F(3, 4)),
        'prior_spectrum_sha256': hashlib.sha256(prior_path.read_bytes()).hexdigest(),
        'scope': 'exact declared mathematical kernel, not nature selection, rate measurement or physical calibration',
    }
    return certificate, [float((x.lo+x.hi)/2) for x in row]


def exact_mixed_check():
    # Deliberately nonconstant diagonal and slot-dependent fixed fields.
    # Test the general finite Gibbs identity independently of A7.
    q, n, gain = 4, 3, F(7, 5)
    mat = [[F((i+1)*(j+1), 11) + (F(i+2, 7) if i == j else 0)
            for j in range(q)] for i in range(q)]
    fields = [[F((slot+1)*(label+2), 13) for label in range(q)] for slot in range(n)]

    def potential(x):
        return gain/(2*n)*sum(mat[a][b] for a in x for b in x) + sum(
            fields[i][a] for i, a in enumerate(x))

    count = 0
    for a, b, c, d, z in itertools.product(range(q), repeat=5):
        if a == b or c == d:
            continue
        mixed = potential((b, d, z))-potential((a, d, z))-potential((b, c, z))+potential((a, c, z))
        target = gain/n*(mat[b][d]-mat[a][d]-mat[b][c]+mat[a][c])
        assert mixed == target
        count += 1
    return {'fixtures': count, 'arithmetic': 'Fraction',
            'scope': 'finite algebra fixtures; all-N identity is proved analytically'}


def numerical_rates(row):
    # Independent finite labelled checks from actual single-site proposals.
    q, gain = 12, 3.0
    a = [[row[(j-i) % q] for j in range(q)] for i in range(q)]

    def rate(x, slot, new, kind):
        old = x[slot]
        assert old != new
        field = [gain/len(x)*sum(a[b][v] for j, v in enumerate(x) if j != slot)
                 for b in range(q)]
        delta = field[new]-field[old]
        if kind in ('HB', 'HB_THINNED'):
            shift = max(field)
            e = [math.exp(v-shift) for v in field]
            value = e[new]/sum(e)
            distance = min((new-old) % q, (old-new) % q)
            if kind == 'HB_THINNED' and distance == 3:
                value /= 1+gain*gain
            return value
        if kind == 'MET':
            return min(1.0, math.exp(delta))/q
        if kind == 'BAR':
            # Overall normalization is immaterial to every ratio here.
            return 1.0/(q*(1+math.exp(-delta)))
        raise ValueError(kind)

    def mixed(d, z, kind):
        x00, xd0, x0d, xdd = (0, 0, z), (d, 0, z), (0, d, z), (d, d, z)
        return (math.log(rate(x0d, 0, d, kind))-math.log(rate(xdd, 0, 0, kind))
                - math.log(rate(x00, 0, d, kind))+math.log(rate(xd0, 0, 0, kind)))

    records = []
    target_ratio = (row[0]-row[3])/(row[0]-row[6])
    for kind in ('HB', 'MET', 'BAR', 'HB_THINNED'):
        for z in (0, 7):
            x3, x6 = mixed(3, z, kind), mixed(6, z, kind)
            target3, target6 = 2*gain/3*(row[0]-row[3]), 2*gain/3*(row[0]-row[6])
            error = max(abs(x3-target3), abs(x6-target6), abs(x3/x6-target_ratio))
            assert error < 3e-14
            records.append({'kinetics': kind, 'spectator': z, 'C3': x3, 'C6': x6,
                            'ratio': x3/x6, 'max_error': error})
    x = (0, 0, 7)
    outgoing_hb = rate(x, 0, 3, 'HB')/rate(x, 0, 6, 'HB')
    outgoing_thinned = rate(x, 0, 3, 'HB_THINNED')/rate(x, 0, 6, 'HB_THINNED')
    assert abs(outgoing_thinned/outgoing_hb-0.1) < 1e-14
    for slot in range(3):
        attempt = sum(rate(x, slot, b, 'HB_THINNED') for b in range(q) if b != x[slot])
        assert 0 < attempt < 1

    # Count-rate negative control: a repeated 0->3 change, N=4.
    x0 = (0, 0, 3, 7)
    x1 = (3, 0, 3, 7)
    x2 = (3, 3, 3, 7)
    forward0 = 2*rate(x0, 0, 3, 'HB')
    reverse0 = 2*rate(x1, 0, 0, 'HB')
    forward1 = rate(x1, 1, 3, 'HB')
    reverse1 = 3*rate(x2, 1, 0, 'HB')
    raw = math.log(forward1/reverse1)-math.log(forward0/reverse0)
    desired = 2*gain/4*(row[0]-row[3])
    corrected = raw + math.log(3)
    assert abs(corrected-desired) < 3e-14
    assert abs(raw-desired+math.log(3)) < 3e-14
    return {
        'status': 'BINARY64_DIAGNOSTICS_NOT_INTERVAL_PROOFS',
        'records': records,
        'same_source_ratio_HB': outgoing_hb,
        'same_source_ratio_thinned': outgoing_thinned,
        'same_source_ratio_change': outgoing_thinned/outgoing_hb,
        'count_multiplicity_negative_control': {
            'raw': raw, 'corrected': corrected, 'target': desired,
            'raw_bias_exact': '-log(3)',
        },
        'thinning_rule': 'h(d,g)=1/(1+g^2) for internal label distance 3, otherwise 1; diagnostic countermodel only',
    }


def main():
    certificate, row = build_certificate()
    algebra = exact_mixed_check()
    numeric = numerical_rates(row)
    (HERE/'CERTIFICATE.json').write_text(json.dumps(certificate, indent=2)+'\n')
    result = {
        'date': '2026-10-16',
        'status': 'PASS_SCOPED_RATE_COUPLING_DERIVATION',
        'scientific_status': 'DERIVED_CONDITIONAL_RATE_RELATION; KINETIC_SOURCE_OPEN',
        'exact_algebra': algebra,
        'certificate': 'CERTIFICATE.json',
        'numerical_checks': numeric,
        'KB_bridge_admission': False,
        'new_microscopic_law_adopted': False,
        'empirical_validation': False,
    }
    (HERE/'RESULTS.json').write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
    print('log(R3)/log(R6):', certificate['conservative_public_ratio'])
    print('Exact rational Gibbs fixtures:', algebra['fixtures'])
    print('Finite kinetic diagnostics:', len(numeric['records']))


if __name__ == '__main__':
    main()
