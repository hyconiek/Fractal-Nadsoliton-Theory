#!/usr/bin/env python3
"""Exact finite controls, not a verifier of PDE theorems. Standard library only."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = []


def check(name, ok, scope):
    checks.append(dict(name=name, passed=bool(ok), scope=scope))
    if not ok:
        raise AssertionError(name)


class Poly:
    """Sparse exact rational polynomial in a fixed nine-variable ring."""
    n = 9

    def __init__(self, terms=None):
        self.t = {k: F(v) for k, v in (terms or {}).items() if v}

    @classmethod
    def const(cls, x):
        return cls({(0,) * cls.n: F(x)})

    @classmethod
    def var(cls, i):
        e = [0] * cls.n
        e[i] = 1
        return cls({tuple(e): 1})

    def __add__(self, other):
        other = other if isinstance(other, Poly) else Poly.const(other)
        t = dict(self.t)
        for k, v in other.t.items():
            t[k] = t.get(k, 0) + v
        return Poly(t)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.t.items()})

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        other = other if isinstance(other, Poly) else Poly.const(other)
        t = {}
        for a, av in self.t.items():
            for b, bv in other.t.items():
                k = tuple(x + y for x, y in zip(a, b))
                t[k] = t.get(k, 0) + av * bv
        return Poly(t)

    __rmul__ = __mul__

    def __pow__(self, n):
        ans = Poly.const(1)
        for _ in range(n):
            ans = ans * self
        return ans

    def diff(self, i):
        t = {}
        for k, v in self.t.items():
            if k[i]:
                new = list(k)
                new[i] -= 1
                t[tuple(new)] = v * k[i]
        return Poly(t)


def zero(p):
    return not p.t


def harmonic_product(a, b):
    out = {}
    for i, av in a.items():
        for j, bv in b.items():
            out[i+j] = out.get(i+j, Poly.const(0)) + av*bv
    return out


def run():
    Z, P, C, R, Az, Bz, Cz, w, ell = [Poly.var(i) for i in range(9)]
    Jz, Jc = P**2+w*Z**2, R**2+ell*C**2
    h = Az*Jz**2+Bz*Jz*Jc+Cz*Jc**2
    hz, hc = 2*Az*Jz+Bz*Jc, Bz*Jz+2*Cz*Jc
    for name, residual in [
        ('actual canonical Z equation', h.diff(1)-2*P*hz),
        ('actual canonical P equation', -h.diff(0)+2*w*Z*hz),
        ('actual canonical C equation', h.diff(3)-2*R*hc),
        ('actual canonical R equation', -h.diff(2)+2*ell*C*hc),
        ('coherent Jz invariant', Jz.diff(0)*h.diff(1)-Jz.diff(1)*h.diff(0)),
        ('coherent Jc invariant', Jc.diff(2)*h.diff(3)-Jc.diff(3)*h.diff(2)),
    ]:
        check(name, zero(residual), 'exact symbolic polynomial identity')
    check('negative: opposite phase sign rejected',
          not zero(-h.diff(0)-2*w*Z*hz), 'nonzero polynomial mutation')
    check('negative: P is not generally Z_t',
          not zero(h.diff(1)), 'nonzero Hamiltonian correction')

    b, bb, chi0, chip, chim, g4, *_ = [Poly.var(i) for i in range(9)]
    carrier = {-1: b*F(1,2), 1: bb*F(1,2)}
    square = harmonic_product(carrier, carrier)
    check('DC carrier normalization', zero(square[0]-b*bb*F(1,2)), 'exact Fourier algebra')
    check('second harmonic normalization', zero(square[-2]-b**2*F(1,4)), 'exact Fourier algebra')
    response = {0: square[0]*chi0, -2: square[-2]*chip, 2: square[2]*chim}
    feedback = harmonic_product(carrier, response)
    cube = harmonic_product(carrier, square)
    actual = 72*feedback[-1]-F(18,35)*g4*cube[-1]
    expected = b**2*bb*(9*(2*chi0+chip)-F(27,140)*g4)
    check('secular cubic coefficient', zero(actual-expected), 'exact harmonic coefficient before i/omega')
    check('negative: missing DC contribution detected',
          not zero(actual-b**2*bb*(9*chip-F(27,140)*g4)), 'omitted DC mutation')
    check('all other stationary frequencies separated',
          {k+1 for k in feedback if k != -1} == {-2,2,4}, 'finite Fourier support')
    check('source 3F normalization yields damping 1', F(72,8)*F(1,9)==1, 'rational Gamma factor')
    check('negative: Gamma for F is not Gamma for 3F', F(72,8)!=1, 'normalization mutation')
    check('half amplitude changes cubic coefficient by four', F(2)**3/F(2)==4, 'coordinate rescaling')

    w, ell, d2, g4, *_ = [Poly.var(i) for i in range(9)]
    delta = ell-4*w
    left = 1260*d2*(2*delta+ell)-27*g4*ell*delta
    right = -27*g4*ell*delta-1260*d2*(8*w-3*ell)
    check('shape susceptibility equals CM shift', zero(left-right),
          'CR3 multiplied by 140 omega ell Delta; exact all-parameter polynomial')
    check('negative: shape resolvent sign error detected',
          not zero(1260*d2*(2*delta-ell)-27*g4*ell*delta-right), 'sign mutation')
    # Exact gap proof: sqrt(3)<7/4, then (sqrt(3)+1/2)^2 <5.
    check('shape versus carrier gap exceeds 1/2',
          F(3)<F(7,4)**2 and F(13,4)+F(7,4)==5, 'strict rational squaring bound')
    check('quartic denominators separated',
          (3-4*F(7,4),3-4*F(5,4))==(-4,-2), 'entire interval, affine endpoints')
    for value in (F(5,4), F(3,2), F(7,4)):
        acm = 6*(2*value-3)/(3*(3-4*value))
        check('shape slip centre/endpoint '+str(value), (acm==0)==(value==F(3,2)),
              'exact coefficient control, not whole-interval sampling proof')

    check('R5 gradient degree four', zero((Z**5).diff(0)-5*Z**4), 'representative degree differentiation')
    check('R5 critical integral relative error vanishes', 4-2-1==1, 'epsilon exponent for C1 remainder')
    check('frozen resonant primitive reaches order epsilon', 3-2==1, 'scaling, not PDE dynamics')
    check('mixed cubic absolute bound does not close', 3-3==0 and 4-4==0,
          'integrals epsilon^3 T^(3/2), epsilon^4 T^2')
    check('mixed quartic bound only order epsilon', 4-3==1 and 5-4==1,
          'integrals epsilon^4 T^(3/2), epsilon^5 T^2')
    check('ST continuum relative bound not small at endpoint', 1-F(2,2)==0 and 2-2==0,
          'excluded endpoint control')
    inputs = json.loads((HERE/'INPUTS.json').read_text())
    for row in inputs['files']:
        data = (ROOT/row['path']).read_bytes()
        check('source SHA256: '+row['path'], hashlib.sha256(data).hexdigest()==row['sha256'], 'input preservation')
        git_id = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        check('pinned Git blob: '+row['path'], git_id==row['git_blob_sha'], 'pinned source identity')
    return inputs['base_commit']


if __name__ == '__main__':
    try:
        commit = run()
        status = 'PASS_FINITE_CONTROLS_CRITICAL_PDE_GATE_OPEN'
    except Exception as exc:
        commit = None
        status = 'FAIL: '+str(exc)
    result = dict(status=status, base_commit=commit, total_checks=len(checks),
                  passed=sum(c['passed'] for c in checks), checks=checks,
                  proof_assistant_verification=False, independent_review=False,
                  actual_critical_damping_law_proved=False)
    (HERE/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'}, indent=2))
    raise SystemExit(0 if status.startswith('PASS') else 1)
