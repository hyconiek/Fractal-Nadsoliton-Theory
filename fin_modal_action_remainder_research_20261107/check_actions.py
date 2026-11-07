"""Exact modal Poisson algebra, mixed-degree audit and source provenance."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
checks = []


def check(name, condition):
    checks.append({'name': name, 'pass': bool(condition)})
    if not condition:
        raise AssertionError(name)


def zero(name, expression):
    check(name, sp.simplify(expression) == 0)


z, p, c, r = sp.symbols('z p c r', real=True)
omega, d = sp.symbols('omega d', positive=True)
lam = sp.sqrt(3)
variables = (z, p, c, r)


def bracket(first, second):
    return (sp.diff(first, z) * sp.diff(second, p) - sp.diff(first, p) * sp.diff(second, z)
            + sp.diff(first, c) * sp.diff(second, r) - sp.diff(first, r) * sp.diff(second, c))


I = (p ** 2 + omega ** 2 * z ** 2) / (2 * omega)
J = (r ** 2 + 3 * c ** 2) / (2 * lam)
H2 = omega * I + lam * J
Delta = 3 - 4 * omega ** 2
A = 2 * (2 * omega ** 2 - 3) / Delta
B = 12 / Delta
C3 = 4 / Delta
E3 = -sp.Rational(4, 3)
F3 = -sp.Rational(8, 27)
G3 = d * (A * z ** 2 * r + B * z * p * c + C3 * p ** 2 * r + E3 * c ** 2 * r + F3 * r ** 3)
H3b = d * (6 * z ** 2 * c + 4 * c ** 3)
zero('CM homological equation and flow sign', bracket(G3, H2) - H3b)
check('opposite flow would retain a cubic term', sp.simplify(H3b - bracket(H2, G3)) != 0)
correction_I = d * (B * omega * z ** 2 * c + (2 * omega * C3 - 2 * A / omega) * z * p * r - B / omega * p ** 2 * c)
correction_J = d * (lam * A * z ** 2 * c + lam * C3 * p ** 2 * c + lam * E3 * c ** 3 + (3 * lam * F3 - 2 * E3 / lam) * c * r ** 2 - B / lam * z * p * r)
zero('explicit cubic I corrector', bracket(I, G3) - correction_I)
zero('explicit cubic J corrector', bracket(J, G3) - correction_J)
zero('energy-weighted correctors', omega * correction_I + lam * correction_J + H3b)
zero('exceptional initial shape slip', A.subs(omega, sp.sqrt(sp.Rational(3, 2))))

fa_pair, fb_pair, mixed_pair = sp.symbols('fa_pair fb_pair mixed_pair', real=True)
Hrad = 6 * z ** 2 * fa_pair + 6 * c ** 2 * fb_pair + 12 * z * c * mixed_pair
zero('all-channel I work coefficients', bracket(I, Hrad) + 12 * p / omega * (z * fa_pair + c * mixed_pair))
zero('all-channel J work coefficients', bracket(J, Hrad) + 12 * r / lam * (c * fb_pair + z * mixed_pair))
az, bz, cz = sp.symbols('az bz cz')
normal_quartic = az * (2 * omega * I) ** 2 + bz * (2 * omega * I) * (2 * lam * J) + cz * (2 * lam * J) ** 2
zero('quartic normal form preserves I', bracket(I, normal_quartic))
zero('quartic normal form preserves J', bracket(J, normal_quartic))

K, f, g, u, v, g4 = sp.symbols('K f g u v g4', real=True)
eta = c * g + v
psi = z * f + u
cubic_density = 2 * K * eta ** 3 + 6 * K * eta * psi ** 2
bound_density = cubic_density.subs({u: 0, v: 0})
rad_density = 6 * z ** 2 * K * f ** 2 * v + 6 * c ** 2 * K * g ** 2 * v + 12 * z * c * K * f * g * u
nonlinear_density = 6 * c * K * g * (v ** 2 + u ** 2) + 12 * z * K * f * v * u
pure_density = 2 * K * v ** 3 + 6 * K * v * u ** 2
zero('complete cubic density expansion', cubic_density - bound_density - rad_density - nonlinear_density - pure_density)
zero('two-continuum I remainder', bracket(I, nonlinear_density) + 12 * p / omega * K * f * v * u)
zero('two-continuum J remainder', bracket(J, nonlinear_density) + 6 * r / lam * K * g * (v ** 2 + u ** 2))
zero('pure-continuum modal bracket', bracket(I, pure_density) + bracket(J, pure_density))
quartic_density = eta ** 4 / 2 + 3 * eta ** 2 * psi ** 2 + g4 * psi ** 4 / 4
mixed_quartic = quartic_density - quartic_density.subs({u: 0, v: 0}) + bracket(cubic_density - bound_density, G3)
degree_counts = {}
for label, action in [('I', I), ('J', J)]:
    polynomial = sp.Poly(sp.expand(bracket(action, mixed_quartic)), z, p, c, r, u, v)
    counts = {}
    localized = True
    allowed = True
    for powers, coefficient in polynomial.terms():
        continuum_degree = powers[4] + powers[5]
        modal_degree = sum(powers[:4])
        counts[continuum_degree] = counts.get(continuum_degree, 0) + 1
        allowed = allowed and 1 <= continuum_degree <= 3 and modal_degree + continuum_degree == 4
        localized = localized and sp.simplify(coefficient.subs({f: 0, g: 0})) == 0
    check(f'{label} quartic degree classification', allowed)
    check(f'{label} quartic coefficient localization', localized)
    degree_counts[label] = counts

# Saturating norm histories test the integrated exponent bookkeeping.
scalings = []
for denominator in [8, 16, 32, 64]:
    epsilon = Fraction(1, denominator)
    T = epsilon ** -2
    local_size = epsilon ** 2
    local_square_integral = T * local_size ** 2
    cubic_error = epsilon * local_square_integral
    quartic_error = epsilon ** 3 * T * local_size + epsilon ** 2 * local_square_integral
    taylor_error = epsilon ** 5 * T
    check(f'critical norm budget epsilon=1/{denominator}', local_square_integral == epsilon ** 2)
    check(f'action remainder order epsilon=1/{denominator}', cubic_error == epsilon ** 3 and taylor_error == epsilon ** 3 and quartic_error <= 2 * epsilon ** 3)
    scalings.append({'epsilon': float(epsilon), 'cubic_error_over_epsilon_squared': float(cubic_error / epsilon ** 2), 'quartic_error_over_epsilon_squared': float(quartic_error / epsilon ** 2), 'taylor_error_over_epsilon_squared': float(taylor_error / epsilon ** 2)})
for source in json.loads((HERE / 'INPUTS.json').read_text())['files']:
    raw = (HERE.parent / source['path']).read_bytes()
    check('source hash ' + source['path'], hashlib.sha256(raw).hexdigest() == source['sha256'])
output = {'scope': 'Exact Poisson and density algebra, integrated-power controls and source provenance', 'sympy_version': sp.__version__, 'mixed_quartic_monomial_counts': degree_counts, 'norm_history_scalings': scalings, 'passed': len(checks), 'failed': 0, 'checks': checks}
(HERE / 'verification.json').write_text(json.dumps(output, indent=2) + '\n')
print(json.dumps({'passed': len(checks), 'failed': 0, 'mixed_quartic_monomial_counts': degree_counts, 'norm_history_scalings': scalings}, indent=2))
