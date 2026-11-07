"""Exact three-channel averaging and finite tail-budget controls."""
from pathlib import Path
import hashlib
import json
import math
import sympy as sp

HERE = Path(__file__).resolve().parent
checks = []


def check(name, condition):
    checks.append({'name': name, 'pass': bool(condition)})
    if not condition:
        raise AssertionError(name)


x, y = sp.symbols('x y', nonzero=True)
amp_z, amp_c, omega, lam = sp.symbols('amp_z amp_c omega lam', positive=True)
da, db, dx = sp.symbols('d_a d_b d_x', real=True)
ra, rb, rab, rba, iab, iba, rs, rd = sp.symbols('r_a r_b r_ab r_ba i_ab i_ba r_sum r_difference', real=True)
dc_a, dc_b, dc_ab, dc_ba = sp.symbols('dc_a dc_b dc_ab dc_ba', real=True)


def cosine(first, second):
    monomial = x ** first * y ** second
    return (monomial + 1 / monomial) / 2


def sine(first, second):
    monomial = x ** first * y ** second
    return (monomial - 1 / monomial) / (2 * sp.I)


pair_fa = -3 * amp_z ** 2 * (dc_a + ra * cosine(2, 0) + da * sine(2, 0))
pair_fa -= 3 * amp_c ** 2 * (dc_ab + rab * cosine(0, 2) + iab * sine(0, 2))
pair_fb = -3 * amp_c ** 2 * (dc_b + rb * cosine(0, 2) + db * sine(0, 2))
pair_fb -= 3 * amp_z ** 2 * (dc_ba + rba * cosine(2, 0) + iba * sine(2, 0))
pair_g = -6 * amp_z * amp_c * (rs * cosine(1, 1) + dx * sine(1, 1) + rd * cosine(-1, 1))
work_i = sp.expand(12 * amp_z * sine(1, 0) * (amp_z * cosine(1, 0) * pair_fa + amp_c * cosine(0, 1) * pair_g))
work_j = sp.expand(12 * amp_c * sine(0, 1) * (amp_c * cosine(0, 1) * pair_fb + amp_z * cosine(1, 0) * pair_g))
mean_i = sp.expand(work_i).coeff(x, 0).coeff(y, 0)
mean_j = sp.expand(work_j).coeff(x, 0).coeff(y, 0)
action_i = omega * amp_z ** 2 / 2
action_j = lam * amp_c ** 2 / 2
expected_i = -36 * da / omega ** 2 * action_i ** 2 - 72 * dx / (omega * lam) * action_i * action_j
expected_j = -36 * db / lam ** 2 * action_j ** 2 - 72 * dx / (omega * lam) * action_i * action_j
check('three-channel mean I work', sp.simplify(mean_i - expected_i) == 0)
check('three-channel mean J work', sp.simplify(mean_j - expected_j) == 0)
check('equal mixed action contributions', sp.diff(mean_i, dx) == sp.diff(mean_j, dx) == -18 * amp_z ** 2 * amp_c ** 2)
check('closed and real responses have zero mean work', not (mean_i.free_symbols | mean_j.free_symbols).intersection({ra, rb, rab, rba, iab, iba, rs, rd, dc_a, dc_b, dc_ab, dc_ba}))
power = 9 * omega * da * amp_z ** 4 + 9 * lam * db * amp_c ** 4 + 18 * (omega + lam) * dx * amp_z ** 2 * amp_c ** 2
check('reference energy-work normalization', sp.simplify(power + omega * mean_i + lam * mean_j) == 0)

allowed = {(0, 0), (2, 0), (-2, 0), (4, 0), (-4, 0), (0, 2), (0, -2), (0, 4), (0, -4), (2, 2), (-2, -2), (2, -2), (-2, 2)}
frequency_sets = {}
for name, expression in [('I', work_i), ('J', work_j)]:
    frequencies = set()
    for term in sp.Add.make_args(expression):
        powers = term.as_powers_dict()
        frequencies.add((int(powers.get(x, 0)), int(powers.get(y, 0))))
    check(name + ' frequency inventory', frequencies <= allowed)
    frequency_sets[name] = sorted(frequencies)
# A rational comparison proves the smallest difference exceeds 1/2.
# sqrt(3)>17/10 and sqrt(7/4)<27/20, so twice the difference exceeds 7/10.
check('uniform frequency gap rational lower bound', sp.Rational(17, 10) ** 2 < 3 and sp.Rational(27, 20) ** 2 > sp.Rational(7, 4) and 2 * (sp.Rational(17, 10) - sp.Rational(27, 20)) > sp.Rational(1, 2))

budgets = []
for denominator in [8, 16, 32, 64, 128]:
    epsilon = 1 / denominator
    horizon = epsilon ** -2
    tail_l1 = 2 * (math.sqrt(1 + horizon) - 1)
    tail_l2 = math.sqrt(math.log(1 + horizon))
    initial_transient = epsilon ** 2 * tail_l2
    history_error = tail_l1 * epsilon ** 3
    work_error = epsilon * (initial_transient + history_error)
    quartic_oscillatory_error = epsilon ** 4
    check(f'finite tail budget epsilon=1/{denominator}', work_error <= 4 * epsilon ** 3 * math.sqrt(1 + math.log(1 / epsilon)))
    budgets.append({'epsilon': epsilon, 'initial_transient_L2_bound': initial_transient, 'history_L2_bound': history_error, 'work_bound_over_epsilon_squared': work_error / epsilon ** 2, 'oscillatory_bound_over_epsilon_squared': quartic_oscillatory_error / epsilon ** 2})
check('little-o comparison ratios decrease', all(first['work_bound_over_epsilon_squared'] > second['work_bound_over_epsilon_squared'] for first, second in zip(budgets, budgets[1:])))
for source in json.loads((HERE / 'INPUTS.json').read_text())['files']:
    raw = (HERE.parent / source['path']).read_bytes()
    check('source hash ' + source['path'], hashlib.sha256(raw).hexdigest() == source['sha256'])
output = {'scope': 'Exact susceptibility work algebra, frequency inventory, finite tail bounds and provenance', 'sympy_version': sp.__version__, 'mean_I_work': str(mean_i), 'mean_J_work': str(mean_j), 'frequency_sets': frequency_sets, 'tail_budgets': budgets, 'passed': len(checks), 'failed': 0, 'checks': checks}
(HERE / 'verification.json').write_text(json.dumps(output, indent=2) + '\n')
print(json.dumps({'passed': len(checks), 'failed': 0, 'mean_I_work': str(mean_i), 'mean_J_work': str(mean_j), 'frequency_sets': frequency_sets, 'tail_budgets': budgets}, indent=2))
