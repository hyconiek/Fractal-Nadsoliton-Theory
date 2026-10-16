#!/usr/bin/env python3
"""Small exact checks for the research submission, not physical validation.

Uses only the Python standard library. The general proofs are in
THEORETICAL_RESULTS.md; finite fixtures do not prove those theorems.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import math

HERE = Path(__file__).resolve().parent


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col))
             for col in zip(*b)] for row in a]


def rank(a):
    a = [[F(x) for x in row] for row in a]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        d = a[r][c]
        a[r] = [x / d for x in a[r]]
        for i in range(len(a)):
            if i != r:
                d = a[i][c]
                a[i] = [x - d * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def shifted_hamiltonian(p, s, center_exp, tilt_exp):
    # exp(center + tilt), with exact rational exponential coordinates.
    z = [a * b for a, b in zip(center_exp, tilt_exp)]
    return sum(p[i] * s[j] * (z[j] / z[i] - 1)
               for i in range(len(p)) for j in range(len(p)))


def centre_checks():
    q = 12
    # Native g=0, off-equilibrium p_i=(i+1)^2/650, s_i=1/12.
    # exp(dF_i/2) is proportional to i+1; the common factor cancels.
    r = [F(i) for i in range(1, q + 1)]
    p = [x*x / sum(y*y for y in r) for x in r]
    s = [F(1, q)] * q
    tilts = [
        [F(2), F(1, 2)] + [F(1)] * 10,
        [F(i + 2, i + 1) for i in range(q)],
        [F(3) if i % 2 else F(1, 3) for i in range(q)],
    ]
    residuals = []
    for z in tilts:
        hp = shifted_hamiltonian(p, s, r, z)
        hm = shifted_hamiltonian(p, s, r, [1 / x for x in z])
        residuals.append(hp - hm)
        assert hp == hm
    weights = [[p[i] * s[j] * r[j] / r[i] for j in range(q)]
               for i in range(q)]
    assert all(weights[i][j] == weights[j][i]
               for i in range(q) for j in range(q))
    hess = [[F(0) for _ in range(q)] for _ in range(q)]
    for i in range(q):
        for j in range(i + 1, q):
            w = 2 * weights[i][j]
            hess[i][i] += w
            hess[j][j] += w
            hess[i][j] -= w
            hess[j][i] -= w
    assert rank(hess) == 11
    # Proposed W=e0-e1, xi=log(2)*(e0-e1) demands a nonzero
    # generalized-symmetry RHS -4*log(2); the computed LHS is zero.
    # log(2)>1/2 follows from integral_1^2 dx/x>1/2.
    return {
        "scope": "exact g=0 q12 fixture; general positive-g proof is analytical",
        "center_symmetry_residuals": [str(x) for x in residuals],
        "hessian_rank": rank(hess),
        "nonzero_W_negative_control": "required RHS=-4 log(2), actual LHS=0",
        "negative_control_absolute_residual_lower_bound": "2 (strict)",
    }


def energy_lift(a):
    degeneracies = [1, 2, 4]
    energies = [2, 1, 0]
    states = [(x, b) for x, d in enumerate(degeneracies) for b in range(d)]
    q = [[F(0) for _ in range(3)] for _ in range(3)]
    for x in range(3):
        for y in range(3):
            if x != y:
                q[x][y] = F(degeneracies[y]) * a[x][y]
        q[x][x] = -sum(q[x])
    g = [[F(0) for _ in states] for _ in states]
    for i, (x, _) in enumerate(states):
        for j, (y, _) in enumerate(states):
            if x != y:
                g[i][j] = a[x][y]
        g[i][i] = -sum(g[i])
    jmat = [[F(int(x == y)) for y in range(3)] for x, _ in states]
    assert matmul(g, jmat) == matmul(jmat, q)
    assert all(g[i][j] == g[j][i] for i in range(7) for j in range(7))
    assert all(sum(row) == 0 for row in g + q)
    pi = [F(d, 7) for d in degeneracies]
    assert all(pi[x] * q[x][y] == pi[y] * q[y][x]
               for x in range(3) for y in range(3))
    assert all(sum(pi[x] * q[x][y] for x in range(3)) == 0 for y in range(3))
    totals = [energies[x] + x for x, _ in states]
    assert totals == [2] * 7
    assert rank(q) == 2
    return q, pi


def energy_checks():
    a = [[F(int(i != j)) for j in range(3)] for i in range(3)]
    b = [row[:] for row in a]
    b[0][1] = b[1][0] = F(2)
    qa, pi = energy_lift(a)
    qb, pi_b = energy_lift(b)
    assert pi == pi_b
    assert qb[0][1] / qa[0][1] == 2
    assert qb[0][2] / qa[0][2] == 1
    assert qa != qb
    return {
        "microstates": 7,
        "microcanonical_stationarity": "uniform",
        "projected_pi": [str(x) for x in pi],
        "system_energies": [2, 1, 0],
        "bath_energies": [0, 1, 2],
        "total_energy_on_selected_shell": 2,
        "Q_A": [[str(x) for x in row] for row in qa],
        "Q_B": [[str(x) for x in row] for row in qb],
        "two_rate_rescaling_ratios": ["2", "1"],
        "intertwining_and_detailed_balance": "exact",
        "scope": "auxiliary CTMC countermodel; no autonomous Hamiltonian or FIN source",
    }


def alternating_enclosure(terms):
    # Terms decrease in absolute value for the two arguments used below.
    partial = F(0)
    previous = None
    for i, term in enumerate(terms):
        if previous is not None:
            assert abs(term) <= previous
        previous = abs(term)
        partial += term
        if i == 20:
            upper = partial
        if i == 21:
            lower = partial
    assert lower <= upper
    return lower, upper


def dilation_checks():
    q = 12
    generator = [[F(2, q) - (2 if i == j else 0)
                  for j in range(q)] for i in range(q)]
    f = [F(q - 1)] + [F(-1)] * (q - 1)
    af = [-sum(generator[i][j] * f[j] for j in range(q)) for i in range(q)]
    assert af == [2 * x for x in f]
    norm2 = sum(x*x for x in f) / q
    a = sum(x*y for x, y in zip(f, af)) / q / norm2
    b = sum(x*x for x in af) / q / norm2
    assert (a, b) == (F(2), F(4))
    sigma2 = F(16)
    t = a / (b + sigma2)
    bound = a*t - (b + sigma2)*t*t/2
    assert t == F(1, 10) and bound == F(1, 10)
    exp_box = alternating_enclosure([
        (-1)**n * F(1, 5)**n / math.factorial(n) for n in range(22)])
    cos_box = alternating_enclosure([
        (-1)**n * F(2, 5)**(2*n) / math.factorial(2*n) for n in range(22)])
    defect = (cos_box[0] - exp_box[1], cos_box[1] - exp_box[0])
    assert defect[0] > bound
    return {
        "reference": "native N2 heat-bath difference eigenmode, rate 2; g=0 fixture",
        "a": str(a), "b": str(b), "sigma_squared": str(sigma2),
        "time": str(t), "universal_defect_lower_bound": str(bound),
        "rotation_countermodel_defect_exact_enclosure": [str(x) for x in defect],
        "rotation_countermodel_defect_float_preview": [float(x) for x in defect],
        "scope": "correlation bound with regular unitary dilation, not empirical power",
    }


def main():
    result = {
        "date": "2026-10-16",
        "status": "PASS_EXACT_FINITE_CHECKS",
        "arithmetic": "Python Fraction; no binary64 decisions",
        "centre_symmetry": centre_checks(),
        "energy_conservation_countermodel": energy_checks(),
        "regular_dilation_bound": dilation_checks(),
        "composition_type_check": {
            "two_q12_state_sets_disjoint_union": 24,
            "two_q12_simultaneous_states_cartesian_product": 144,
            "meaning": "different composition objects; cardinalities alone are not a locality theorem",
        },
        "scientific_bridge": "OPEN_INPUT",
        "KB_001": "NOT_EXECUTED_AS_BRIDGE_ADMISSION",
        "KB_002": "NOT_ACTIVATED",
        "KB_003": "NOT_ACTIVATED",
        "physical_validation": False,
    }
    (HERE / "RESULTS.json").write_text(json.dumps(result, indent=2) + "\n")
    print(result["status"])
    print("centre symmetry: exact; energy lift: exact; regular dilation bound: exact")
    print("Bridge: OPEN_INPUT; no microscopic law adopted")


if __name__ == "__main__":
    main()
