#!/usr/bin/env python3
"""Exact coefficient, sparse-jump and integrity controls, not a path proof."""
import argparse
import ast
import hashlib
import json
import re
import subprocess
from fractions import Fraction as F
from itertools import product
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "fin_sqrt_local_regime_research_20261020/NEXT_CRITICAL_ORDER.md",
    "fin_sqrt_local_regime_research_20261020/CONTRACT.md",
    "fin_sqrt_local_regime_research_20261020/LOCAL_WINDOW_THEOREM.md",
    "fin_sqrt_local_regime_research_20261020/FINITE_SOURCE_DIRECT_THEOREM.md",
    "fin_spectator_limit_research_20261020/SCALING_AND_GENERATOR_PROOF.md",
    "fin_spectator_limit_research_20261020/PATH_LIMIT_AND_TV_SCOPE.md",
    "fin_reference_transport_research_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md",
    "fin_hyp_history_entropy_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_hyp_history_entropy_review_20261020/PROOF_AUDIT_AND_DUAL.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_critical_metastability_campaign_20261014/MT-003_CONDITIONAL_GIBBS.md",
    "fin_uniform_critical_campaign_20261014/UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md",
    "fin_uniform_critical_campaign_20261014/UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md",
    "fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md",
    "fin_stationary_kinetic_rigidity_review_20261016/PROOF_COMPLETIONS_AND_SCOPE.md",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
]
PROTECTED = [
    "fin_sqrt_local_regime_research_20261020",
    "fin_spectator_limit_research_20261020",
    "fin_spectator_composition_research_20261020",
    "fin_reference_transport_research_20261020",
    "fin_hyp_history_entropy_20261020",
    "fin_hyp_history_entropy_review_20261020",
    "fin_critical_metastability_review_20261015",
    "fin_critical_metastability_campaign_20261014",
    "fin_uniform_critical_campaign_20261014",
    "fin_uniform_critical_review_20261014",
    "fin_intrinsic_observation_review_20261015",
    "fin_stationary_kinetic_rigidity_review_20261016",
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder):
    seen = set()
    for line in (folder / "MANIFEST.sha256").read_text().splitlines():
        h, name = line.split("  ", 1)
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder / name) == h, (folder.name, name)
        seen.add(name)
    actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
    actual |= {n for n in seen if "__pycache__" in Path(n).parts or n.endswith(".pyc")}
    assert actual == seen | {"MANIFEST.sha256"}, (folder.name, actual - seen, seen - actual)
    return {"entries": len(seen), "manifest_sha256": sha(folder / "MANIFEST.sha256"), "coverage": "complete_noncache"}


def mul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), F(0)) for j in range(len(b[0]))] for i in range(len(a))]


def matvec(a, x):
    return [sum((t * s for t, s in zip(row, x)), F(0)) for row in a]


def rowvec(x, a):
    return [sum((x[i] * a[i][j] for i in range(len(x))), F(0)) for j in range(len(x))]


def add_term(d, key, value):
    d[key] = d.get(key, F(0)) + value


def clean(d):
    return {k: v for k, v in d.items() if v}


def line_laurent(q, sigma, jump_power):
    out = {}
    for i, j in product(range(q), repeat=2):
        ds = sigma[j] - sigma[i]
        add_term(out, (ds, 0), F(ds ** jump_power, q))
        add_term(out, (ds, 1), F(sigma[i] * ds ** jump_power, q))
    return clean(out)


def fast_laurent(a, sigma, w, g, order):
    q = len(sigma)
    aw = matvec(a, w)
    out = {}
    for i, j in product(range(q), repeat=2):
        ds = sigma[j] - sigma[i]
        derivative = g * (aw[j] - aw[i]) / 2
        if order == 1:
            const = w[i] + derivative / q
            linear = sigma[i] * derivative / q
        else:
            const = w[i] * derivative + derivative ** 2 / (2 * q)
            linear = sigma[i] * derivative ** 2 / (2 * q)
        add_term(out, (ds, 0), ds * const)
        add_term(out, (ds, 1), ds * linear)
    return clean(out)


def line_series(q, k, bracket=False):
    # Direct expansion of the two opposite-parity edge groups, not sinh code.
    out = [F(0)] * 8
    for power in range(8):
        for s in (-1, 1):
            pref = F(q) if bracket else F(-q * s, 2)
            out[power] += pref * F((-s * k) ** power, factorial(power))
            if power:
                out[power] += pref * s * F((-s * k) ** (power - 1), factorial(power - 1))
    return out


def block_law(mu, p, jumps):
    out = {}
    for word in product(range(len(mu)), repeat=jumps + 1):
        weight = mu[word[0]]
        for i, j in zip(word, word[1:]):
            weight *= p[i][j]
        out[word] = weight
    return out


def warp(t, knots):
    for (x, y), (xx, yy) in zip(knots, knots[1:]):
        if x <= t <= xx:
            return y + (t-x) * (yy-y) / (xx-x)
    raise ValueError(t)


def at(t, times, values):
    return values[sum(s <= t for s in times)]


def controls():
    checks = []

    def ck(name, ok):
        assert name not in checks and ok, name
        checks.append(name)

    q = 12
    sig = [(-1) ** j for j in range(q)]
    p6 = [[F(sig[i] * sig[j], q) for j in range(q)] for i in range(q)]
    cos3 = [F(1), F(0), F(-1), F(0)]
    cos4 = [F(1), F(-1, 2), F(-1, 2)]
    p3 = [[F(2, q) * cos3[(i-j) % 4] for j in range(q)] for i in range(q)]
    p4 = [[F(2, q) * cos4[(i-j) % 3] for j in range(q)] for i in range(q)]
    a = [[p6[i][j] + p3[i][j] / 2 + p4[i][j] / 3 for j in range(q)] for i in range(q)]
    zero = [[F(0)] * q for _ in range(q)]
    ck("rational_symmetry_fixture_exact_projectors", all(mul(p, p) == p for p in (p6, p3, p4)) and all(mul(p, r) == zero for p, r in ((p6, p3), (p6, p4), (p3, p4))))
    ck("rational_fixture_rank_five_not_A7", [sum(p[i][i] for i in range(q)) for p in (p6, p3, p4)] == [1, 2, 2])
    ck("rational_fixture_circulant_zero_mean_top_parity", matvec(a, [F(1)] * q) == [0] * q and matvec(a, sig) == sig and all(a[(i+2) % q][(j+2) % q] == a[i][j] for i, j in product(range(q), repeat=2)))
    expected_b = {(2, 0): F(q, 2), (2, 1): F(-q, 2), (-2, 0): F(-q, 2), (-2, 1): F(-q, 2)}
    expected_g = {(2, 0): F(q), (2, 1): F(-q), (-2, 0): F(q), (-2, 1): F(q)}
    ck("exact_two_group_hyperbolic_drift_identity", line_laurent(q, sig, 1) == expected_b)
    ck("exact_two_group_hyperbolic_bracket_identity", line_laurent(q, sig, 2) == expected_g)
    frame = []
    for parity in (0, 1):
        for j in range(parity + 2, q, 2):
            w = [F(0)] * q
            w[parity], w[j] = F(1), F(-1)
            frame.append(w)
    ck("all_ten_fast_directions_have_zero_shift_average", len(frame) == 10 and all(sum(w) == 0 and sum(s * x for s, x in zip(sig, w)) == 0 and all(sum(w[(i+2*k) % q] for k in range(6)) == 0 for i in range(q)) for w in frame))
    ck("first_fast_Laurent_derivative_vanishes_identically", all(fast_laurent(a, sig, w, g, 1) == {} for w in frame for g in (F(12), F(131, 10))))
    w = [cos3[j % 4] for j in range(q)]
    second = fast_laurent(a, sig, w, F(12), 2)
    at_origin = sum(value for (power, mdegree), value in second.items() if mdegree == 0)
    ck("nonzero_fast_quadratic_term_not_discarded", at_origin == 216)
    ck("fast_fixture_eigenvalue_is_transverse", matvec(a, w) == [x / 2 for x in w])
    for k in (F(0), F(1), F(3, 2)):
        b, bracket = line_series(q, k), line_series(q, k, True)
        expected = [F(0)] * 8
        expected[1] = q * (k - 1)
        for power in (3, 5, 7):
            expected[power] = q * (k ** power / factorial(power) - k ** (power-1) / factorial(power-1))
        ck(f"parity_line_Taylor_coefficients_k{k}", b == expected)
        ck(f"parity_bracket_even_coefficients_k{k}", bracket[0] == 2*q and all(bracket[j] == 0 for j in (1, 3, 5, 7)) and bracket[2] == 2*q*(k*k/2-k))
    ck("critical_cubic_and_noise_keep_actual_clock", line_series(q, F(1))[3] == -4 and line_series(q, F(1))[5] == F(-2, 5) and line_series(q, F(1), True)[0] == 24)
    # Exact current/divergence coefficients for rho'/rho=a*z-z^3/3.
    aa = F(2, 5)
    logderiv = {1: aa, 3: F(-1, 3)}
    drift = {degree: q * coef for degree, coef in logderiv.items()}
    ck("quartic_zero_current_matches_diffusion_q", drift == {1: F(24, 5), 3: F(-4)})
    ck("lambda_zero_stationary_fourth_moment_identity", -F(2*q, 3) * 3 + 2*q == 0)
    for x, j in ((F(2, 3), F(1, 7)), (F(-3, 5), F(-2, 11))):
        ck(f"fourth_jump_compensator_algebra_{x}_{j}", (x+j)**4-x**4-4*x**3*j == 6*x*x*j*j+4*x*j**3+j**4)
    # A genuine labelled/cavity/count jump fixture, not a critical domain.
    states = list(product(range(q), repeat=2))
    amp = [2 if x == y else 1 for x, y in states]
    index = {x: i for i, x in enumerate(states)}
    counts = lambda x: tuple(x.count(j) for j in range(q))
    rows, drifts, brackets = {}, {}, {}
    total_rates, label_balance = [], True
    for i, x in enumerate(states):
        row, bval, gval = {}, F(0), F(0)
        for slot in range(2):
            for new in range(q):
                if new == x[slot]:
                    continue
                y = list(x)
                old = y[slot]
                y[slot] = new
                y = tuple(y)
                j = index[y]
                rate = F(amp[j], amp[i])
                key = counts(y)
                row[key] = row.get(key, F(0)) + rate
                dm = F(sig[new]-sig[old], 2)
                bval += rate * dm
                gval += rate * dm*dm
                label_balance &= F(amp[i]**2) * rate == F(amp[j]**2) / rate
        n = counts(x)
        if n in rows:
            assert rows[n] == row and drifts[n] == bval and brackets[n] == gval
        rows[n], drifts[n], brackets[n] = row, bval, gval
        total_rates.append(sum(row.values()))
    ck("q12_N2_exact_label_balance_and_144_to_78_count_lumping", label_balance and len(rows) == 78)
    correct = True
    for n, row in rows.items():
        bval, bbracket = F(0), F(0)
        for old in range(q):
            if not n[old]:
                continue
            cavity = [n[j] - (j == old) for j in range(q)]
            for new in range(q):
                if new == old:
                    continue
                r = F(2) ** (cavity[new]-cavity[old])
                nn = list(n)
                nn[old] -= 1
                nn[new] += 1
                correct &= row[tuple(nn)] == n[old]*r
                bval += F(n[old], 2)*r*(sig[new]-sig[old])
                bbracket += F(n[old], 2)*r*(sig[new]-sig[old])**2
        correct &= bval == drifts[n] and bbracket == 2*brackets[n]
    ck("q12_N2_actual_cavity_drift_and_bracket_sums", correct)
    same = (2,) + (0,)*11
    ck("small_N_noise_not_silently_set_to_2q", drifts[same] == -6 and 2*brackets[same] == 12 and max(total_rates) == 24)
    n = (1, 1) + (0,)*10
    removed_b, removed_g = F(0), F(0)
    for nn, rate in rows[n].items():
        if nn[0] > 1:
            dm = sum(sig[i]*(nn[i]-n[i]) for i in range(q))/F(2)
            removed_b += rate * dm
            removed_g += rate * dm*dm
    ck("reflection_exact_outgoing_drift_sign", drifts[n] == 0 and removed_b == 2 and drifts[n]-removed_b == -2)
    ck("reflection_positive_quadratic_term_retained", removed_g == 2 and 2*brackets[n]-2*removed_g == 24)
    ceil = lambda x, d: (x+d-1)//d
    resource = []
    for root in (12, 24, 120):
        nn = root**4
        space, clock = root, root**2
        dz = F(2*space, nn)
        ck(f"critical_bracket_clock_cancellation_root{root}", F(clock*space*space, nn) == 1 and dz == F(2, root**3))
        ck(f"third_and_fourth_jump_scale_root{root}", clock*nn*dz**3 == F(8, root**3) and dz**2 == F(4, root**6))
        c = 11*ceil(nn**3, 11)
        rcap, ss = F(3), F(2)
        amax = 1+rcap*nn
        budget = (F(7, 2)*amax*amax*clock*ss + 2*amax)/c
        uniform = F(7, 2)*(1+rcap)**2*ss + 2*(1+rcap)
        ck(f"finite_source_N3_horizon_budget_root{root}", nn**3 <= c < nn**3+11 and budget*root**2 <= uniform)
        resource.append({"N": nn, "c": c, "L": c//11, "fixture_R": str(rcap), "fixture_S": str(ss), "budget": str(budget), "scope": "symbolic budget normalization only, not native finite-N confidence"})
    # Scaled time-warp geometry, with large state values agreeing exactly.
    grid = [F(1, 16), F(1, 8)]
    actual = [F(1, 20), F(2, 17)]
    knots = [(F(0), F(0)), *zip(grid, actual), (F(1), F(1))]
    scaled = [(x/4, y/4) for x, y in knots]
    times = sorted({x for x, y in scaled} | {(scaled[i][0]+scaled[i+1][0])/2 for i in range(3)})
    ck("scaled_time_warp_no_amplitude_error_multiplier", max(abs(x-y) for x, y in scaled) <= F(1, 64) and all(at(t, [x/4 for x in grid], [0, 12, -4]) == at(warp(t, scaled), [x/4 for x in actual], [0, 12, -4]) for t in times))
    p = [[F(2, 3), F(1, 3)], [F(1, 4), F(3, 4)]]
    one = block_law([F(3, 4), F(1, 4)], p, 3)
    two = block_law([F(1, 2)]*2, p, 3)
    ck("whole_future_initial_TV_paid_once", sum(abs(one[w]-two[w]) for w in one)/2 == F(1, 4))
    # Original c=0 source, actual q3/N2 Gibbs model: not a spectator limit.
    s3 = list(product(range(3), repeat=2))
    cls = [x == y for x, y in s3]
    pi3 = [F(2, 9) if v else F(1, 18) for v in cls]
    z2 = [F(2, 15) if v else F(1, 270) for v in cls]
    uv = F(1, 45)
    p0 = [[F(0)]*9 for _ in range(9)]
    factors = set()
    for i, x in enumerate(s3):
        for j, y in enumerate(s3):
            if sum(a != b for a, b in zip(x, y)) > 1:
                continue
            zz = z2[i] if i == j else uv if cls[i] != cls[j] else F(1, 270)
            p0[i][j] = zz/pi3[i]
            if i != j:
                qrate = F(2 if cls[j] else 1, 2 if cls[i] else 1)
                factors.add(5*p0[i][j]/qrate)
    ck("unaugmented_source_actual_positive_scaling", all(sum(row) == 1 and min(row) >= 0 for row in p0) and all(pi3[i]*p0[i][j] == pi3[j]*p0[j][i] for i, j in product(range(9), repeat=2)))
    ck("unaugmented_source_not_equal_to_SQRT_after_degree_clock", factors == {F(1), F(1, 3)} and 5*z2[0]/pi3[0] == 3 and 5*z2[1]/pi3[1] == F(1, 3))
    lower, degree = F(1, 3), 5
    eps = lower*lower/(degree*degree)
    transformed = [[(p0[i][j]-(eps if i == j else 0))/(1-eps) for j in range(9)] for i in range(9)]
    ck("unaugmented_IF_lower_bound_pays_polynomial_self_loops", min(p0[i][i] for i in range(9)) >= eps and all(sum(row) == 1 and min(row) >= 0 for row in transformed))
    sos = True
    for f in ([F(i-4) for i in range(9)], [F(i*i % 7, 5) for i in range(9)]):
        pf = matvec(p0, f)
        lhs = sum(pi3[i]*f[i]*(pf[i]+(1-2*eps)*f[i]) for i in range(9))
        rhs = (1-eps)*sum(pi3[i]*transformed[i][j]*(f[i]+f[j])**2 for i, j in product(range(9), repeat=2))/2
        sos &= lhs == rhs and rhs >= 0
    ck("unaugmented_negative_spectrum_sum_of_squares_identity", sos)
    mu3, weights = [F(1, 9)]*9, [F(16, 9) if v else F(4, 9) for v in cls]
    law = mu3[:]
    dom = rowvec(weights, p0) == weights
    for _ in range(4):
        law = rowvec(law, p0)
        dom &= all(law[i] <= weights[i] for i in range(9))
    ck("unaugmented_source_weight_domination_in_intrinsic_steps", dom and all(sum(row[j] for j in range(9) if j != i) <= 1 for i, row in enumerate(p0)))
    geom = json.loads((ROOT / "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json").read_text())
    ck("actual_MT_spectral_and_threshold_inputs_bound", sha(ROOT/geom["spectral_path"]) == geom["spectral_sha256"] and F(geom["b"]) == F(1, 2*10**10) and F(geom["epsilon"]) == F(1, 100))
    return checks, {"nonzero_fast_quadratic_coefficient": str(at_origin), "q12_N2_drift_at_double_even": str(drifts[same]), "q12_N2_NGamma_at_double_even": str(2*brackets[same]), "unaugmented_normalized_edge_products": sorted(map(str, factors)), "scaling_budget_controls": resource}


def check_links():
    for p in HERE.glob("*.md"):
        lines, fence = [], False
        for line in p.read_text().splitlines():
            assert line == line.rstrip(), p.name
            if line.lstrip().startswith("```"):
                fence = not fence
            elif not fence and not line.startswith(("    ", "\t")):
                lines.append(line)
        text = re.sub(r"`[^`\n]*`", "", "\n".join(lines))
        for name in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", text):
            name = name.strip("<>").split("#")[0]
            if name and not name.startswith(("http:", "https:")):
                assert (p.parent/name).exists(), (p.name, name)


def duplicate_search():
    names = sorted({p.relative_to(ROOT).as_posix() for d in PROTECTED for p in (ROOT/d).glob("*.md")} | {n for n in INPUTS if n.endswith(".md") and n != "AGENTS.md"})
    pattern = r"SQRT|square.root|sinh|smooth.{0,40}fast|unaugmented|intrinsic.{0,35}scal|first.{0,20}fast.{0,25}deriv"
    r = subprocess.run(["rg", "-n", "-i", "--max-columns", "360", "--", pattern, *names], cwd=ROOT, capture_output=True, text=True, timeout=60)
    assert r.returncode in (0, 1), r.stderr
    return {"selected_files": names, "files": len(names), "pattern": pattern, "matches": sorted(r.stdout.splitlines()), "scope": "bounded predecessor search, not exhaustive novelty/absence", "reuse": "native UF/OBS/KR methods and paid SL/SR transfer; no claim of a new generic diffusion method", "new_dependency": "actual SQRT process and N3 finite source; distinct full global scaling of unaugmented hypothesis"}


def main():
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    inputs = {n: sha(ROOT/n) for n in INPUTS}
    protected = {d: manifest(ROOT/d) for d in PROTECTED}
    checks, fixtures = controls()
    duplicate = duplicate_search()
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    assert protected == {d: manifest(ROOT/d) for d in PROTECTED}
    assert inputs == {n: sha(ROOT/n) for n in INPUTS}
    result = {"date": "2026-10-20", "status": "PASS_SCOPED_SQ_CONTROLS", "bounded_exact_controls": len(checks), "checks": checks, "fixtures": fixtures, "fixture_scope": "q12 rational rank-five symmetry matrix and A=I q12/N2 rational Gibbs jumps, not exact A7/MT domain; q3 c0 source control", "arithmetic": "Fraction only, no binary64", "input_hashes": inputs, "protected_packages": protected, "protected_unchanged": True, "AGENTS_full_sha256": inputs["AGENTS.md"], "old_scientific_writers_executed": [], "independent_review": False, "internet_acquisition": False, "strict_F01": "OPEN_SOURCE", "broad_goal": "ACTIVE_NOT_COMPLETE"}
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result, indent=2)+"\n")
        (HERE/"DUPLICATE_SEARCH.json").write_text(json.dumps(duplicate, indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text("".join(f"{h}  {n}\n" for n, h in inputs.items()))
        check_links()
        files = sorted(p for p in HERE.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        manifest(HERE)
        assert json.loads((HERE/"RESULTS.json").read_text()) == result
        assert json.loads((HERE/"DUPLICATE_SEARCH.json").read_text()) == duplicate
        assert (HERE/"INPUTS.sha256").read_text() == "".join(f"{h}  {n}\n" for n, h in inputs.items())
        check_links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "bounded_exact_controls": len(checks), "bound_inputs": len(inputs), "protected_packages": len(protected), "protected_entries": sum(v["entries"] for v in protected.values()), "duplicate_search_files": duplicate["files"], "AGENTS_or_source_modified": False}, indent=2))


if __name__ == "__main__":
    main()
