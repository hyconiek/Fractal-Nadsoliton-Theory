#!/usr/bin/env python3
"""Exact bounded AD controls. Not an independent or mechanized path proof."""
import argparse
import ast
import hashlib
import json
import re
import subprocess
from fractions import Fraction as F
from itertools import product
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "fin_intrinsic_scaling_research_20261020/NEXT_GLOBAL_ACTIVITY_DEFECT_ORDER.md",
    "fin_intrinsic_scaling_research_20261020/GLOBAL_SCALING_THEOREM.md",
    "fin_intrinsic_scaling_research_20261020/ACTUAL_LOCAL_WINDOW.md",
    "fin_intrinsic_scaling_research_20261020/DISCRETE_CRITICAL_GATE.md",
    "fin_intrinsic_scaling_research_20261020/DUPLICATE_SEARCH.json",
    "fin_sqrt_critical_process_research_20261020/MICROSCOPIC_CHARACTERISTICS.md",
    "fin_sqrt_critical_process_research_20261020/STATIONARY_PATH_THEOREM.md",
    "fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md",
    "fin_hyp_history_entropy_review_20261020/PROOF_AUDIT_AND_DUAL.md",
    "fin_reference_transport_research_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
]
PROTECTED = [
    "fin_intrinsic_scaling_research_20261020",
    "fin_sqrt_critical_process_research_20261020",
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
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        h, name = line.split("  ", 1)
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name) == h, (folder.name, name)
        seen.add(name)
    actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
    actual |= {n for n in seen if "__pycache__" in Path(n).parts or n.endswith(".pyc")}
    assert actual == seen | {"MANIFEST.sha256"}, (folder.name, actual-seen, seen-actual)
    return {"entries": len(seen), "sha256": sha(folder/"MANIFEST.sha256"), "coverage": "complete_noncache"}


def log_box(x, terms=80):
    """Outward exact log enclosure using the positive atanh series."""
    assert x > 0
    if x == 1:
        return F(0), F(0)
    if x < 1:
        lo, hi = log_box(1/x, terms)
        return -hi, -lo
    z = (x-1)/(x+1)
    lower = 2*sum((z**(2*k+1)/F(2*k+1) for k in range(terms)), F(0))
    tail = 2*z**(2*terms+1)/F(2*terms+1)/(1-z*z)
    return lower, lower+tail


def times_interval(a, interval):
    lo, hi = interval
    return (a*lo, a*hi) if a >= 0 else (a*hi, a*lo)


def recorded_box(lo, hi, digits=30):
    """Keep exact scientific arithmetic; round only recorded endpoints outward."""
    scale = 10**digits
    lower = F((lo.numerator*scale)//lo.denominator, scale)
    upper = F(-((-hi.numerator*scale)//hi.denominator), scale)
    assert lower <= lo <= hi <= upper
    return {"lower": str(lower), "upper": str(upper), "endpoint_grid": f"1/10^{digits}", "rounding": "exact outward, output only"}


def controls():
    names = []

    def ck(name, ok):
        assert name not in names and ok, name
        names.append(name)

    jet_data = []
    # Formal first-order row algebra, with a nonzero transport term retained.
    for name, r, a0, kap, transport in (
        ("jet1", F(2), F(3, 2), F(5, 7), F(-4, 5)),
        ("jet2", F(3), F(2, 3), F(-1, 5), F(7, 9)),
    ):
        tau = r/(a0*a0)
        a1 = (1-(1+kap)*a0*a0-a0*transport)/(2*tau*a0)
        coefficient = 2*tau*a0*a1+(1+kap)*a0*a0+a0*transport
        wrong_kap = (1-a0*a0-a0*transport)/(2*tau*a0)
        wrong_transport = (1-(1+kap)*a0*a0)/(2*tau*a0)
        ck(name+"_first_corrector_exact_row_coefficient", tau*a0*a0 == r and coefficient == 1)
        ck(name+"_reject_omitted_cavity_correction", 2*tau*a0*wrong_kap+(1+kap)*a0*a0+a0*transport-1 == kap*a0*a0 != 0)
        ck(name+"_reject_omitted_transport_gradient", 2*tau*a0*wrong_transport+(1+kap)*a0*a0+a0*transport-1 == a0*transport != 0)
        jet_data.append({"name": name, "r": str(r), "tau": str(tau), "a0": str(a0), "a1": str(a1), "scope": "generic exact Taylor-jet algebra, not an A7 numerical rate"})
    # Actual field derivatives on a separate rational q12 circulant rank-five
    # symmetry fixture; no floating A7 substitution.
    q, g = 12, F(12)
    sig = [(-1)**j for j in range(q)]
    cos3 = [F(1), F(0), F(-1), F(0)]
    cos4 = [F(1), F(-1, 2), F(-1, 2)]
    a = [[F(sig[i]*sig[j], q)+F(1, 12)*cos3[(i-j) % 4]+F(1, 18)*cos4[(i-j) % 3] for j in range(q)] for i in range(q)]
    ck("rational_field_fixture_row_zero_and_top_parity", all(sum(row) == 0 for row in a) and [sum(a[i][j]*sig[j] for j in range(q)) for i in range(q)] == sig)
    tr = sum(a[i][i] for i in range(q))
    tr2 = sum(a[i][j]*a[j][i] for i, j in product(range(q), repeat=2))
    ck("rational_field_fixture_not_A7_scope", tr == F(8, 3) and tr2 == F(31, 18))
    r = F(q-1)
    ell = [[g*(a[i][i]-a[j][i])/2 for j in range(q)] for i in range(q)]
    kap = sum(ell[i][j]/q for i, j in product(range(q), repeat=2) if i != j)
    a1 = -kap/(2*r)
    ck("uniform_cavity_trace_correction_exact", kap == g*tr/2 == 16 and a1 == F(-8, 11))
    sum_a02, sum_ell2, sum_sec = F(0), F(0), F(0)
    derivatives = True
    for old, new in product(range(q), repeat=2):
        if old == new:
            continue
        direction = [F((j == new)-(j == old)) for j in range(q)]
        ad = [sum(a[i][j]*direction[j] for j in range(q)) for i in range(q)]
        tau2, kap1, tau1 = F(0), F(0), F(0)
        for i, j in product(range(q), repeat=2):
            if i == j:
                continue
            e1 = g*(ad[j]-ad[i])/2
            dc = direction[i]+e1/q
            tau1 += dc
            kap1 += dc*ell[i][j]
            tau2 += direction[i]*e1+e1*e1/(2*q)
        derivatives &= tau1 == kap1 == 0
        a02 = -tau2/(2*r)
        sum_a02 += a02/q
        sum_ell2 += ell[old][new]**2/(2*q)
        sum_sec += (ell[old][new]**2/2+a1*ell[old][new]+a02)/q
    ck("uniform_neighbour_first_derivatives_vanish_exactly", derivatives)
    ck("uniform_second_direction_sum_independent_trace_route", sum_a02 == -(-q*g*tr+g*g*tr2/2)/(2*r))
    ck("uniform_LOO_second_sum_independent_trace_route", sum_ell2 == g*g*(tr*tr+tr2)/(8*q))
    first = (kap+2*r*a1)/r
    second = (sum_sec+a1+a1*(1+kap+r*a1))/r
    ck("uniform_actual_row_residual_order_two_not_order_one", first == 0 and second != 0)
    ck("leading_only_trial_fails_order_two_requirement", kap/r == F(16, 11) > 0)
    ck("zero_gain_uniform_trial_exact_control", F(1)-(1+0)*F(1)**2 == 0 and (r+F(1, 37))/(r+F(1, 37)) == 1)
    jet_data.append({"name": "rational_q12_field", "rank": 5, "g": str(g), "trace": str(tr), "kappa": str(kap), "a1_uniform": str(a1), "corrected_residual_first_coefficient": str(first), "corrected_residual_second_coefficient": str(second), "scope": "formal true LOO Taylor coefficients of a separate rational field fixture, not actual A7"})
    # Exact source flows in q3/N2 with A=I, g=4log2, then a positive
    # auxiliary perturbation. Flow-localization algebra does not require
    # that this finite auxiliary was built by the asymptotic trial above.
    states = list(product(range(3), repeat=2))
    size = len(states)
    z = [F(6) if x == y else F(1) for x, y in states]
    support = [[sum(x != y for x, y in zip(states[i], states[j])) <= 1 for j in range(size)] for i in range(size)]
    w = [z[i]*sum(z[j] for j in range(size) if support[i][j]) for i in range(size)]
    total = sum(w)
    pi = [value/total for value in w]
    c = [[z[i]*z[j]/total if support[i][j] else F(0) for j in range(size)] for i in range(size)]
    t = [F(200+(i % 3)-1, 200) for i in range(size)]
    u = [[c[i][j]*(t[i]*t[j])**2 for j in range(size)] for i in range(size)]
    margins = [sum(row) for row in u]
    delta = max(abs(margins[i]/pi[i]-1) for i in range(size))
    lmax = max(max(abs(lo), abs(hi)) for lo, hi in (log_box(1/(x*x)) for x in t))
    ck("finite_exact_source_and_auxiliary_margins", sum(pi) == 1 and all(sum(c[i]) == pi[i] for i in range(size)) and 0 < delta < F(1, 2))
    p = [[c[i][j]/pi[i] for j in range(size)] for i in range(size)]
    papp = [[u[i][j]/pi[i] if i != j else F(0) for j in range(size)] for i in range(size)]
    for i in range(size):
        papp[i][i] = 1-sum(papp[i])
    ck("conservative_auxiliary_diagonal_charges_margin_error", all(papp[i][i] == u[i][i]/pi[i]-(margins[i]/pi[i]-1) > 0 for i in range(size)) and all(sum(row) == 1 for row in papp))
    ck("conservative_analysis_table_reversible_not_selected_source", all(pi[i]*papp[i][j] == pi[j]*papp[j][i] for i, j in product(range(size), repeat=2)) and papp != p)
    flow_data = []
    for name, inside in (("global", list(range(size))), ("local", [0, 1, 3, 4])):
        outside = [i for i in range(size) if i not in inside]
        mass = sum(pi[i] for i in inside)
        boundary = [i for i in inside if any(support[i][j] for j in outside)]
        bmass = sum(pi[i] for i in boundary)
        left = [F(0)]*size
        cross = [F(0)]*size
        interior = [F(0)]*size
        slo, shi = F(0), F(0)
        h2 = F(0)
        variation, weighted_square = F(0), F(0)
        for i in inside:
            left[i] = pi[i]-margins[i]
            cross[i] = sum(c[i][j]-u[i][j] for j in outside)
            weight = 1+F(i, 3)
            for j in inside:
                if not support[i][j]:
                    continue
                coef = (c[i][j]-u[i][j])/2
                interior[i] += coef
                interior[j] += coef
                lo, hi = times_interval(coef, log_box(c[i][j]/u[i][j]))
                slo += lo
                shi += hi
                h2 += c[i][j]*(1-t[i]*t[j])**2
                variation += weight*abs(c[i][j]-u[i][j])
                weighted_square += weight*weight*c[i][j]*(1+t[i]*t[j])**2
        ck(name+"_local_log_identity_exact_coefficients", all(interior[i] == left[i]-cross[i] for i in range(size)))
        if outside:
            ck("reject_reversed_boundary_sign", any(interior[i] != left[i]+cross[i] for i in range(size)))
        ck(name+"_scalar_log_Hellinger_inequality_outward", 0 < 2*h2 <= slo <= shi)
        upper = lmax*(delta*mass+(2+delta)*bmass)
        ck(name+"_localized_margin_plus_boundary_bound", shi <= upper and h2 <= upper/2)
        ck(name+"_weighted_Hellinger_Cauchy_control", variation*variation <= h2*weighted_square)
        moment = sum(pi[i]*(1+F(i, 3))**2 for i in inside)
        ck(name+"_weighted_row_moment_payment", weighted_square <= 6*moment)
        flow_data.append({"subset": name, "pi_mass": str(mass), "departure_mass": str(bmass), "margin_defect_sup": str(delta), "Hellinger_squared": str(h2), "recorded_log_form_box": recorded_box(slo, shi), "scope": "finite exact positive-flow identity, not an all-N A7 simulation"})
    # Boundary omission is genuinely false: complete q3,N1 uniform source,
    # raw auxiliary zapp=(1/6,1/6,5/3). Margins inside {0,1} are unchanged.
    cc = [[F(1, 9)]*3 for _ in range(3)]
    zz = [F(1, 6), F(1, 6), F(5, 3)]
    uu = [[x*y for y in zz] for x in zz]
    ck("boundary_countercontrol_zero_inside_margin_defect", sum(uu[0]) == sum(uu[1]) == F(1, 3))
    log2lo, log2hi = log_box(F(2))
    ck("boundary_countercontrol_positive_internal_log_form", log2lo/3 > 0 and log2lo/3 >= F(2, 9))
    ck("boundary_countercontrol_wrong_omission_rejected", sum(cc[i][j]-uu[i][j] for i in (0, 1) for j in (2,)) == F(-1, 3) and log2lo/3 > 0)
    # Exact centered-step variance identity for the conservative source/table.
    f = [F(x+y, 2) for x, y in states]
    h = F(1, 17)
    for name, matrix in (("exact_source", p), ("analysis_table", papp)):
        valid, strict = True, False
        for i in range(size):
            mean = sum(matrix[i][j]*(f[j]-f[i]) for j in range(size))
            second_raw = sum(matrix[i][j]*(f[j]-f[i])**2 for j in range(size))
            variance = sum(matrix[i][j]*(f[j]-f[i]-mean)**2 for j in range(size))
            beta, raw = mean/h, second_raw/h
            valid &= raw-h*beta*beta == variance/h >= 0
            strict |= mean != 0 and variance < second_raw
        ck(name+"_actual_discrete_centering_not_CTMC_Gamma", valid and strict)
    for root in (12, 24):
        n = root**4
        degree = 1+11*n
        step = F(1, degree*root*root)
        maxbeta = F(2*degree, root)
        ck(f"critical_centering_scale_root{root}", step*maxbeta == F(2, root**3))
        # Only a symbolic defect budget, not calibrated probabilities.
        ck(f"weighted_drift_prefactor_times_second_order_flow_root{root}", F(degree, root*n) <= F(12, root))
    geom = json.loads((ROOT/"fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json").read_text())
    ck("actual_spectral_and_true_boundary_inputs_unchanged", sha(ROOT/geom["spectral_path"]) == geom["spectral_sha256"] and F(geom["b"]) == F(1, 2*10**10))
    return names, jet_data, flow_data


def links():
    for p in HERE.glob("*.md"):
        prose, fence = [], False
        for line in p.read_text().splitlines():
            assert line == line.rstrip(), p.name
            if line.lstrip().startswith("```"):
                fence = not fence
            elif not fence and not line.startswith(("    ", "\t")):
                prose.append(line)
        text = re.sub(r"`[^`\n]*`", "", "\n".join(prose))
        for name in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", text):
            name = name.strip("<>").split("#")[0]
            if name and not name.startswith(("http:", "https:")):
                assert (p.parent/name).exists(), (p.name, name)


def duplicate_search():
    old = json.loads((ROOT/"fin_intrinsic_scaling_research_20261020/DUPLICATE_SEARCH.json").read_text())
    names = sorted(set(old["selected_files"]) | {p.relative_to(ROOT).as_posix() for d in PROTECTED for p in (ROOT/d).glob("*.md")})
    pat = r"Hellinger|unequal.{0,25}margin|scaling.{0,35}(residual|correct)|logarithmic.{0,25}flow|localiz.{0,25}entropy|cavity"
    matched = []
    for i in range(0, len(names), 90):
        r = subprocess.run(["rg", "-n", "-i", "--max-columns", "350", "--", pat, *names[i:i+90]], cwd=ROOT, capture_output=True, text=True, timeout=60)
        assert r.returncode in (0, 1), r.stderr
        matched.extend(r.stdout.splitlines())
    return {"selected_files": names, "files": len(names), "pattern": pat, "matching_lines": sorted(matched), "scope": "bounded predecessor/root search; no exhaustive absence or mathematical priority claim", "reuse": "known entropy/log monotonicity and martingale methods; actual HYP/RT/IS/MT/SQ chain", "new_dependency": "global first-corrected actual row trial and localized unequal-margin flow bound pay original-source critical defects"}


def main():
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    inputs = {n: sha(ROOT/n) for n in INPUTS}
    before = {d: manifest(ROOT/d) for d in PROTECTED}
    names, jets, flows = controls()
    duplicate = duplicate_search()
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    assert before == {d: manifest(ROOT/d) for d in PROTECTED}
    assert inputs == {n: sha(ROOT/n) for n in INPUTS}
    result = {"date": "2026-10-20", "status": "PASS_SCOPED_AD_CONTROLS", "bounded_exact_controls": len(names), "checks": names, "Taylor_controls": jets, "flow_controls": flows, "arithmetic": "Fraction plus exact outward atanh-series log bounds; no binary64", "scope": "bounded algebra and integrity, not mechanical all-N proof or independent review", "input_hashes": inputs, "protected_packages": before, "protected_unchanged": True, "AGENTS_full_sha256": inputs["AGENTS.md"], "old_scientific_writers_executed": [], "independent_review": False, "network_acquisition": False, "strict_F01": "OPEN_SOURCE", "broad_goal": "ACTIVE_NOT_COMPLETE"}
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result, indent=2)+"\n")
        (HERE/"DUPLICATE_SEARCH.json").write_text(json.dumps(duplicate, indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text("".join(f"{h}  {n}\n" for n, h in inputs.items()))
        links()
        files = sorted(p for p in HERE.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        manifest(HERE)
        assert json.loads((HERE/"RESULTS.json").read_text()) == result
        assert json.loads((HERE/"DUPLICATE_SEARCH.json").read_text()) == duplicate
        assert (HERE/"INPUTS.sha256").read_text() == "".join(f"{h}  {n}\n" for n, h in inputs.items())
        links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "bounded_exact_controls": len(names), "inputs": len(inputs), "protected_packages": len(before), "protected_entries": sum(v["entries"] for v in before.values()), "duplicate_search_files": duplicate["files"], "AGENTS_or_source_modified": False}, indent=2))


if __name__ == "__main__":
    main()
