#!/usr/bin/env python3
"""Exact finite controls, not an asymptotic MT proof or a physical test."""
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
    "AGENTS.md",
    "research_tasks/00_EXECUTION_RULES.md",
    "fin_spectator_limit_research_20261020/NEXT_FINITE_N_ORDER.md",
    "fin_spectator_limit_research_20261020/ACTIVITY_PROVIDER_CONTRACT.json",
    "fin_spectator_limit_research_20261020/ACTIVITY_AND_COUNT_CONSEQUENCE.md",
    "fin_spectator_limit_research_20261020/SCALING_AND_GENERATOR_PROOF.md",
    "fin_spectator_limit_research_20261020/PATH_LIMIT_AND_TV_SCOPE.md",
    "fin_spectator_composition_research_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_reference_transport_research_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md",
    "fin_critical_metastability_review_20261015/INTAKE_20261015.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
    "fin_critical_metastability_campaign_20261014/MT-001_COERCIVITY_CONNECTIVITY.md",
    "fin_critical_metastability_campaign_20261014/MT-002_IID_TRAPPING.md",
    "fin_critical_metastability_campaign_20261014/MT-003_CONDITIONAL_GIBBS.md",
    "fin_critical_metastability_campaign_20261014/MT-004_REFLECTED_MIXING_WINDOW.md",
    "fin_uniform_critical_campaign_20261014/UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md",
    "fin_uniform_critical_review_20261014/INTAKE_20261014.md",
    "fin_kinetic_invariants_campaign_20261014/KI-002/REVERSIBILITY_AND_MOBILITY.md",
    "fin_rank7_mathphysics_next/proofs/MP7-031_finite_N_variational_model.md",
    "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py",
]
PROTECTED = [
    "fin_spectator_limit_research_20261020",
    "fin_spectator_composition_research_20261020",
    "fin_reference_transport_research_20261020",
    "fin_hyp_history_entropy_20261020",
    "fin_hyp_history_entropy_review_20261020",
    "fin_critical_metastability_review_20261015",
    "fin_critical_metastability_campaign_20261014",
    "fin_uniform_critical_campaign_20261014",
    "fin_uniform_critical_review_20261014",
    "fin_kinetic_invariants_campaign_20261014",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def covered_manifest(folder):
    seen = set()
    for line in (folder / "MANIFEST.sha256").read_text().splitlines():
        h, name = line.split("  ", 1)
        p = Path(name)
        assert name not in seen and not p.is_absolute() and ".." not in p.parts
        assert sha(folder / name) == h, (folder.name, name)
        seen.add(name)
    actual = {
        p.relative_to(folder).as_posix()
        for p in folder.rglob("*")
        if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"
    }
    actual |= {n for n in seen if "__pycache__" in Path(n).parts or n.endswith(".pyc")}
    assert actual == seen | {"MANIFEST.sha256"}, (folder.name, actual - seen, seen - actual)
    return {"entries": len(seen), "sha256": sha(folder / "MANIFEST.sha256"), "coverage": "complete_noncache"}


def matvec(mu, p):
    return [sum((mu[i] * p[i][j] for i in range(len(mu))), F(0)) for j in range(len(mu))]


def generator(off):
    q = [row[:] for row in off]
    for i in range(len(q)):
        assert q[i][i] == 0
        q[i][i] = -sum(q[i])
    return q


def stochastic(p):
    return all(min(row) >= 0 and sum(row) == 1 for row in p)


def balance(q, pi):
    return all(pi[i] * q[i][j] == pi[j] * q[j][i] for i, j in product(range(len(pi)), repeat=2))


def reflected(q, keep, discrete=False):
    out = [[q[i][j] if i != j else F(0) for j in keep] for i in keep]
    for a in range(len(keep)):
        out[a][a] = (F(1) if discrete else F(0)) - sum(out[a])
    return out


def energy(q, pi, f):
    return sum((pi[i] * q[i][j] * (f[i] - f[j]) ** 2 for i, j in product(range(len(pi)), repeat=2)), F(0)) / 2


def count_quotient(states, matrix):
    q = max(max(s) for s in states) + 1
    counts = sorted({tuple(s.count(a) for a in range(q)) for s in states})
    groups = [[i for i, s in enumerate(states) if tuple(s.count(a) for a in range(q)) == n] for n in counts]
    rows = []
    for group in groups:
        candidate = [[sum((matrix[i][j] for j in target), F(0)) for target in groups] for i in group]
        assert all(row == candidate[0] for row in candidate)
        rows.append(candidate[0])
    return counts, groups, rows


def read_links():
    for p in HERE.glob("*.md"):
        prose, fenced = [], False
        for line in p.read_text().splitlines():
            assert line == line.rstrip(), p.name
            if line.lstrip().startswith("```"):
                fenced = not fenced
            elif not fenced and not line.startswith(("    ", "\t")):
                prose.append(line)
        text = re.sub(r"`[^`\n]*`", "", "\n".join(prose))
        for name in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", text):
            name = name.strip("<>").split("#")[0]
            if name and not name.startswith(("https:", "http:")):
                assert (p.parent / name).exists(), (p.name, name)


def duplicate_search():
    names = []
    for directory in PROTECTED:
        names.extend(p.relative_to(ROOT).as_posix() for p in (ROOT / directory).glob("*.md"))
    names.extend(INPUTS)
    names = sorted({n for n in names if n.endswith(".md") and n != "AGENTS.md"})
    pattern = r"SQRT|square.root|Dirichlet.{0,60}compar|spectator.{0,55}(mix|local|trap)|rate.{0,40}compar"
    result = subprocess.run(
        ["rg", "-n", "-i", "--max-columns", "350", "--", pattern, *names],
        cwd=ROOT, capture_output=True, text=True, timeout=60, check=False,
    )
    assert result.returncode in (0, 1), result.stderr
    return {
        "searched_files": names, "files": len(names), "pattern": pattern,
        "matching_lines": sorted(result.stdout.splitlines()),
        "scope": "Bounded local predecessor search, not exhaustive global absence or novelty proof",
        "reused": "MT method and generic form comparison; SQRT shape and source scaling already SL inputs",
        "new_dependency": "Actual SQRT local-window application; direct finite discrete source window with laziness",
    }


def controls():
    names = []

    def ck(name, value):
        assert name not in names
        assert value, name
        names.append(name)

    # Exact actual Gibbs rates in q3/N2, A=I, g=4 log(2). This is NOT A7/MT.
    states = list(product(range(3), repeat=2))
    amp = [2 if x == y else 1 for x, y in states]
    pi = [F(a * a, 18) for a in amp]
    mu = [F(1, 9)] * 9
    w = [F(4 * a * a, 9) for a in amp]
    qs, bs = [[F(0)] * 9 for _ in states], [[F(0)] * 9 for _ in states]
    for i, x in enumerate(states):
        for j, y in enumerate(states):
            changes = [a for a in range(2) if x[a] != y[a]]
            if len(changes) == 1:
                qs[i][j] = F(amp[j], amp[i])
                a = changes[0]
                bs[i][j] = F(4 if y[a] == x[1-a] else 1, 6)
    qs, bs = generator(qs), generator(bs)
    ck("q3_labelled_Gibbs_balance", balance(qs, pi) and balance(bs, pi) and sum(pi) == 1)
    ck("q3_labelled_true_SQRT_activity", all(qs[i][j] * qs[j][i] == 1 for i, j in product(range(9), repeat=2) if i != j and qs[i][j]))
    ratios = {qs[i][j] / bs[i][j] for i, j in product(range(9), repeat=2) if i != j and qs[i][j]}
    ck("q3_exact_field_ratio_not_common_clock", ratios == {F(3), F(6)})
    ck("q3_uniform_comparison_global_bound", min(ratios) >= F(3, 256) and max(ratios) <= 768)
    ck("q3_actual_exit_cap", max(-qs[i][i] for i in range(9)) == 6 and 6 <= 2 * 32)
    ck("q3_unnormalized_weight_dominates_iid", all(mu[i] <= w[i] for i in range(9)))
    ck("normalized_Gibbs_not_dominating_iid_control", any(mu[i] > pi[i] for i in range(9)))
    types, groups, q = count_quotient(states, qs)
    _, _, hb = count_quotient(states, bs)
    pc = [sum((pi[i] for i in group), F(0)) for group in groups]
    wc = [sum((w[i] for i in group), F(0)) for group in groups]
    mc = [sum((mu[i] for i in group), F(0)) for group in groups]
    ck("q3_count_lumping_both_laws", len(types) == 6 and balance(q, pc) and balance(hb, pc))
    cap = True
    multiplicity = True
    for i, n in enumerate(types):
        for j, m in enumerate(types):
            if i == j or not q[i][j]:
                continue
            a = next(a for a in range(3) if n[a] == m[a] + 1)
            b = next(a for a in range(3) if n[a] + 1 == m[a])
            cavity = [n[t] - (t == a) for t in range(3)]
            multiplicity &= q[i][j] == n[a] * F(2) ** (cavity[b] - cavity[a])
            multiplicity &= q[i][j] * q[j][i] == n[a] * (n[b] + 1)
            cap &= q[i][j] <= n[a] * 16
    ck("q3_exact_cavity_departure_reverse_multiplicity", multiplicity and cap)
    ck("q3_full_W_generator_invariance", matvec(wc, q) == [0] * 6 and matvec(wc, hb) == [0] * 6)
    keep = [i for i, n in enumerate(types) if n[2] <= 1]
    mass = sum(pc[i] for i in keep)
    pd = [pc[i] / mass for i in keep]
    qd, hd = reflected(q, keep), reflected(hb, keep)
    ck("q3_reflected_invariant_measure", len(keep) == 5 and balance(qd, pd) and balance(hd, pd))
    probes = [[F(0), F(1), F(-2), F(3), F(-1)], [F(i * i, 7) for i in range(5)]]
    ck("q3_actual_weighted_Dirichlet_comparison", all(3 * energy(hd, pd, f) <= energy(qd, pd, f) <= 6 * energy(hd, pd, f) for f in probes))
    f = probes[0]
    gamma = [sum((qd[i][j] * (f[j] - f[i]) ** 2 for j in range(5)), F(0)) for i in range(5)]
    ck("q3_Dirichlet_bracket_half_factor", 2 * energy(qd, pd, f) == sum(pd[i] * gamma[i] for i in range(5)))
    initial_in = sum(mc[i] for i in keep)
    nu = [mc[i] / initial_in for i in keep]
    density = [nu[i] / pd[i] for i in range(5)]
    bound = sum(wc[i] for i in keep) / initial_in
    norm2 = sum(pd[i] * (density[i] - 1) ** 2 for i in range(5))
    ck("q3_actual_iid_initial_density_bound", max(density) <= bound and norm2 <= bound - 1)
    unif = [[F(i == j) + q[i][j] / 8 for j in range(6)] for i in range(6)]
    ck("q3_uniformization_is_only_analysis_kernel", stochastic(unif) and matvec(wc, unif) == wc)
    law = mc[:]
    dominance = True
    for _ in range(4):
        law = matvec(law, unif)
        dominance &= all(law[i] <= wc[i] for i in range(6))
    ck("q3_CTMC_uniformization_weight_domination", dominance)
    summaries = []
    for c in (16, 32):
        # Closed positive scaling: u/v=2(c+3)/(c+1), pi_same/pi_diff=4.
        ratio = F(2 * (c + 3), c + 1)
        v2 = F(c + 1, 18 * (c + 3) * (c + 5))
        u2, uv = ratio * ratio * v2, ratio * v2
        source = [[F(0)] * 9 for _ in states]
        for i in range(9):
            for j in range(9):
                if i == j or qs[i][j]:
                    zprod = u2 if amp[i] == amp[j] == 2 else v2 if amp[i] == amp[j] == 1 else uv
                    source[i][j] = ((c + 1) if i == j else 1) * zprod / pi[i]
        ck(f"finite_source_exact_scaling_and_balance_c{c}", stochastic(source) and balance(source, pi))
        _, _, source_count = count_quotient(states, source)
        ck(f"finite_source_exact_count_projection_c{c}", stochastic(source_count) and balance(source_count, pc))
        ck(f"finite_source_laziness_c{c}", min(source_count[i][i] for i in range(6)) >= F(1, 2))
        r_same2 = F(c * (c + 3), (c + 1) * (c + 5))
        r_diff2 = F(c * (c + 1), (c + 3) * (c + 5))
        ck(f"finite_source_true_scaling_lower_bound_c{c}", (1 - F(7, c)) ** 2 <= min(r_same2, r_diff2) <= max(r_same2, r_diff2) <= 1)
        exact = True
        source_ratios = []
        for i in range(9):
            for j in range(9):
                if i == j or not qs[i][j]:
                    continue
                factor = F(c, c + 5) if amp[i] != amp[j] else r_diff2
                exact &= c * source[i][j] == qs[i][j] * factor
                source_ratios.append(c * source[i][j] / qs[i][j])
        ck(f"finite_source_exact_corrected_rates_c{c}", exact and F(1, 4) <= min(source_ratios) <= max(source_ratios) <= 1)
        ck(f"finite_source_not_equal_to_SQRT_c{c}", max(source_ratios) < 1 and len(set(source_ratios)) > 1)
        source_d = reflected(source_count, keep, discrete=True)
        ck(f"finite_source_reflected_laziness_and_pi_c{c}", stochastic(source_d) and balance(source_d, pd) and min(source_d[i][i] for i in range(5)) >= F(1, 2))
        source_g = [[c * (source_d[i][j] - F(i == j)) for j in range(5)] for i in range(5)]
        ck(f"finite_source_reflected_form_comparison_c{c}", all(F(1, 4) * energy(qd, pd, f) <= energy(source_g, pd, f) <= energy(qd, pd, f) for f in probes))
        law = mc[:]
        dom = matvec(wc, source_count) == wc
        for _ in range(4):
            law = matvec(law, source_count)
            dom &= all(law[i] <= wc[i] for i in range(6))
        ck(f"finite_source_weight_domination_c{c}", dom)
        # Exact absorbing computation of the actual exit event; initial outside included.
        alive = [mc[i] for i in keep]
        expected_exits = F(0)
        min_out = sum(mc[i] for i in range(6) if i not in keep)
        for _ in range(4):
            expected_exits += sum(alive[a] * sum(source_count[i][j] for j in range(6) if j not in keep) for a, i in enumerate(keep))
            alive = [sum((alive[a] * source_count[i][j] for a, i in enumerate(keep)), F(0)) for j in keep]
        exit_probability = 1 - sum(alive)
        ck(f"finite_source_exit_compensator_initial_once_c{c}", exit_probability == min_out + expected_exits and exit_probability > min_out)
        layer = [i for i in keep if any(q[i][j] for j in range(6) if j not in keep)]
        ck(f"finite_source_boundary_flux_bound_c{c}", exit_probability <= min_out + F(4 * 6, c) * sum(wc[i] for i in layer))
        summaries.append({"c": c, "L": c // 2, "scaling_u_over_v": str(ratio), "r_same_squared": str(r_same2), "r_diff_squared": str(r_diff2), "actual_four_step_exit": str(exit_probability)})
    flip = [[F(0), F(1)], [F(1), F(0)]]
    periodic = [F(1), F(0)]
    tv = []
    for _ in range(4):
        periodic = matvec(periodic, flip)
        tv.append(sum(abs(x - F(1, 2)) for x in periodic) / 2)
    ck("continuous_gap_without_discrete_laziness_countercontrol", tv == [F(1, 2)] * 4)
    # Exact q12 zero-g labelled quotient, with the true 12-vs-1 rate convention.
    states12 = list(product(range(12), repeat=2))
    count12 = lambda x: tuple(x.count(a) for a in range(12))
    rows = {}
    for x in states12:
        record = {}
        for slot in range(2):
            for new in range(12):
                if new == x[slot]:
                    continue
                y = list(x)
                y[slot] = new
                n = count12(tuple(y))
                record[n] = record.get(n, 0) + 1
        key = count12(x)
        if key in rows:
            assert rows[key] == record
        rows[key] = record
    ck("q12_g0_actual_label_count_quotient_144_to_78", len(rows) == 78 and len(states12) == 144)
    multiplicities = True
    clock = True
    finite = True
    c = 66  # q12 hidden L=6; c >= 2 A_N=46 in this exact g=0 fixture.
    for n, record in rows.items():
        multiplicities &= sum(record.values()) == 22
        for m, rate in record.items():
            a = next(a for a in range(12) if n[a] == m[a] + 1)
            b = next(a for a in range(12) if n[a] + 1 == m[a])
            multiplicities &= rate == n[a]
            clock &= F(rate) / F(rate, 12) == 12
            finite &= F(c * rate, c + 23) >= F(rate, 4)
            finite &= F(c * rate, c + 23) <= rate
    ck("q12_g0_departure_multiplicity_and_true_cap", multiplicities)
    ck("q12_native_clock_not_silently_rescaled", clock)
    ck("q12_uniform_finite_source_correction", finite and F(c + 1, c + 23) >= F(1, 2))
    # Exact integer choices; no floating roots/ceil on enormous integers.
    r = F(32)
    floor_ok, resource_ok = True, True
    ceil = lambda a, b: (a + b - 1) // b
    for n in (12, 120, 12000, 12 * 10**30):
        ell = ceil(2 * (1 + int(r) * n), 11)
        c = 11 * ell
        floor_ok &= c >= 2 * (1 + r * n) and c < 2 * (1 + r * n) + 11
        resource_ok &= ell <= 2 * (1 + r * n) / 11 + 1
    ck("integer_spectator_count_pays_global_cap", floor_ok and resource_ok)
    ck("next_N3_sequence_exact_divisibility", all(11 * ceil(n ** 3, 11) >= n ** 3 and 11 * ceil(n ** 3, 11) < n ** 3 + 11 for n in (12, 120, 12000)))
    geom = json.loads((ROOT / "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json").read_text())
    assert sha(ROOT / geom["spectral_path"]) == geom["spectral_sha256"]
    spectra = json.loads((ROOT / geom["spectral_path"]).read_text())["exact"]["laplacian_intervals"]
    boxes = {k: tuple(map(F, spectra[k])) for k in (3, 4, 5, 6)}
    ck("actual_spectral_boxes_bind_MT_gap_premise", F(2) < boxes[6][0] <= boxes[6][1] < F(12, 5) and 1 - max(boxes[k][1] for k in (3, 4, 5)) / boxes[6][0] > F(93, 5000))
    ck("actual_MT_thresholds_not_enlarged", F(geom["h"]) == F(1, 10**10) and F(geom["b"]) == F(1, 2 * 10**10) and F(geom["epsilon"]) == F(1, 100) and F(geom["delta"]) == F(1, 10000))
    ck("actual_MT_positive_fast_margin_is_inherited", F(geom["fast_hessian_floor"]) > F(8, 100) and F(geom["gradient_fast_margin"]) > F(793, 10000))
    return names, summaries


def main():
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    before = {n: covered_manifest(ROOT / n) for n in PROTECTED}
    inputs = {n: sha(ROOT / n) for n in INPUTS}
    checks, fixtures = controls()
    duplicate = duplicate_search()
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    after = {n: covered_manifest(ROOT / n) for n in PROTECTED}
    assert before == after
    assert inputs == {n: sha(ROOT / n) for n in INPUTS}
    result = {
        "date": "2026-10-20", "status": "PASS_SCOPED_SR_CONTROLS",
        "scientific_control_count": len(checks), "scientific_controls": checks,
        "scope": "Exact small rational controls and integrity; not an all-N stochastic proof or independent review",
        "q3_fixture": {"q": 3, "N": 2, "A": "identity", "g": "4 log(2)", "count_states": 6, "reflected_states": 5, "exact_A_N": 7, "finite_sources": fixtures, "condition_scope": "Exact A_N=7 pays c>=2A_N in these controls; the coarse universal field cap is deliberately looser"},
        "q12_fixture": {"q": 12, "N": 2, "g": 0, "labelled_states": 144, "count_states": 78, "scope": "multiplicity and clock only, not MT domain"},
        "protected_packages": before, "protected_unchanged_after_checks": before == after,
        "input_hashes": inputs, "AGENTS_full_sha256": inputs["AGENTS.md"],
        "scripts_replayed_from_old_sources": [], "scientific_arithmetic": "Fraction only; no binary64",
        "new_agent_or_independent_review": False, "internet_acquisition": False,
        "strict_F01": "OPEN_SOURCE", "broad_goal": "ACTIVE_NOT_COMPLETE",
    }
    if args.record:
        (HERE / "RESULTS.json").write_text(json.dumps(result, indent=2) + "\n")
        (HERE / "DUPLICATE_SEARCH.json").write_text(json.dumps(duplicate, indent=2) + "\n")
        (HERE / "INPUTS.sha256").write_text("".join(f"{h}  {n}\n" for n, h in inputs.items()))
        read_links()
        files = sorted(p for p in HERE.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE / "MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        covered_manifest(HERE)
        assert json.loads((HERE / "RESULTS.json").read_text()) == result
        assert json.loads((HERE / "DUPLICATE_SEARCH.json").read_text()) == duplicate
        assert (HERE / "INPUTS.sha256").read_text() == "".join(f"{h}  {n}\n" for n, h in inputs.items())
        read_links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents, "Use a /tmp replay output, not source evidence"
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "bounded_exact_controls": len(checks), "source_inputs": len(inputs), "protected_packages": len(before), "duplicate_files": duplicate["files"], "source_or_AGENTS_modified": False}, indent=2))


if __name__ == "__main__":
    main()
