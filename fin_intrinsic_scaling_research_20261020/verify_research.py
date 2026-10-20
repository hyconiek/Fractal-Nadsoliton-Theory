#!/usr/bin/env python3
"""Exact finite controls for IS. Not an all-N A7 path or formal proof."""
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
    "fin_sqrt_critical_process_research_20261020/NEXT_INTRINSIC_SCALING_ORDER.md",
    "fin_sqrt_critical_process_research_20261020/UNAUGMENTED_SOURCE_FRONTIER.md",
    "fin_sqrt_critical_process_research_20261020/STATIONARY_PATH_THEOREM.md",
    "fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md",
    "fin_hyp_history_entropy_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_hyp_history_entropy_review_20261020/PROOF_AUDIT_AND_DUAL.md",
    "fin_reference_transport_research_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_sqrt_local_regime_research_20261020/RATE_AND_FORM_COMPARISON.md",
    "fin_sqrt_local_regime_research_20261020/LOCAL_WINDOW_THEOREM.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
    "fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json",
]
PROTECTED = [
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
    return {"entries": len(seen), "manifest_sha256": sha(folder/"MANIFEST.sha256"), "coverage": "complete_noncache"}


def neighbors(x, q):
    out = [x]
    for slot in range(len(x)):
        for a in range(q):
            if a != x[slot]:
                y = list(x)
                y[slot] = a
                out.append(tuple(y))
    return out


def source(q, n, weight):
    states = list(product(range(q), repeat=n))
    index = {x: i for i, x in enumerate(states)}
    z = [F(weight(x)) for x in states]
    adj = [[index[y] for y in neighbors(x, q)] for x in states]
    assert all(len(set(row)) == 1+n*(q-1) for row in adj)
    w = [z[i]*sum(z[j] for j in row) for i, row in enumerate(adj)]
    total = sum(w)
    pi = [x/total for x in w]
    p = [[F(0)]*len(states) for _ in states]
    for i, row in enumerate(adj):
        for j in row:
            p[i][j] = z[i]*z[j]/w[i]
    return states, index, z, adj, w, total, pi, p


def rowvec(mu, p):
    return [sum((mu[i]*p[i][j] for i in range(len(mu))), F(0)) for j in range(len(mu))]


def energy(p, pi, f):
    return sum((pi[i]*p[i][j]*(f[i]-f[j])**2 for i, j in product(range(len(pi)), repeat=2)), F(0))/2


def checks():
    names, summaries = [], []

    def ck(name, good):
        assert name not in names and good, name
        names.append(name)

    cases = [
        ("nonexchangeable_binary", 2, 3, lambda x: 1+sum((i+1)*(a+1)**2 for i, a in enumerate(x))),
        ("exchangeable_ternary", 3, 3, lambda x: 1+sum(a*a for a in x)+3*(len(set(x)) == 1)),
        ("uniform", 3, 2, lambda x: 1),
        ("q3_Gibbs_N2", 3, 2, lambda x: 6 if x[0] == x[1] else 1),
    ]
    for name, q, n, weight in cases:
        states, index, z, adj, w, total, pi, p = source(q, n, weight)
        edges = [(i, j) for i in range(len(z)) for j in adj[i] if i != j]
        c = max(pi[i]/pi[j] for i, j in edges)
        m = max(z[i]/z[j] for i, j in edges)
        mbar = c+(c-1)*(n-1)*(q-1)
        d, dout = 1+n*(q-1), (n-1)*(q-1)
        bound = 1+n*(q-1)*mbar
        ck(name+"_source_exact_scaling_and_stationarity", all(sum(row) == 1 and min(row) >= 0 for row in p) and rowvec(pi, p) == pi)
        ck(name+"_all_oriented_neighbour_ratios", 1 <= m <= mbar and all(1/mbar <= z[i]/z[j] <= mbar for i, j in edges))
        i, j = next((i, j) for i, j in edges if z[i]/z[j] == m)
        x, y = states[i], states[j]
        slot = next(k for k in range(n) if x[k] != y[k])
        clique = []
        for a in range(q):
            xx = list(x)
            xx[slot] = a
            clique.append(index[tuple(xx)])
        sx = sum(z[k] for k in clique)
        hx = sum(z[k] for k in adj[i] if k not in clique)
        hy = sum(z[k] for k in adj[j] if k not in clique)
        ck(name+"_whole_clique_pairing_and_maximum_chain", len(clique) == q and hx >= hy/m and hy <= dout*m*z[j] and hy <= dout*sx and pi[i]/pi[j] == m*(sx+hx)/(sx+hy) and c >= pi[i]/pi[j] >= 1+(m-1)*sx/(sx+hy) >= 1+(m-1)/(1+dout))
        self_p = [z[i]*z[i]/w[i] for i in range(len(z))]
        ck(name+"_global_self_probability_lower_bound", all(p[i][i] == self_p[i] >= 1/bound for i in range(len(z))))
        ck(name+"_normalized_activity_lower_bound_squared", all(d*d*self_p[i]*self_p[j] >= (d/bound)**2 for i, j in edges))
        ck(name+"_polynomial_constant_arithmetic", bound <= q*(1+(q-1)*(c-1))*n*n and bound >= d >= 2)
        eps = 1/bound
        s = [[(p[i][j]-(eps if i == j else 0))/(1-eps) for j in range(len(z))] for i in range(len(z))]
        ck(name+"_negative_spectrum_stochastic_remainder", all(sum(row) == 1 and min(row) >= 0 for row in s) and rowvec(pi, s) == pi)
        f = [F(sum((k+1)*a for k, a in enumerate(x)), n) for x in states]
        lhs = sum(pi[i]*f[i]*(sum(p[i][j]*f[j] for j in range(len(z)))+(1-2*eps)*f[i]) for i in range(len(z)))
        rhs = (1-eps)*sum(pi[i]*s[i][j]*(f[i]+f[j])**2 for i, j in product(range(len(z)), repeat=2))/2
        ck(name+"_negative_spectrum_sum_of_squares", lhs == rhs >= 0)
        if name == "uniform":
            ck("uniform_exact_constant_scaling_control", c == m == mbar == 1 and bound == d and all(p[i][j] == F(j in adj[i], d) for i, j in product(range(len(z)), repeat=2)))
        if name.startswith("exchangeable") or name == "q3_Gibbs_N2":
            count = lambda x: tuple(x.count(a) for a in range(q))
            types = sorted({count(x) for x in states})
            groups = [[i for i, x in enumerate(states) if count(x) == nn] for nn in types]
            pm = [sum(pi[i] for i in group) for group in groups]
            pcm = []
            for group in groups:
                rows = [[sum(p[i][j] for j in target) for target in groups] for i in group]
                assert all(row == rows[0] for row in rows)
                pcm.append(rows[0])
            correct = True
            for ii, nn in enumerate(types):
                i = groups[ii][0]
                correct &= pcm[ii][ii] == p[i][i]
                for jj, mm in enumerate(types):
                    if ii == jj or not pcm[ii][jj]:
                        continue
                    old = next(k for k in range(q) if nn[k] == mm[k]+1)
                    new = next(k for k in range(q) if nn[k]+1 == mm[k])
                    mult = nn[old]*(nn[new]+1)
                    correct &= pcm[ii][jj]*pcm[jj][ii]/mult == self_p[i]*self_p[groups[jj][0]]
                    q_squared = mult*pm[jj]/pm[ii]
                    correct &= pcm[ii][jj]**2 >= q_squared/(bound*bound)
            ck(name+"_count_transport_multiplicities_not_dropped", correct and rowvec(pm, pcm) == pm)
            if name == "q3_Gibbs_N2":
                ck("Gibbs_fixture_actual_potential_ratios", set(w) == {15, 60} and c == 4 and m == 6 and mbar == 10 and bound == 41)
                mu = [F(1, 9)]*9
                gibbs_w = [F(4, 135)*t for t in w]
                law = mu[:]
                dom = all(mu[i] <= gibbs_w[i] for i in range(9)) and rowvec(gibbs_w, p) == gibbs_w
                for _ in range(4):
                    law = rowvec(law, p)
                    dom &= all(law[i] <= gibbs_w[i] for i in range(9))
                ck("actual_unaugmented_Gibbs_iid_weight_domination", dom)
        # A general scalar-step identity; h is a rational algebraic fixture.
        h = F(1, 29)
        centered, strict = True, False
        correction, expected_abs = F(0), F(0)
        maxjump = max(abs(f[j]-f[i]) for i, j in edges)
        for i in range(len(z)):
            mean = sum(p[i][j]*(f[j]-f[i]) for j in adj[i])
            second = sum(p[i][j]*(f[j]-f[i])**2 for j in adj[i])
            variance = sum(p[i][j]*(f[j]-f[i]-mean)**2 for j in adj[i])
            beta, raw = mean/h, second/h
            centered &= raw-h*beta*beta == variance/h >= 0
            strict |= mean != 0
            correction += pi[i]*h*beta*beta
            expected_abs += pi[i]*abs(beta)
        ck(name+"_actual_discrete_centered_variance_identity", centered and strict)
        ck(name+"_centering_paid_by_L1_drift_algebra", correction <= maxjump*expected_abs)
        summaries.append({"case": name, "q": q, "N": n, "labelled_states": len(z), "C_actual": str(c), "max_z_neighbour_ratio": str(m), "Mbar": str(mbar), "B": str(bound), "scope": "small exact source fixture, not actual A7/MT"})
    # Exact all-N colour formula at a few declared arguments, no large carrier.
    colour_data = []
    for n in (12, 24, 120):
        q = 12
        high, low, d = q*n*n, n*n+(q-2)*n+1, 1+n*(q-1)
        w = [F(high)] + [F(low)]*(q-1)
        pi = [x/sum(w) for x in w]
        p = [[F(1, q) for _ in range(q)]]
        for i in range(1, q):
            p.append([F(n*n, low) if j == 0 else F(1, low) if j == i else F(n, low) for j in range(q)])
        ck(f"colour_family_rows_balance_N{n}", all(sum(row) == 1 for row in p) and rowvec(pi, p) == pi and all(pi[i]*p[i][j] == pi[j]*p[j][i] for i, j in product(range(q), repeat=2)))
        ck(f"colour_family_fixed_ratio_bound_N{n}", 1 < F(high, low) < q)
        ck(f"colour_family_exact_small_activity_N{n}", F(d, low) < F(q-1, n) and F(d*n, low) >= F(q-1, 3))
        # The low-colour self loop is exactly 1/low; off-labelled low-low
        # activity is d/low, while the colour quotient has N parallel moves.
        ck(f"colour_family_distinguishes_quotient_multiplicity_N{n}", p[1][2] == F(n, low) and d*p[1][2]/n == F(d, low))
        mb = q+(q-1)*(n-1)*(q-1)
        bb = 1+n*(q-1)*mb
        ck(f"colour_family_obeys_polynomial_theorem_N{n}", n <= mb and low <= bb and F(d, low) >= F(d, bb))
        colour_data.append({"N": n, "high_weight": high, "low_weight": low, "normalized_low_low_activity": str(F(d, low)), "scope": "exact quotient of a distinct diagnostic stationary input, not A7 Gibbs"})
    n, q = 24, 12
    balanced = [n//q]*q
    x, y = balanced[:], balanced[:]
    x[0] -= 1
    x[1] += 1
    y[0] -= 1
    y[2] += 1
    colour = lambda nn: sum(a*t for a, t in enumerate(nn)) % q
    ck("colour_control_exact_near_uniform_count_geometry", colour(balanced) == 0 and colour(x) == 1 and colour(y) == 2 and sum(abs(a-b) for a, b in zip(x, y)) == 2 and min(x+y) >= 0)
    ck("colour_control_D12_and_slot_scope", all(sum(((a+s) % q)*t for a, t in enumerate(x)) % q == colour(x) for s in range(q)) and sum(((-a) % q)*t for a, t in enumerate(x)) % q == (-colour(x)) % q)
    ck("colour_input_order_opposite_to_A7_parity_lower_bound", F(n*n+10*n+1, 12*n*n) < 1 and F(12, 1)*F(4, 12)/(2*n) == F(2, n) > 0)
    # Verify actual colour quotient counting once on a tiny complete carrier.
    states = list(product(range(3), repeat=3))
    ck("colour_quotient_N_moves_per_other_colour", all([sum(sum(y) % 3 == c for y in neighbors(x, 3)) for c in range(3)] == [1 if c == sum(x) % 3 else 3 for c in range(3)] for x in states))
    for r in (12, 24):
        n, q = r**4, 12
        d = 1+n*(q-1)
        h, jump = F(1, d*r*r), F(2, r**3)
        beta_max = F(2*d, r)
        ck(f"actual_critical_grid_centering_scale_root{r}", h*beta_max == jump and jump*jump/h == F(4*d, n) <= 4*q)
    geom = json.loads((ROOT/"fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json").read_text())
    ck("actual_spectral_boundary_inputs_retained", sha(ROOT/geom["spectral_path"]) == geom["spectral_sha256"] and F(geom["b"]) == F(1, 2*10**10))
    return names, summaries, colour_data


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
    names = {p.relative_to(ROOT).as_posix() for d in PROTECTED for p in (ROOT/d).glob("*.md")}
    prior = json.loads((ROOT/"fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json").read_text())
    for name in prior["searched_files"]:
        p = Path(name)
        if len(p.parts) == 1 and re.search(r"(Entropy|Source|Hamming|Configuration|Generator|Sinkhorn|Scaling|Matrix|Reference)", p.name, re.I):
            names.add(name)
    names = sorted(names)
    pat = r"Hamming|neighbou?r.{0,30}ratio|one.coordinate.{0,25}clique|global.{0,35}scal|polynomial.{0,35}scal|Sinkhorn|unaugmented"
    matches = []
    for i in range(0, len(names), 90):
        r = subprocess.run(["rg", "-n", "-i", "--max-columns", "350", "--", pat, *names[i:i+90]], cwd=ROOT, capture_output=True, text=True, timeout=60)
        assert r.returncode in (0, 1), r.stderr
        matches.extend(r.stdout.splitlines())
    return {"selected_files": names, "files": len(names), "pattern": pat, "matching_lines": sorted(matches), "scope": "bounded selected local search, no exhaustive absence or global literature-priority claim", "reused": "entropy/matrix-scaling source, MT paths, SR rate/density, SQ interface", "new_dependency": "global polynomial scaling of the actual unaugmented source; stronger native constant remains open"}


def main():
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    inputs = {n: sha(ROOT/n) for n in INPUTS}
    before = {d: manifest(ROOT/d) for d in PROTECTED}
    checked, fixtures, colour = checks()
    duplicate = duplicate_search()
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    assert inputs == {n: sha(ROOT/n) for n in INPUTS}
    assert before == {d: manifest(ROOT/d) for d in PROTECTED}
    result = {"date": "2026-10-20", "status": "PASS_SCOPED_IS_CONTROLS", "bounded_exact_controls": len(checked), "checks": checked, "finite_source_fixtures": fixtures, "colour_countercontrols": colour, "arithmetic": "Fraction only; no binary64", "source_inputs": inputs, "protected_packages": before, "protected_unchanged": True, "AGENTS_full_sha256": inputs["AGENTS.md"], "original_A7_constant_target": "OPEN_NOT_REFUTED_BY_DIFFERENT_STATIONARY_INPUT", "actual_A7_critical_defects": "OPEN", "scope": "finite exact controls and integrity, not mechanical all-N proof or independent review", "old_scientific_writers_run": [], "internet_acquisition": False, "strict_F01": "OPEN_SOURCE", "broad_goal": "ACTIVE_NOT_COMPLETE"}
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
    print(json.dumps({"status": result["status"], "bounded_exact_controls": len(checked), "inputs": len(inputs), "protected_packages": len(before), "protected_entries": sum(v["entries"] for v in before.values()), "duplicate_search_files": duplicate["files"], "AGENTS_or_source_changed": False}, indent=2))


if __name__ == "__main__":
    main()
