#!/usr/bin/env python3
"""Tiny exact controls; no claim to enumerate native FIN or prove limits."""
import argparse
import ast
from collections import Counter
from functools import lru_cache
from fractions import Fraction as F
import hashlib
from itertools import permutations, product
import json
from math import factorial, prod
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "FIN_RESEARCH_EXECUTION_PLAN.md",
    "fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md",
    "fin_hyp_history_entropy_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_hyp_history_entropy_review_20261020/PROOF_AUDIT_AND_DUAL.md",
    "fin_hyp_history_entropy_review_20261020/ENTROPY_DEFICIT_PATH_THEOREM.md",
    "fin_reference_transport_research_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_intrinsic_scaling_research_20261020/GLOBAL_SCALING_THEOREM.md",
    "fin_global_activity_defect_research_20261020/ORIGINAL_SOURCE_PATH_THEOREM.md",
    "fin_global_activity_defect_research_20261020/ACTUAL_DRIFT_NOISE_PAYMENTS.md",
    "fin_global_activity_defect_self_stress_20261020/USER_DIRECTION.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "FIN_RESEARCH_AFTER_256_294_FULL_HANDOFF_295_299_20260927/01_REPORTS/MICROSCOPIC_PERIODIC_BRIDGE_296.md",
]
PROTECTED = [
    "fin_hyp_history_entropy_20261020",
    "fin_hyp_history_entropy_review_20261020",
    "fin_reference_transport_research_20261020",
    "fin_intrinsic_scaling_research_20261020",
    "fin_global_activity_defect_research_20261020",
    "fin_global_activity_defect_self_stress_20261020",
    "fin_critical_metastability_review_20261015",
    "fin_uniform_critical_review_20261014",
    "fin_intrinsic_observation_review_20261015",
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder):
    seen = set()
    for line in (folder / "MANIFEST.sha256").read_text().splitlines():
        h, name = line.split("  ", 1)
        assert not Path(name).is_absolute() and ".." not in Path(name).parts
        assert name not in seen and sha(folder / name) == h, (folder, name)
        seen.add(name)
    actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*")
              if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
    actual |= {n for n in seen if "__pycache__" in Path(n).parts or n.endswith(".pyc")}
    assert actual == seen | {"MANIFEST.sha256"}, folder
    return {"entries": len(seen), "sha256": sha(folder / "MANIFEST.sha256")}


def validate(K, n):
    m = len(n)
    if not m or any(type(x) is not int or x < 0 for x in n) or not sum(n):
        raise ValueError("nonempty nonnegative integer occupation with L>0 required")
    if len(K) != m or any(len(row) != m for row in K):
        raise ValueError("square carrier reference required")
    if any(type(K[i][j]) is not int or K[i][j] not in (0, 1)
           or K[i][j] != K[j][i] for i in range(m) for j in range(m)):
        raise ValueError("symmetric unweighted reference required")
    if any(K[i][i] != 1 for i in range(m)):
        raise ValueError("one self continuation required")


def words(K, n):
    validate(K, n)
    m, L = len(n), sum(n)
    assert m**L <= 20000, "bounded diagnostic only"
    return [w for w in product(range(m), repeat=L)
            if tuple(w.count(i) for i in range(m)) == tuple(n)
            and all(K[w[t]][w[(t+1) % L]] for t in range(L))]


def completions(K, root):
    m = len(K)
    if type(root) is not int or not 0 <= root < m:
        raise ValueError("invalid root")

    @lru_cache(None)
    def C(x, r):
        if sum(r) == 0:
            return int(x == root)
        if not r[x]:
            return 0
        rem = tuple(r[i] - (i == x) for i in range(m))
        return sum(C(y, rem) for y in range(m) if K[x][y])

    def next_law(x, r):
        if not sum(r) or not C(x, tuple(r)):
            return None
        rem = tuple(r[i] - (i == x) for i in range(m))
        den = C(x, tuple(r))
        return tuple(F(K[x][y] * C(y, rem), den) for y in range(m))
    return C, next_law


def source(K, a):
    m = len(a)
    W = [a[i] * sum(K[i][j]*a[j] for j in range(m)) for i in range(m)]
    total = sum(W)
    pi = tuple(F(x, total) for x in W)
    c = tuple(tuple(F(K[i][j]*a[i]*a[j], total) for j in range(m)) for i in range(m))
    P = tuple(tuple(c[i][j]/pi[i] for j in range(m)) for i in range(m))
    return pi, c, P


def matrix_inverse(A):
    n = len(A)
    a = [[F(x) for x in A[i]] + [F(i == j) for j in range(n)] for i in range(n)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        d = a[j][j]
        a[j] = [x/d for x in a[j]]
        for i in range(n):
            if i != j:
                d = a[i][j]
                a[i] = [x-d*y for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def empirical_table(w, m):
    out = [[0]*m for _ in range(m)]
    for t, x in enumerate(w):
        out[x][w[(t+1) % len(w)]] += 1
    return tuple(tuple(row) for row in out)


def last_exit_words(M, root, parents):
    m, L = len(M), sum(map(sum, M))
    lists = []
    for x in range(m):
        destinations = tuple(y for y in range(m) for _ in range(M[x][y]))
        if x == root:
            lists.append(sorted(set(permutations(destinations))))
        else:
            prefix = list(destinations)
            prefix.remove(parents[x])
            lists.append([p+(parents[x],) for p in sorted(set(permutations(prefix)))])
    result = []
    for rows in product(*lists):
        used, w, x = [0]*m, [], root
        for _ in range(L):
            assert used[x] < len(rows[x]), "stalled row stack"
            w.append(x)
            nxt = rows[x][used[x]]
            used[x] += 1
            x = nxt
        assert x == root and used == list(map(sum, M))
        assert empirical_table(w, m) == M
        result.append(tuple(w))
    return result


def _unit_log(x, count=60):
    assert 1 <= x <= 2
    w = (x-1)/(x+1)
    value = 2*sum((w**(2*j+1)/F(2*j+1) for j in range(count)), F(0))
    tail = 2*w**(2*count+1)/(F(2*count+1)*(1-w*w))
    return value, value+tail


@lru_cache(None)
def log_box(x):
    """Exact atanh series and tail, public outward 80-bit dyadic endpoints."""
    x = F(x)
    assert x > 0
    k = 0
    while x < 1:
        x *= 2
        k -= 1
    while x >= 2:
        x /= 2
        k += 1
    lo, hi = _unit_log(x)
    l2, u2 = _unit_log(F(2))
    lo, hi = (lo+k*l2, hi+k*u2) if k >= 0 else (lo+k*u2, hi+k*l2)
    scale = 1 << 80
    return F((lo*scale).__floor__(), scale), F((hi*scale).__ceil__(), scale)


def kl_box(p, q):
    ratios = Counter()
    for key, mass in p.items():
        if mass:
            assert q.get(key, 0) > 0
            ratios[mass/q[key]] += mass
    lo, hi = F(0), F(0)
    for ratio, mass in ratios.items():
        l, u = log_box(ratio)
        lo += mass*l
        hi += mass*u
    return lo, hi


def box_json(box):
    scale = 1 << 80
    return {"bits": 80, "lower": int((box[0]*scale).__floor__()),
            "upper": int((box[1]*scale).__ceil__())}


def block(Uwords, k):
    L = len(Uwords[0])
    count = Counter(tuple(w[t % L] for t in range(k+1)) for w in Uwords)
    return {w: F(n, len(Uwords)) for w, n in count.items()}


def markov_block(pi, P, k):
    m = len(pi)
    return {w: pi[w[0]] * prod(P[w[t]][w[t+1]] for t in range(k))
            for w in product(range(m), repeat=k+1)}


def check_links():
    for p in HERE.glob("*.md"):
        prose, fence = [], False
        for line in p.read_text().splitlines():
            assert line == line.rstrip(), p.name
            if line.lstrip().startswith("```"):
                fence = not fence
            elif not fence and not line.startswith(("    ", "\t")):
                prose.append(line)
        txt = re.sub(r"`[^`\n]*`", "", "\n".join(prose))
        for name in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", txt):
            name = name.strip("<>").split("#")[0]
            if name and not name.startswith(("https:", "http:")):
                assert (p.parent/name).exists(), (p.name, name)


def run():
    checks, details = [], {}
    def ck(name, condition):
        assert name not in checks and condition, name
        checks.append(name)
    K = ((1, 1, 0), (1, 1, 1), (0, 1, 1))
    n, L = (2, 3, 2), 7
    pi, c, P = source(K, (1, 1, 1))
    W = words(K, n)
    ck("noncomplete_reference_exact_nonuniform_source", pi == (F(2,7), F(3,7), F(2,7)) and
       all(sum(row) == 1 for row in P) and c[0][2] == 0)
    ck("frozen_occupation_class_nonempty", bool(W))
    trace_weight = prod(P[i][i]**n[i] for i in range(3))
    ck("all_transition_products_equal", all(prod(P[w[t]][w[(t+1)%L]] for t in range(L)) == trace_weight for w in W))
    periodic = words(((1,1),(1,1)), (2,2))
    orbit_sizes = {len({w[t:]+w[:t] for t in range(4)}) for w in periodic}
    ck("rotation_including_periodic_stabilizers", {w[1:]+w[:1] for w in W} == set(W)
       and len(periodic) == 6 and orbit_sizes == {2,4})
    ck("reversal_preserves_same_class", {tuple(reversed(w)) for w in W} == set(W))
    counts = Counter(w[0] for w in W)
    ck("root_marginal_exact_occupation", all(F(counts[x], len(W)) == pi[x] for x in range(3)))
    for a in range(3):
        C, law = completions(K, a)
        ck(f"completion_recursion_root_{a}", C(a, n) == counts[a])
        probabilities = []
        for w in W:
            if w[0] != a:
                continue
            rem, probability = n, F(1)
            for t in range(L):
                row = law(w[t], rem)
                assert row is not None and sum(row) == 1
                probability *= row[w[(t+1)%L]]
                rem = tuple(rem[i]-(i == w[t]) for i in range(3))
            probabilities.append(probability)
        ck(f"causal_probabilities_telescope_root_{a}", all(p == F(1, counts[a]) for p in probabilities))
    C, law = completions(K, 0)
    ck("impossible_support_returns_none", C(0, (1,0,1)) == 0 and law(0, (1,0,1)) is None)
    ck("terminal_zero_has_no_next_transition", C(0, (0,0,0)) == 1 and law(0, (0,0,0)) is None)
    for label, bad in [("noninteger", (1,F(3,2),1)), ("zero_length", (0,0,0)), ("negative", (1,-1,2)), ("boolean", (True,1,1))]:
        rejected = False
        try:
            validate(K, bad)
        except ValueError:
            rejected = True
        ck("reject_"+label, rejected)
    C2, law2 = completions(((1,1),(1,1)), 0)
    binary_memory = words(((1,1),(1,1)), (3,2))
    prefix_a = [w for w in binary_memory if w[:4] == (0,0,1,0)]
    prefix_b = [w for w in binary_memory if w[:4] == (0,1,1,0)]
    ck("finite_record_not_state_only_markov", law2(0, (1,1)) == (F(0),F(1))
       and law2(0, (2,0)) == (F(1),F(0)) and prefix_a and prefix_b
       and {w[4] for w in prefix_a} == {1} and {w[4] for w in prefix_b} == {0})
    pi2, _, P2 = source(K, (1,2,1))
    path = (0,1,1)
    product_open = P2[0][1]*P2[1][1]
    diag_product = P2[0][0]*P2[1][1]
    ck("open_path_endpoint_ratio_required", product_open == 2*diag_product and product_open != diag_product)
    M = tuple(tuple(K[i][j] for j in range(3)) for i in range(3))
    stacks = last_exit_words(M, 0, {1:0, 2:1})
    ck("last_exit_tree_consumes_and_is_injective", len(stacks) == len(set(stacks)) == 4)
    table_words = [w for w in W if w[0] == 0 and empirical_table(w, 3) == M]
    ck("last_exit_lists_match_all_path_tree_words", set(stacks) == set(table_words))
    lists_formula = prod(F(factorial(sum(row)), prod(factorial(x) for x in row)) for row in M)
    lists_formula *= prod(F(M[x][{1:0,2:1}[x]], sum(M[x])) for x in (1,2))
    ck("row_multinomial_tree_factor_exact", lists_formula == len(stacks))
    badM = ((2,0),(0,2))
    binary_words = words(((1,1),(1,1)), (2,2))
    ck("balanced_table_without_connectivity_is_not_word", not any(empirical_table(w, 2) == badM for w in binary_words))
    ck("rounding_table_exact_native_margins_fixture", M == tuple(tuple(int(L*x) for x in row) for row in c) and tuple(map(sum,M)) == n)
    # A second table, with genuine floor errors, verifies the universal squared bound.
    LL = 45
    cc = ((F(1,5),F(1,10),F(1,30)), (F(1,10),F(1,5),F(1,30)), (F(1,30),F(1,30),F(4,15)))
    assert sum(map(sum,cc)) == 1
    nn = tuple(LL*sum(row) for row in cc)
    # Its exact row margins are thirds; no approximate drift/row repair.
    assert all(x.denominator == 1 for x in nn)
    MM = [[0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(i+1,3):
            MM[i][j] = MM[j][i] = int(LL*cc[i][j])
    for i in range(3):
        MM[i][i] = int(nn[i])-sum(MM[i])
    err2 = sum((F(MM[i][j],LL)-cc[i][j])**2 for i in range(3) for j in range(3))
    ck("symmetric_floor_and_squared_error_bound", tuple(map(sum,MM)) == tuple(nn) and err2 <= F(18,LL**2))
    # Log boxes are only used to enclose finite diagnostic divergences.
    ck("log_zero_and_outward_reciprocal", log_box(F(1)) == (0,0) and
       log_box(F(3,2))[0]+log_box(F(2,3))[0] <= 0 <= log_box(F(3,2))[1]+log_box(F(2,3))[1])
    Ufull = {w+(w[0],): F(1,len(W)) for w in W}
    Qfull = {w: pi[w[0]]*trace_weight for w in Ufull}
    direct = kl_box(Ufull, Qfull)
    hbox = tuple(-sum(pi[i]*log_box(P[i][i])[1-j] for i in range(3)) for j in range(2))
    Hbox = tuple(-sum(pi[i]*log_box(pi[i])[1-j] for i in range(3)) for j in range(2))
    zbox = log_box(F(len(W)))
    identity = (L*hbox[0]-zbox[1]+Hbox[0], L*hbox[1]-zbox[0]+Hbox[1])
    ck("full_KL_identity_with_root_entropy", max(direct[0],identity[0]) <= min(direct[1],identity[1]))
    Aupper = 11*log_box(F(L+1))[1]+F(18,L*F(1,7))+Hbox[1]
    ck("finite_counting_upper_KL_payment", direct[1] <= Aupper)
    blocks = []
    for k in (1,2,3):
        u, q = block(W,k), markov_block(pi,P,k)
        db = kl_box(u,q)
        ck(f"whole_block_chain_rule_k_{k}", (L//k)*db[1] <= direct[0])
        pD = pi[0]
        uc = {w:p/pD for w,p in u.items() if w[0] == 0}
        qc = {w:p/pD for w,p in q.items() if w[0] == 0}
        cond = kl_box(uc,qc)
        ck(f"rare_initial_partition_KL_k_{k}", pD*cond[1] <= db[0])
        blocks.append({"transitions": k, "KL_dyadic_box":box_json(db),
                       "conditional_root0_KL_dyadic_box":box_json(cond)})
    complete = ((1,1,1),)*3
    nn = (3,2,1)
    WW = words(complete, nn)
    pp = tuple(F(x,6) for x in nn)
    bias = [sum(pp[w[0]] for w in WW if w[0] == x) for x in range(3)]
    total = sum(bias)
    bias = tuple(x/total for x in bias)
    ck("stationary_closing_root_bias_counterexample", bias == (F(9,14),F(2,7),F(1,14)) and bias != pp)
    # Exact inverse validates the continuity operator algebra, not an all-N theorem.
    inv = matrix_inverse([[F(i == j)+P[i][j] for j in range(3)] for i in range(3)])
    ck("I_plus_P_inverse_residual", all(sum((F(i==t)+P[i][t])*inv[t][j] for t in range(3)) == F(i==j) for i in range(3) for j in range(3)))
    ck("self_loop_inverse_infinity_bound", max(sum(abs(x) for x in row) for row in inv) <= F(3,2))
    # N=1 Hamming control; synthetic rational pi, not an A7 critical fixture.
    pih, ch, Ph = source(complete,(1,2,3))
    Cedge, B = F(3), 7
    ck("hamming_IS_lower_self_fixture", min(Ph[i][i] for i in range(3)) >= F(1,B))
    ck("hamming_IS_lower_flow_weaker_rational_fixture", min(map(min,ch)) >= min(pih)/(B*Cedge))
    # Rational largest remainder: exact ties broken by fixed state index.
    occupation_target = (F(1,7),F(2,7),F(4,7))
    length = 20
    floors = [int(length*x) for x in occupation_target]
    order = sorted(range(3), key=lambda i:(-(length*occupation_target[i]-floors[i]),i))
    for i in order[:length-sum(floors)]:
        floors[i] += 1
    ck("largest_remainder_exact_total_and_coordinate_error", sum(floors) == length and all(abs(F(floors[i],length)-occupation_target[i]) <= F(1,length) for i in range(3)))
    ck("native_diagonal_length_exact_integer", 12**(12**2) > (12**12)**2 and 12**(12**2) == pow(12,144))
    details.update({"carrier":"three-state path with one self at each state; not native A7",
                    "occupation":list(n),"cycle_length":L,"rooted_words":len(W),
                    "root_counts":{str(x):value for x,value in sorted(counts.items())},"last_exit_words":len(stacks),
                    "trace_weight":str(trace_weight),"full_KL_dyadic_box":box_json(direct),
                    "block_controls":blocks,"biased_root_control":[str(x) for x in bias],
                    "largest_remainder_control":floors,
                    "finite_native_graph_enumerated":False})
    return checks, details


def main():
    ap = argparse.ArgumentParser()
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--record",action="store_true")
    group.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path)
    args = ap.parse_args()
    if args.check:
        manifest(HERE)
    inputs = {name:sha(ROOT/name) for name in INPUTS}
    protected = {d:manifest(ROOT/d) for d in PROTECTED}
    checks, details = run()
    assert inputs == {name:sha(ROOT/name) for name in INPUTS}
    assert protected == {d:manifest(ROOT/d) for d in PROTECTED}
    result = {"date":"2026-10-20","status":"CONDITIONAL_FINITE_HISTORY_AND_NATIVE_DIAGONAL_THEOREMS_WITH_BOUNDED_CONTROLS",
              "bounded_checks":len(checks),"checks":checks,"finite_controls":details,
              "input_hashes":inputs,"protected_manifests":protected,"sources_unchanged":True,
              "same_author":True,"independent_peer_review":False,"subagents_used":False,
              "new_physical_postulate_adopted":False,"strict_F01":"OPEN_SOURCE",
              "AGENTS_modified":False,"broad_goal":"ACTIVE_NOT_COMPLETE"}
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    check_links()
    input_text = "".join(f"{h}  {n}\n" for n,h in inputs.items())
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result,indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(input_text)
        files = sorted(p for p in HERE.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        assert json.loads((HERE/"RESULTS.json").read_text()) == result
        assert (HERE/"INPUTS.sha256").read_text() == input_text
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"bounded_checks":len(checks),
                      "rooted_words":details["rooted_words"],"inputs":len(inputs),
                      "protected_manifests":len(protected),"protected_entries":sum(v["entries"] for v in protected.values()),
                      "sources_and_AGENTS_unchanged":True,"independent_audit":False},indent=2))


if __name__ == "__main__":
    main()
