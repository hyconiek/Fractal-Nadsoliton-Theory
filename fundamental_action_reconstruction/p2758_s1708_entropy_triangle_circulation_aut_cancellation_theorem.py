#!/usr/bin/env python3
"""P2758/S1708: entropy triangle-circulation Aut-cancellation theorem.

P2757 closed intrinsic sign-blind selectors inside local entropy pair-current
observables.  This pivot tests a genuinely different typed object: a
translation-summed oriented three-point/triangle entropy circulation on Z12.
For every alternating local law B(a,b,c) on the five-level four-quanta entropy
alphabet, the oriented triangle feature vector for step -u is the exact negative
of the feature vector for step u.  Therefore selector-free Aut(Z12) averaging
cancels the entire class unless an independent strict orientation/polarity law
and an explicit P2721 coupling theorem are exported.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from itertools import combinations
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

N = 12
QUANTA = 4
UNITS = (1, 5, 7, 11)
OPPOSITE = {1: 11, 11: 1, 5: 7, 7: 5}
GEN = ROOT / "generated"
P2757 = GEN / "p2757_s1707_entropy_pair_current_polarity_selector_invariant_no_go.json"
OUT = GEN / "p2758_s1708_entropy_triangle_circulation_aut_cancellation_theorem.json"
MD = GEN / "p2758_s1708_entropy_triangle_circulation_aut_cancellation_theorem.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "post_p2757_pivot_obligation": r"P2757|outside local entropy pair-current|strict step/polarity source|P2697-P2757",
    "selector_boundary": r"P2721|lambda/P2721|QW-2191|orientation torsor|selector closure",
    "aut_boundary": r"Aut\(Z12\)|opposite step|1/11|5/7|polarity",
    "closure_forbidden": r"role transfer|L_total|ToE closure|bridge closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "entropy_triangle_circulation_selector_exported",
    "p2721_triangle_circulation_coupling_theorem_exported",
    "strict_orientation_or_polarity_law_exported",
    "lambda_p2721_fixed",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def run_rg(pattern: str) -> list[str]:
    cmd = ["rg", "-n", "--glob", "!generated/*.json", pattern, "AGENTS.md", "fundamental_action_reconstruction"]
    proc = subprocess.run(cmd, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed for {pattern!r}: {proc.stderr}")
    return [line for line in proc.stdout.splitlines() if line.strip()]


def evidence_scan() -> dict[str, Any]:
    rows = []
    for name, pattern in CONTENT_PATTERNS.items():
        hits = run_rg(pattern)
        rows.append({"content_lane": name, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return {"content_pattern_count": len(rows), "rows": rows, "hit_counts": {r["content_lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def triangle_basis() -> list[tuple[int, int, int]]:
    return list(combinations(range(QUANTA + 1), 3))


def permutation_sign_to_sorted(triple: tuple[int, int, int]) -> tuple[tuple[int, int, int] | None, int]:
    if len(set(triple)) < 3:
        return None, 0
    sorted_triple = tuple(sorted(triple))
    inversions = 0
    vals = list(triple)
    for i in range(3):
        for j in range(i + 1, 3):
            if vals[i] > vals[j]:
                inversions += 1
    return sorted_triple, -1 if inversions % 2 else 1


def feature_vector(levels: tuple[int, ...], step: int, basis: list[tuple[int, int, int]]) -> tuple[int, ...]:
    index = {triple: k for k, triple in enumerate(basis)}
    vec = [0] * len(basis)
    for i, a in enumerate(levels):
        triple = (a, levels[(i + step) % N], levels[(i + 2 * step) % N])
        key, sign = permutation_sign_to_sorted(triple)
        if key is not None:
            vec[index[key]] += sign
    return tuple(vec)


def rank_int(rows: list[tuple[int, ...]]) -> int:
    from fractions import Fraction
    mat = [[Fraction(x) for x in row] for row in rows if any(row)]
    if not mat:
        return 0
    m, n = len(mat), len(mat[0])
    rank = col = 0
    while rank < m and col < n:
        pivot = next((r for r in range(rank, m) if mat[r][col] != 0), None)
        if pivot is None:
            col += 1
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        pv = mat[rank][col]
        mat[rank] = [x / pv for x in mat[rank]]
        for r in range(m):
            if r != rank and mat[r][col] != 0:
                factor = mat[r][col]
                mat[r] = [x - factor * y for x, y in zip(mat[r], mat[rank])]
        rank += 1
        col += 1
    return rank


def triangle_circulation_audit() -> dict[str, Any]:
    basis = triangle_basis()
    all_directed_rows: list[tuple[int, ...]] = []
    all_aut_rows: list[tuple[int, ...]] = []
    nonzero_by_step = {u: 0 for u in UNITS}
    opposite_failures = []
    aut_failures = []
    sample_nonzero = []
    for counts in compositions(QUANTA, N):
        levels = tuple(counts)
        vectors = {u: feature_vector(levels, u, basis) for u in UNITS}
        for u, vec in vectors.items():
            all_directed_rows.append(vec)
            if any(vec):
                nonzero_by_step[u] += 1
                if len(sample_nonzero) < 8:
                    sample_nonzero.append({"counts": list(counts), "step": u, "feature_vector": list(vec), "opposite_step": OPPOSITE[u], "opposite_feature_vector": list(vectors[OPPOSITE[u]])})
        for u in (1, 5):
            paired_sum = tuple(a + b for a, b in zip(vectors[u], vectors[OPPOSITE[u]]))
            if any(paired_sum) and len(opposite_failures) < 8:
                opposite_failures.append({"counts": list(counts), "step": u, "paired_sum": list(paired_sum)})
        aut_sum = tuple(sum(vectors[u][j] for u in UNITS) for j in range(len(basis)))
        all_aut_rows.append(aut_sum)
        if any(aut_sum) and len(aut_failures) < 8:
            aut_failures.append({"counts": list(counts), "aut_sum": list(aut_sum), "vectors": {str(k): list(v) for k, v in vectors.items()}})
    return {
        "modulus": N,
        "quanta": QUANTA,
        "composition_count": math.comb(QUANTA + N - 1, N - 1),
        "basis_dimension": len(basis),
        "basis": [list(t) for t in basis],
        "nonzero_by_step": nonzero_by_step,
        "directed_feature_rank": rank_int(all_directed_rows),
        "aut_summed_rank": rank_int(all_aut_rows),
        "opposite_failure_count": len(opposite_failures),
        "aut_sum_failure_count": len(aut_failures),
        "sample_nonzero_rows": sample_nonzero,
        "sample_opposite_failures": opposite_failures,
        "sample_aut_failures": aut_failures,
        "theorem_statement": "For every translation-summed oriented three-point entropy circulation with alternating local law B(a,b,c), replacing u by -u reverses the cyclic order of each triangle after reindexing, hence the full feature vector satisfies C(-u)=-C(u).  Aut(Z12) pairs 1/11 and 5/7, so selector-free Aut averaging cancels the entire alternating triangle-circulation class.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_post_p2757_pivot_obligation": scan["all_patterns_have_hits"],
        "triangle_circulation_basis_is_nontrivial": audit["basis_dimension"] == 10 and audit["directed_feature_rank"] > 0,
        "opposite_steps_cancel_exactly": audit["opposite_failure_count"] == 0,
        "selector_free_aut_sum_rank_is_zero": audit["aut_summed_rank"] == 0 and audit["aut_sum_failure_count"] == 0,
        "strict_orientation_or_polarity_law_exported": False,
        "p2721_triangle_circulation_coupling_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_strict_selector_source": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The oriented triangle entropy-circulation class is nontrivial when a direction is supplied, but Aut(Z12) pairs opposite directions and cancels the whole alternating three-point basis; no strict orientation/polarity law or P2721 coupling theorem is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["entropy_triangle_circulation_audit"]
    lines = [
        "# P2758/S1708 entropy triangle-circulation Aut-cancellation theorem",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite triangle-circulation audit",
        f"- modulus={a['modulus']}",
        f"- quanta={a['quanta']}",
        f"- composition_count={a['composition_count']}",
        f"- basis_dimension={a['basis_dimension']}",
        f"- directed_feature_rank={a['directed_feature_rank']}",
        f"- aut_summed_rank={a['aut_summed_rank']}",
        f"- opposite_failure_count={a['opposite_failure_count']}",
        f"- aut_sum_failure_count={a['aut_sum_failure_count']}",
        f"- nonzero_by_step={a['nonzero_by_step']}",
        "",
        "## Theorem",
        a["theorem_statement"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2757 = read_json(P2757)
    scan = evidence_scan()
    audit = triangle_circulation_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2758_ENTROPY_TRIANGLE_CIRCULATION_AUT_CANCELLATION_THEOREM",
        "input_hashes": {"P2757_ENTROPY_PAIR_CURRENT_POLARITY_SELECTOR_INVARIANT_NO_GO": sha(P2757)},
        "input_statuses": {"P2757_ENTROPY_PAIR_CURRENT_POLARITY_SELECTOR_INVARIANT_NO_GO": p2757.get("status")},
        "audited_candidate_class": "translation-summed oriented three-point alternating Shannon entropy triangle circulations on Z12",
        "content_evidence_scan": scan,
        "entropy_triangle_circulation_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not continue the entropy-selector lane by escalating from pair-currents to local oriented triangle/three-point entropy circulations: P2758 proves the full alternating triangle basis is nontrivial only after a direction is supplied, and selector-free Aut(Z12) averaging cancels it exactly.  The next proof-grade move must either export an independent strict orientation/polarity law with an explicit P2721 coupling theorem, or pivot to a genuinely new typed object outside local entropy observables; otherwise preserve the P2697-P2758 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2758/S1708 entropy triangle-circulation Aut-cancellation theorem", "## P2758/S1708 entropy triangle-circulation Aut-cancellation theorem\n\n`P2758/S1708` pivots outside local pair-current observables and audits translation-summed oriented three-point Shannon entropy triangle circulations on the four-quanta `Z12` entropy alphabet.  The full alternating local triangle basis has dimension `10` and directed rank is nonzero, so the typed object is not empty.  However, replacing a step `u` by `-u` reverses the triangle orientation after reindexing, giving exact vector negation for the `1/11` and `5/7` pairs; the selector-free `Aut(Z12)`-summed triangle-circulation rank is `0` with zero opposite-pair and Aut-sum failures.  Thus no strict orientation/polarity law, `P2721` coupling theorem, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, bridge closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2758/S1708 entropy triangle-circulation Ltotal guard", "## P2758/S1708 entropy triangle-circulation Ltotal guard\n\n`P2758/S1708` adds no variational source term.  It shows that the entire local alternating entropy triangle-circulation basis is cancelled by selector-free `Aut(Z12)` handling, despite being nontrivial once a directed step is supplied.  Therefore it cannot promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure without an independent strict orientation/polarity law and explicit `P2721` coupling theorem.\n")
    append_once(AGENTS, "Current entropy triangle-circulation Aut-cancellation guardrail (P2758/S1708, 2026-06-15)", "## Current entropy triangle-circulation Aut-cancellation guardrail (P2758/S1708, 2026-06-15)\n\n- P2758 pivots beyond local entropy pair-current observables by auditing translation-summed oriented three-point Shannon entropy triangle circulations on the four-quanta `Z12` entropy alphabet.\n- The alternating triangle basis is nontrivial when a direction is supplied, but opposite directed steps are exact negatives and the selector-free `Aut(Z12)`-summed triangle-circulation basis has rank `0` with zero opposite-pair and Aut-sum failures.\n- Do not continue the entropy-selector lane by escalating to local oriented triangle/three-point entropy circulations as `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure.  A next admissible move must export an independent strict orientation/polarity law with explicit `P2721` coupling, pivot to a genuinely new typed object outside local entropy observables, or preserve the P2697-P2758 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
