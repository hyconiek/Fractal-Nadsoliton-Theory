#!/usr/bin/env python3
"""P2755/S1705: entropy-gradient current Aut-cancellation audit.

P2754 verified the scalar four-bit Shannon fact but showed scalar entropy has
zero equivariant maps to the orientation torsor.  This follow-up tests the next
more dynamical entropy idea: an inversion-odd entropy current/gradient/flux on
Z12.  The audit deliberately constructs a nonzero directed entropy current, then
checks whether the current becomes strict/selector-free after Aut(Z12) handling.
It does not: choosing a single directed edge step is already a selector premise,
and Aut-averaging over the exported units pairs each step with its opposite and
cancels the signed current.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

N = 12
QUANTA = 4
UNITS = (1, 5, 7, 11)
OPPOSITE = {1: 11, 11: 1, 5: 7, 7: 5}
TOL = 1e-12
GEN = ROOT / "generated"
P2754 = GEN / "p2754_s1704_shannon_entropy_four_bit_selector_audit.json"
OUT = GEN / "p2755_s1705_entropy_gradient_current_aut_cancellation_audit.json"
MD = GEN / "p2755_s1705_entropy_gradient_current_aut_cancellation_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "post_p2754_entropy_current_obligation": r"P2754|inversion-odd entropy current|entropy current|entropy.*flux|P2697-P2754",
    "selector_boundary": r"P2721|lambda/P2721|QW-2191|orientation torsor|selector closure",
    "aut_z12_boundary": r"Aut\(Z12\)|Aut\(Z_12\)|inversion|units `7` and `11`|orientation-reversing",
    "closure_forbidden": r"role transfer|L_total|ToE closure|bridge closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "entropy_current_promoted_to_selector_source",
    "directed_step_promoted_without_selector_premise",
    "aut_averaged_entropy_current_nonzero",
    "p2721_entropy_current_coupling_theorem_exported",
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


def entropy_density(counts: tuple[int, ...]) -> tuple[float, ...]:
    total = sum(counts)
    if total == 0:
        return tuple(0.0 for _ in counts)
    out = []
    for c in counts:
        if c == 0:
            out.append(0.0)
        else:
            p = c / total
            out.append(-p * math.log(p))
    return tuple(out)


def directed_entropy_current(h: tuple[float, ...], step: int) -> float:
    """A concrete antisymmetric local entropy flux along a directed Z12 step.

    J_step(h)=sum_i h_i*h_{i+step}*(h_i-h_{i+step}).
    It is intentionally odd under step -> -step, so it is the right kind of
    toy current to test.  If even this object cannot survive Aut handling
    without selecting a direction, scalar-entropy replay has not solved P2721.
    """
    return sum(h[i] * h[(i + step) % N] * (h[i] - h[(i + step) % N]) for i in range(N))


def sign(x: float) -> int:
    return 1 if x > TOL else -1 if x < -TOL else 0


def current_audit() -> dict[str, Any]:
    rows = []
    nonzero_by_step = {u: 0 for u in UNITS}
    positive_by_step = {u: 0 for u in UNITS}
    negative_by_step = {u: 0 for u in UNITS}
    opposite_failures = []
    aut_average_failures = []
    sample_nonzero = []
    for counts in compositions(QUANTA, N):
        h = entropy_density(counts)
        currents = {u: directed_entropy_current(h, u) for u in UNITS}
        for u, value in currents.items():
            s = sign(value)
            if s:
                nonzero_by_step[u] += 1
                if s > 0:
                    positive_by_step[u] += 1
                else:
                    negative_by_step[u] += 1
                if len(sample_nonzero) < 8:
                    sample_nonzero.append({"counts": list(counts), "step": u, "current": value, "sign": s, "opposite_step": OPPOSITE[u], "opposite_current": currents[OPPOSITE[u]]})
        for u in (1, 5):
            if abs(currents[u] + currents[OPPOSITE[u]]) > TOL and len(opposite_failures) < 8:
                opposite_failures.append({"counts": list(counts), "step": u, "current": currents[u], "opposite_step": OPPOSITE[u], "opposite_current": currents[OPPOSITE[u]]})
        aut_average = sum(currents.values()) / len(UNITS)
        if abs(aut_average) > TOL and len(aut_average_failures) < 8:
            aut_average_failures.append({"counts": list(counts), "currents": currents, "aut_average": aut_average})
    total = math.comb(QUANTA + N - 1, N - 1)
    for u in UNITS:
        rows.append({
            "step": u,
            "opposite_step": OPPOSITE[u],
            "nonzero_current_count": nonzero_by_step[u],
            "positive_current_count": positive_by_step[u],
            "negative_current_count": negative_by_step[u],
        })
    return {
        "modulus": N,
        "quanta": QUANTA,
        "composition_count": total,
        "current_formula": "J_u(h)=sum_i h_i*h_{i+u}*(h_i-h_{i+u}) for Shannon density h_i=-p_i log p_i",
        "directed_step_rows": rows,
        "has_nonzero_directed_entropy_currents": any(nonzero_by_step.values()),
        "opposite_pair_failure_count": len(opposite_failures),
        "sample_opposite_pair_failures": opposite_failures,
        "aut_average_failure_count": len(aut_average_failures),
        "sample_aut_average_failures": aut_average_failures,
        "sample_nonzero_directed_currents": sample_nonzero,
        "all_opposite_steps_cancel": len(opposite_failures) == 0,
        "aut_averaged_current_identically_zero": len(aut_average_failures) == 0,
        "proof_obstruction": "The concrete entropy current is nonzero only after choosing a directed unit step.  Aut(Z12) contains opposite step pairs (1,11) and (5,7); the current is antisymmetric under u -> -u, so the Aut-averaged selector-free current is identically zero.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_post_p2754_entropy_current_obligation": scan["all_patterns_have_hits"],
        "nonzero_directed_entropy_current_exists": audit["has_nonzero_directed_entropy_currents"],
        "opposite_step_pairing_verified": audit["all_opposite_steps_cancel"],
        "aut_averaged_entropy_current_identically_zero": audit["aut_averaged_current_identically_zero"],
        "strict_step_or_polarity_selector_exported": False,
        "p2721_entropy_current_coupling_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_strict_entropy_current_selector": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "Directed entropy currents can be nonzero, but the nonzero sign requires a chosen directed step/polarity.  Aut(Z12)-symmetric handling pairs each step with its opposite and cancels the current, and no P2721 coupling theorem is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["entropy_gradient_current_audit"]
    lines = [
        "# P2755/S1705 entropy-gradient current Aut-cancellation audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Directed current scan",
        f"- modulus={a['modulus']}",
        f"- quanta={a['quanta']}",
        f"- composition_count={a['composition_count']}",
        f"- has_nonzero_directed_entropy_currents={a['has_nonzero_directed_entropy_currents']}",
        f"- opposite_pair_failure_count={a['opposite_pair_failure_count']}",
        f"- aut_average_failure_count={a['aut_average_failure_count']}",
        f"- aut_averaged_current_identically_zero={a['aut_averaged_current_identically_zero']}",
        "",
        "## Step rows",
    ]
    for row in a["directed_step_rows"]:
        lines.append(f"- step={row['step']} opposite={row['opposite_step']}: nonzero={row['nonzero_current_count']}, +={row['positive_current_count']}, -={row['negative_current_count']}")
    lines.extend(["", "## Obstruction", a["proof_obstruction"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2754 = read_json(P2754)
    scan = evidence_scan()
    audit = current_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2755_ENTROPY_GRADIENT_CURRENT_AUT_CANCELLATION_AUDIT_NO_GO",
        "input_hashes": {"P2754_SHANNON_ENTROPY_FOUR_BIT_SELECTOR_AUDIT_NO_GO": sha(P2754)},
        "input_statuses": {"P2754_SHANNON_ENTROPY_FOUR_BIT_SELECTOR_AUDIT_NO_GO": p2754.get("status")},
        "audited_candidate_class": "inversion-odd Shannon entropy gradient/current/flux on Z12",
        "content_evidence_scan": scan,
        "entropy_gradient_current_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not replay scalar entropy, but also do not claim the first directed entropy current as closure.  P2755 constructs a genuinely nonzero directed Shannon entropy current, then proves computationally that selector-free Aut(Z12) handling pairs steps 1/11 and 5/7 and cancels the current identically.  The next proof-grade move must either export a strict law selecting a directed entropy-current step/polarity together with an explicit P2721 coupling theorem, or pivot to a different typed object outside scalar entropy and Aut-averaged entropy-current replay; otherwise preserve the P2697-P2755 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2755/S1705 entropy-gradient current Aut-cancellation audit", "## P2755/S1705 entropy-gradient current Aut-cancellation audit\n\n`P2755/S1705` tests the next entropy move after P2754: an inversion-odd Shannon entropy gradient/current/flux.  The concrete directed current `J_u(h)=sum_i h_i h_{i+u}(h_i-h_{i+u})` is nonzero on the finite `Z12` four-quanta entropy-density scan, so the candidate is more dynamical than scalar entropy.  However, the sign requires choosing a directed unit step; `Aut(Z12)` pairs `1` with `11` and `5` with `7`, and the finite audit finds zero opposite-pair failures and zero Aut-average failures.  Thus selector-free Aut handling cancels the current, and no `P2721` coupling theorem, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, bridge closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2755/S1705 entropy-current Ltotal guard", "## P2755/S1705 entropy-current Ltotal guard\n\n`P2755/S1705` constructs a nonzero directed entropy-current toy term but does not add it as a variational source: without a strict law selecting a directed step/polarity and an explicit `P2721` coupling theorem, Aut-symmetric handling cancels the current.  Therefore it does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current entropy-gradient current Aut-cancellation guardrail (P2755/S1705, 2026-06-15)", "## Current entropy-gradient current Aut-cancellation guardrail (P2755/S1705, 2026-06-15)\n\n- P2755 follows the P2754 recommendation by testing a concrete inversion-odd Shannon entropy gradient/current/flux rather than replaying scalar entropy.\n- The directed entropy current is genuinely nonzero on the finite `Z12` four-quanta entropy-density scan, but its sign requires a chosen directed unit step; `Aut(Z12)` pairs steps `1/11` and `5/7`, giving zero opposite-pair failures and identically zero Aut-averaged current.\n- Do not promote a chosen entropy-current direction to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without a strict law selecting the step/polarity and an explicit `P2721` coupling theorem.  Otherwise pivot to a different typed object or preserve the P2697-P2755 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
