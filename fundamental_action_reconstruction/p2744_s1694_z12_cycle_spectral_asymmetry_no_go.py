#!/usr/bin/env python3
"""P2744/S1694: Z12 cycle spectral-asymmetry signed observable no-go.

P2743 says the next admissible move must either source a transition-unit
character polarity or pivot outside finite Z12 sign-character/frame observables.
This audit takes that pivot: a Hermitian spectral-asymmetry observable for the
Z12 cycle derivative across all exported integer character twists.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2744_s1694_z12_cycle_spectral_asymmetry_no_go.json"
MD = GEN / "p2744_s1694_z12_cycle_spectral_asymmetry_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {"P2743_FRAME_CHARACTER_NO_GO": GEN / "p2743_s1693_affine_frame_transition_character_no_go.json"}

CONTENT_PATTERNS = {
    "post_p2743_pivot_obligation": r"P2743|pivot outside finite `?Z12`? sign-character/frame|strict transition-unit source|no-new-live-frontier",
    "spectral_asymmetry_boundary": r"spectral asymmetry|eta invariant|eta/spectral asymmetry|spectral signed|eigenvalue",
    "cycle_derivative_boundary": r"Z12 cycle|cycle derivative|Hermitian|skew adjacency|character twist",
    "p2721_closure_boundary": r"P2721 polarity|lambda/P2721|QW-2191|selector closure|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "spectral_asymmetry_exported_as_strict_signed_source",
    "nonzero_eta_value_exported",
    "twist_sector_selected_nonpremise",
    "p2721_polarity_selected",
    "lambda_fixed",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]
N = 12
TOL = 1e-10


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


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
    return {
        "content_pattern_count": len(CONTENT_PATTERNS),
        "rows": rows,
        "hit_counts": {row["content_lane"]: row["hit_count"] for row in rows},
        "all_patterns_have_hits": all(row["hit_count"] > 0 for row in rows),
    }


def eigenvalue(twist: int, k: int) -> float:
    # Hermitian form of the oriented nearest-neighbour derivative: -i(T-T^-1).
    # Integer twists are the exported Z12 character sectors; they permute k.
    return 2.0 * math.sin(2.0 * math.pi * ((k + twist) % N) / N)


def sign(value: float) -> int:
    if value > TOL:
        return 1
    if value < -TOL:
        return -1
    return 0


def spectral_asymmetry_audit() -> dict[str, Any]:
    rows = []
    pairing_failures = []
    for twist in range(N):
        values = [eigenvalue(twist, k) for k in range(N)]
        signs = [sign(value) for value in values]
        pair_ok = True
        for k in range(N):
            paired = (-2 * twist - k) % N
            if sign(values[k] + values[paired]) != 0:
                pair_ok = False
                pairing_failures.append({"twist": twist, "k": k, "paired": paired, "sum": values[k] + values[paired]})
        rows.append(
            {
                "twist_sector": twist,
                "positive_eigenvalues": signs.count(1),
                "negative_eigenvalues": signs.count(-1),
                "zero_eigenvalues": signs.count(0),
                "eta_sign_sum": sum(signs),
                "pairing_involution_k_to_minus_2twist_minus_k_passes": pair_ok,
                "rounded_eigenvalues": [round(value, 12) for value in values],
            }
        )
    return {
        "typed_candidate": "Hermitian spectral asymmetry eta-sign-sum of the oriented Z12 cycle derivative over integer character twists",
        "cycle_size": N,
        "twist_sector_count": N,
        "sectors_with_nonzero_eta_sign_sum": sum(1 for row in rows if row["eta_sign_sum"] != 0),
        "sectors_with_balanced_positive_negative_counts": sum(1 for row in rows if row["positive_eigenvalues"] == row["negative_eigenvalues"]),
        "common_positive_count": sorted({row["positive_eigenvalues"] for row in rows}),
        "common_negative_count": sorted({row["negative_eigenvalues"] for row in rows}),
        "common_zero_count": sorted({row["zero_eigenvalues"] for row in rows}),
        "pairing_failure_count": len(pairing_failures),
        "twist_rows": rows,
        "pairing_failures": pairing_failures[:8],
        "finite_theorem": "For every exported integer Z12 character twist, the Hermitian cycle-derivative spectrum is paired by k -> -2*twist-k mod 12.  Each sector has five positive, five negative, and two zero eigenvalues, so the eta sign-sum is zero in all 12 sectors.  Thus this spectral-asymmetry pivot supplies no nonzero strict signed value and cannot fix lambda/P2721 without importing a non-exported twist/source premise.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_spectral_pivot": scan["all_patterns_have_hits"],
        "candidate_spectral_family_computed": audit["twist_sector_count"] == 12,
        "nonzero_eta_value_exported": audit["sectors_with_nonzero_eta_sign_sum"] > 0,
        "strict_twist_or_spectral_source_exported": False,
        "p2721_polarity_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_strict_signed_source": all(facts.values()),
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The Z12 cycle derivative has zero eta sign-sum in every exported integer twist sector, and current artifacts export no strict twist/source premise or P2721 coupling theorem.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["spectral_asymmetry_audit"]
    lines = [
        "# P2744/S1694 Z12 cycle spectral-asymmetry no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite spectral audit",
        f"- cycle_size={audit['cycle_size']}",
        f"- twist_sector_count={audit['twist_sector_count']}",
        f"- sectors_with_nonzero_eta_sign_sum={audit['sectors_with_nonzero_eta_sign_sum']}",
        f"- sectors_with_balanced_positive_negative_counts={audit['sectors_with_balanced_positive_negative_counts']}",
        f"- common_positive_count={audit['common_positive_count']}",
        f"- common_negative_count={audit['common_negative_count']}",
        f"- common_zero_count={audit['common_zero_count']}",
        f"- pairing_failure_count={audit['pairing_failure_count']}",
        "",
        "## Theorem statement",
        audit["finite_theorem"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    scan = evidence_scan()
    audit = spectral_asymmetry_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2744_Z12_CYCLE_SPECTRAL_ASYMMETRY_NO_GO" if not acceptance["accepted_as_strict_signed_source"] else "P2744_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "Z12 cycle spectral-asymmetry signed observable after P2743 frame-character no-go",
        "content_evidence_scan": scan,
        "spectral_asymmetry_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote the Z12 cycle spectral-asymmetry pivot: P2744 proves all 12 exported integer twist sectors have eta sign-sum zero.  The next proof-grade move must either introduce a genuinely new strict source that creates a nonzero spectral asymmetry and couples it to P2721, or pivot outside finite Z12 sign-character/frame/cycle-spectrum observables; otherwise preserve the P2697-P2744 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2744/S1694 Z12 cycle spectral-asymmetry no-go", "## P2744/S1694 Z12 cycle spectral-asymmetry no-go\n\n`P2744/S1694` pivots outside finite `Z12` sign-character/frame observables to the Hermitian spectral asymmetry of the oriented `Z12` cycle derivative across all `12` exported integer character twists.  In every sector the spectrum pairs by `k -> -2*twist-k mod 12`, giving `5` positive, `5` negative, and `2` zero eigenvalues; the eta sign-sum is `0` in all sectors.  No nonzero spectral signed value, strict twist/source theorem, `P2721` coupling theorem, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2744/S1694 spectral-asymmetry Ltotal guard", "## P2744/S1694 spectral-asymmetry Ltotal guard\n\n`P2744/S1694` adds no variational source term: the finite `Z12` cycle spectral-asymmetry audit has eta sign-sum `0` in every exported integer twist sector and lacks a strict twist/source theorem plus `P2721` coupling.  Therefore this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current Z12 cycle spectral-asymmetry no-go guardrail (P2744/S1694, 2026-06-14)", "## Current Z12 cycle spectral-asymmetry no-go guardrail (P2744/S1694, 2026-06-14)\n\n- P2744 pivots outside finite `Z12` sign-character/frame observables to the Hermitian spectral asymmetry of the oriented `Z12` cycle derivative across all `12` exported integer character twists.\n- The finite computation finds `5` positive, `5` negative, and `2` zero eigenvalues in every twist sector, with pairing `k -> -2*twist-k mod 12`; eta sign-sum is `0` in all sectors.\n- Do not promote this spectral-asymmetry pivot to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without a genuinely new strict source creating nonzero spectral asymmetry and a `P2721` coupling theorem.  A next admissible move must supply such a source, pivot outside finite `Z12` sign-character/frame/cycle-spectrum observables, or preserve the P2697-P2744 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
