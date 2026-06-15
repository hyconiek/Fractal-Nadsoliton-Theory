#!/usr/bin/env python3
"""P2746/S1696: Gauss-phase orbit selector/source-law no-go.

P2745 found nonzero affine-orbit signed coefficients for quadratic Gauss phases,
but left the exact missing premise: a strict law selecting one coefficient orbit
and one polarity, with a P2721 coupling theorem.  This audit tests the internal
polarity-blind invariants of those Gauss orbits before any external coefficient
label or sign convention is imported.
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
P2745 = GEN / "p2745_s1695_z12_quadratic_gauss_phase_signed_observable_audit.json"
OUT = GEN / "p2746_s1696_gauss_phase_orbit_selector_source_law_no_go.json"
MD = GEN / "p2746_s1696_gauss_phase_orbit_selector_source_law_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "p2745_missing_premise": r"P2745|Gauss-phase|coefficient orbit|strict law selecting one nonzero",
    "orbit_selector_boundary": r"selector|source law|coefficient-orbit|orbit/polarity|polarity-paired",
    "affine_gauss_boundary": r"affine|quadratic Gauss|orbit-safe signed|signed-sum coefficients",
    "p2721_closure_boundary": r"P2721 polarity|lambda/P2721|QW-2191|selector closure|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "gauss_orbit_selector_exported",
    "gauss_polarity_source_law_exported",
    "p2721_polarity_coupling_exported",
    "lambda_fixed",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


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
    return {"content_pattern_count": len(rows), "rows": rows, "hit_counts": {r["content_lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def gcd_multiset(reps: list[list[int]]) -> list[int]:
    return sorted(math.gcd(a, 12) for a, _b in reps)


def b_residue_multiset(reps: list[list[int]], modulus: int) -> list[int]:
    return sorted(b % modulus for _a, b in reps)


def polarity_blind_signature(orbit: dict[str, Any]) -> tuple[Any, ...]:
    reps = orbit["representatives"]
    return (
        orbit["size"],
        abs(orbit["signed_sum_coefficient"]),
        orbit["positive_rows"] + orbit["negative_rows"],
        orbit["zero_rows"],
        tuple(sorted((orbit["positive_rows"], orbit["negative_rows"]))),
        tuple(gcd_multiset(reps)),
        tuple(b_residue_multiset(reps, 2)),
        tuple(b_residue_multiset(reps, 4)),
    )


def selector_audit(p2745: dict[str, Any]) -> dict[str, Any]:
    nonzero = [o for o in p2745["gauss_phase_audit"]["orbit_rows"] if o["signed_sum_coefficient"] != 0]
    rows = []
    by_signature: dict[str, list[dict[str, Any]]] = {}
    for orbit in nonzero:
        sig = polarity_blind_signature(orbit)
        key = repr(sig)
        row = {
            "orbit_id": orbit["orbit_id"],
            "size": orbit["size"],
            "signed_sum_coefficient": orbit["signed_sum_coefficient"],
            "polarity": 1 if orbit["signed_sum_coefficient"] > 0 else -1,
            "abs_coefficient": abs(orbit["signed_sum_coefficient"]),
            "polarity_blind_signature": key,
            "representatives": orbit["representatives"],
        }
        rows.append(row)
        by_signature.setdefault(key, []).append(row)
    signature_rows = []
    unpaired = []
    for key, members in sorted(by_signature.items()):
        polarities = sorted({m["polarity"] for m in members})
        coeffs = sorted(m["signed_sum_coefficient"] for m in members)
        item = {"signature": key, "member_count": len(members), "orbit_ids": [m["orbit_id"] for m in members], "coefficients": coeffs, "polarities": polarities}
        signature_rows.append(item)
        if polarities != [-1, 1]:
            unpaired.append(item)
    return {
        "typed_candidate": "strict source law selecting one P2745 nonzero quadratic-Gauss affine coefficient orbit and polarity from polarity-blind internal orbit invariants",
        "nonzero_orbit_count": len(nonzero),
        "signature_class_count": len(signature_rows),
        "signature_classes_with_both_polarities": sum(1 for row in signature_rows if row["polarities"] == [-1, 1]),
        "unpaired_signature_class_count": len(unpaired),
        "candidate_unique_selector_count": 0 if not unpaired else len(unpaired),
        "nonzero_orbit_rows": rows,
        "signature_rows": signature_rows,
        "unpaired_signature_rows": unpaired,
        "finite_theorem": "For the 8 nonzero P2745 Gauss-phase affine coefficient orbits, every tested polarity-blind internal signature class contains both positive and negative signed-sum coefficients.  The three signature classes pair coefficients as [-2,-2,2,2], [-1,1], and [-1,1].  Therefore these intrinsic Gauss-orbit data do not select a unique orbit or polarity; a lambda/P2721 source would still need an extra strict sign law not present in current artifacts.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_p2745_missing_premise": scan["all_patterns_have_hits"],
        "nonzero_gauss_orbits_loaded": audit["nonzero_orbit_count"] == 8,
        "all_signature_classes_polarity_paired": audit["unpaired_signature_class_count"] == 0,
        "strict_unique_orbit_selector_exported": False,
        "strict_polarity_source_law_exported": False,
        "p2721_polarity_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_gauss_selector_source": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The P2745 nonzero Gauss coefficients survive affine quotienting, but their polarity-blind internal signatures are exactly polarity-paired; no current artifact exports the extra strict sign law needed to choose one orbit/polarity or couple it to P2721.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["selector_audit"]
    lines = [
        "# P2746/S1696 Gauss-phase orbit selector/source-law no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite selector audit",
        f"- nonzero_orbit_count={a['nonzero_orbit_count']}",
        f"- signature_class_count={a['signature_class_count']}",
        f"- signature_classes_with_both_polarities={a['signature_classes_with_both_polarities']}",
        f"- unpaired_signature_class_count={a['unpaired_signature_class_count']}",
        f"- candidate_unique_selector_count={a['candidate_unique_selector_count']}",
        "",
        "## Signature rows",
    ]
    for row in a["signature_rows"]:
        lines.append(f"- orbit_ids={row['orbit_ids']}; coefficients={row['coefficients']}; polarities={row['polarities']}")
    lines += ["", "## Theorem statement", a["finite_theorem"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2745 = read_json(P2745)
    scan = evidence_scan()
    audit = selector_audit(p2745)
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2746_GAUSS_PHASE_ORBIT_SELECTOR_SOURCE_LAW_NO_GO",
        "input_hashes": {"P2745_GAUSS_PHASE_AUDIT": sha(P2745)},
        "input_statuses": {"P2745_GAUSS_PHASE_AUDIT": p2745.get("status")},
        "audited_candidate_class": "P2745 missing-premise audit: strict Gauss coefficient-orbit/polarity selector source law",
        "content_evidence_scan": scan,
        "selector_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not continue P2745 by choosing a Gauss coefficient orbit or polarity from the current polarity-blind orbit data.  P2746 shows the nonzero Gauss coefficients are real but every intrinsic signature class is paired across both polarities.  The next proof-grade move must either provide a genuinely new strict sign law that breaks one of these paired Gauss signature classes and proves the P2721 coupling, or pivot away from Gauss-phase selectors; otherwise preserve the P2697-P2746 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2746/S1696 Gauss-phase orbit selector/source-law no-go", "## P2746/S1696 Gauss-phase orbit selector/source-law no-go\n\n`P2746/S1696` audits the exact missing premise left by `P2745/S1695`: a strict law selecting one nonzero quadratic-Gauss coefficient orbit and polarity.  On the `8` nonzero affine coefficient orbits, every tested polarity-blind internal signature class contains both signs, pairing coefficients as `[-2,-2,2,2]`, `[-1,1]`, and `[-1,1]`.  Thus current Gauss-orbit data do not export a unique orbit/polarity selector, `P2721` coupling theorem, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2746/S1696 Gauss selector Ltotal guard", "## P2746/S1696 Gauss selector Ltotal guard\n\n`P2746/S1696` adds no variational source term: the nonzero P2745 Gauss coefficients remain paired by polarity-blind orbit signatures and no strict sign law or `P2721` coupling theorem is exported.  Therefore this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current Gauss-phase orbit selector/source-law no-go guardrail (P2746/S1696, 2026-06-14)", "## Current Gauss-phase orbit selector/source-law no-go guardrail (P2746/S1696, 2026-06-14)\n\n- P2746 audits the exact missing premise left by P2745: a strict law selecting one nonzero quadratic-Gauss coefficient orbit and polarity.\n- The finite computation checks the `8` nonzero affine coefficient orbits and finds every polarity-blind internal signature class contains both signs, pairing coefficients as `[-2,-2,2,2]`, `[-1,1]`, and `[-1,1]`; no unique orbit/polarity selector is exported.\n- Do not promote Gauss-phase coefficient choice to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without a genuinely new strict sign law breaking one paired signature class and an explicit `P2721` coupling theorem.  Otherwise pivot away from Gauss-phase selectors or preserve the P2697-P2746 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
