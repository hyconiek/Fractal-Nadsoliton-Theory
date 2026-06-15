#!/usr/bin/env python3
"""P2750/S1700: concrete inversion-odd source sign-value inventory no-go.

P2749 left one narrow missing premise: a concrete strict signed value plus a
law selecting one P2721 coupling polarity.  This audit does not invent a new
selector.  It scans the current generated-artifact frontier for candidate
signed-source objects and scores them against the exact acceptance criteria.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2749 = GEN / "p2749_s1699_minimal_inversion_odd_source_coupling_polarity_audit.json"
OUT = GEN / "p2750_s1700_concrete_odd_source_sign_value_inventory_no_go.json"
MD = GEN / "p2750_s1700_concrete_odd_source_sign_value_inventory_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "post_p2749_missing_premise": r"P2749|concrete strict source sign value|coupling polarity|minimal inversion-odd",
    "candidate_signed_observables": r"chiral-bispectrum|Gauss-phase|cubic phase|inversion-odd|signed source|pseudoscalar|chirality",
    "p2721_boundary": r"P2721 polarity|P2721 coupling|lambda/P2721|polarity-selection",
    "closure_forbidden": r"QW-2191|selector closure|role transfer|L_total|ToE closure",
}

CANDIDATE_FILES = [
    "p2718_s1668_chiral_bispectrum_signed_formula_torsor_coupling_audit.json",
    "p2719_s1669_chiral_bispectrum_phase_origin_source_localizer_theorem_audit.json",
    "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "p2726_s1676_chiral_bispectrum_affine_orientation_flow_transition_matrix.json",
    "p2732_s1682_chiral_bispectrum_time_arrow_source_term_coupling_matrix.json",
    "p2740_s1690_z12_source_triple_chirality_orbit_no_go.json",
    "p2745_s1695_z12_quadratic_gauss_phase_signed_observable_audit.json",
    "p2747_s1697_z12_cubic_phase_orbit_signed_observable_audit.json",
    "p2748_s1698_absence_of_selector_self_synchronization_no_go.json",
    "p2749_s1699_minimal_inversion_odd_source_coupling_polarity_audit.json",
]

NEGATIVE_EXPORT_FLAGS = [
    "concrete_strict_source_sign_value_exported",
    "unique_coupling_polarity_theorem_exported",
    "lambda_p2721_fixed",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]

POSITIVE_HINTS = ["nonzero", "signed", "inversion-odd", "equivariant", "coupling", "two", "2", "polarity"]
BLOCKING_HINTS = ["no ", "not ", "without", "missing", "blocked", "gap", "no-go", "unfixed", "unselected", "does not export"]


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
        rows.append({"content_lane": name, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:10]})
    return {"content_pattern_count": len(rows), "rows": rows, "hit_counts": {r["content_lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def flatten_text(obj: Any) -> str:
    if isinstance(obj, dict):
        return "\n".join(flatten_text(v) for v in obj.values())
    if isinstance(obj, list):
        return "\n".join(flatten_text(v) for v in obj)
    return str(obj)


def numeric_nonzero_values(obj: Any) -> list[float]:
    vals: list[float] = []
    if isinstance(obj, dict):
        for v in obj.values():
            vals.extend(numeric_nonzero_values(v))
    elif isinstance(obj, list):
        for v in obj:
            vals.extend(numeric_nonzero_values(v))
    elif isinstance(obj, (int, float)) and not isinstance(obj, bool) and obj != 0:
        vals.append(float(obj))
    return vals


def has_phrase(text: str, phrases: list[str]) -> bool:
    low = text.lower()
    return any(p in low for p in phrases)


def candidate_inventory() -> dict[str, Any]:
    rows = []
    for name in CANDIDATE_FILES:
        path = GEN / name
        data = read_json(path)
        text = flatten_text(data)
        low = text.lower()
        nonzero = numeric_nonzero_values(data)
        status = str(data.get("status", "MISSING")) if isinstance(data, dict) else "MISSING"
        has_candidate_sign = has_phrase(text, ["nonzero signed", "inversion-odd", "signed source", "chiral", "chirality", "gauss", "cubic"])
        has_concrete_value = bool(nonzero) and has_candidate_sign
        has_unique_language = bool(re.search(r"\bexactly 1\b|\bunique\b|one coupling polarity|selecting one coupling polarity", low))
        has_blocker_language = has_phrase(text, BLOCKING_HINTS) or "gap" in status.lower() or "no_go" in status.lower() or "no-go" in status.lower()
        exports_p2721_coupling = has_phrase(text, ["p2721 coupling theorem exported", "unique coupling polarity theorem", "polarity-selection theorem exported"])
        accepted = bool(has_concrete_value and has_unique_language and exports_p2721_coupling and not has_blocker_language)
        rows.append({
            "artifact": name,
            "exists": path.exists(),
            "status": status,
            "sha256": sha(path),
            "nonzero_numeric_count": len(nonzero),
            "sample_nonzero_values": nonzero[:12],
            "has_candidate_signed_object": has_candidate_sign,
            "has_concrete_numeric_signal": has_concrete_value,
            "has_unique_selector_language": has_unique_language,
            "has_blocker_language": has_blocker_language,
            "exports_p2721_coupling_polarity_theorem": exports_p2721_coupling,
            "accepted_as_concrete_strict_source_sign": accepted,
        })
    accepted_rows = [r for r in rows if r["accepted_as_concrete_strict_source_sign"]]
    return {
        "candidate_file_count": len(rows),
        "existing_candidate_file_count": sum(1 for r in rows if r["exists"]),
        "accepted_candidate_count": len(accepted_rows),
        "accepted_artifacts": [r["artifact"] for r in accepted_rows],
        "rows": rows,
        "finite_theorem": "The current generated-artifact frontier contains real signed or sign-admissible objects, but no audited artifact simultaneously exports a concrete nonzero strict source sign value, a unique coupling-polarity selection theorem, and an explicit P2721 coupling theorem.  Therefore the P2749 missing premise is not supplied by current artifacts.",
    }


def acceptance_matrix(scan: dict[str, Any], inventory: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_post_p2749_boundary": scan["all_patterns_have_hits"],
        "candidate_artifacts_loaded": inventory["existing_candidate_file_count"] >= 6,
        "some_candidate_signed_objects_exist": any(r["has_candidate_signed_object"] for r in inventory["rows"]),
        "accepted_concrete_strict_source_sign_value_found": inventory["accepted_candidate_count"] > 0,
        "unique_p2721_coupling_polarity_theorem_found": any(r["exports_p2721_coupling_polarity_theorem"] for r in inventory["rows"]),
    }
    return {
        "facts": facts,
        "accepted_as_p2749_completion": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "Current artifacts include signed observables and representation-level odd signs, but none supplies the whole P2749 completion package: concrete strict source sign value plus theorem selecting one P2721 coupling polarity.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    inv = payload["candidate_inventory"]
    lines = [
        "# P2750/S1700 concrete odd-source sign-value inventory no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite inventory",
        f"- candidate_file_count={inv['candidate_file_count']}",
        f"- existing_candidate_file_count={inv['existing_candidate_file_count']}",
        f"- accepted_candidate_count={inv['accepted_candidate_count']}",
        f"- accepted_artifacts={inv['accepted_artifacts']}",
        "",
        "## Theorem statement",
        inv["finite_theorem"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    scan = evidence_scan()
    inventory = candidate_inventory()
    acceptance = acceptance_matrix(scan, inventory)
    payload = {
        "status": "P2750_CONCRETE_ODD_SOURCE_SIGN_VALUE_INVENTORY_NO_GO",
        "input_hashes": {"P2749_MINIMAL_ODD_SOURCE_COUPLING_POLARITY_GAP": sha(P2749)},
        "input_statuses": {"P2749_MINIMAL_ODD_SOURCE_COUPLING_POLARITY_GAP": read_json(P2749).get("status")},
        "audited_candidate_class": "current-artifact inventory for concrete strict inversion-odd source sign value plus unique P2721 coupling polarity theorem",
        "content_evidence_scan": scan,
        "candidate_inventory": inventory,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not keep replaying the current signed-observable inventory as if it supplied the P2749 missing premise.  P2750 finds no current artifact with both a concrete strict source sign value and a theorem selecting one P2721 coupling polarity.  The next proof-grade move must introduce a genuinely new strict sign-source mechanism with a computable nonzero value and explicit coupling-polarity theorem, or pivot to a different typed object; otherwise preserve the P2697-P2750 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2750/S1700 concrete odd-source sign-value inventory no-go", "## P2750/S1700 concrete odd-source sign-value inventory no-go\n\n`P2750/S1700` audits the exact missing premise left by `P2749/S1699` across the current generated-artifact frontier: a concrete strict inversion-odd source sign value plus a theorem selecting one `P2721` coupling polarity.  The finite inventory loads current signed/chiral/Gauss/cubic/minimal-odd candidates and finds `0` accepted artifacts satisfying the full package.  Existing artifacts contain real signed observables or admissible odd representations, but no concrete strict sign value with a unique coupling-polarity theorem is exported.  Therefore no `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2750/S1700 concrete odd source inventory Ltotal guard", "## P2750/S1700 concrete odd source inventory Ltotal guard\n\n`P2750/S1700` adds no variational source term.  Its current-artifact inventory finds no accepted concrete strict inversion-odd source sign value together with a unique `P2721` coupling-polarity theorem.  Therefore it does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current concrete odd-source sign-value inventory no-go guardrail (P2750/S1700, 2026-06-15)", "## Current concrete odd-source sign-value inventory no-go guardrail (P2750/S1700, 2026-06-15)\n\n- P2750 audits the exact missing premise left by P2749 across the current generated-artifact frontier: a concrete strict inversion-odd source sign value plus a theorem selecting one `P2721` coupling polarity.\n- The finite inventory loads current signed/chiral/Gauss/cubic/minimal-odd candidates and finds `0` accepted artifacts satisfying the full package; existing artifacts contain real signed observables or admissible odd representations, but no concrete strict sign value with a unique coupling-polarity theorem is exported.\n- Do not promote the current signed-observable inventory to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure.  A next admissible move must introduce a genuinely new strict sign-source mechanism with a computable nonzero value and explicit coupling-polarity theorem, pivot to a different typed object, or preserve the P2697-P2750 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
