#!/usr/bin/env python3
"""P1969 S919 strict ToE progress dashboard after selector-axis chain.

This packet is intentionally a synthesis/checkpoint, not a theorem closure.  It
is the next honest move after grepping the repository: P1966-P1968 made real
progress on QW-2191 by separating kernel-only obstruction, axis-only strict
sources, and residual-sign safety for projective OS observables.  The question
now is what that progress means for ToE closure and what remains the real
strict bottleneck.

The answer exported here is: projective/quadratic OS continuation is safe, but
strict-core/global ToE closure remains open because full PO2 EOM extraction,
dressed Cutkosky unitarity, and sign-sensitive/global selector closure are not
closed.
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1969_s919_strict_toe_progress_dashboard_after_selector_axis_chain.json"

GREP_FILES = [
    "P706_CURRENT_RELEASE_7_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_CLOSURE_DASHBOARD_PROBE.md",
    "N327_CURRENT_FIRST_STRICT_TOE_CLOSURE_DOMINANT_MISSING_INGREDIENT_CLASS_THEOREM.md",
    "P1954_STRICT_DRESSED_AMPLITUDE_NONAVAILABILITY_THEOREM_PACKET.md",
    "P1965_STRICT_DELTA_BG_YF_EOM_NORMAL_FORM_EXTRACTION_NONAVAILABILITY_PACKET.md",
    "P1968_STRICT_P1967_AXES_RESIDUAL_Z2_DOWNSTREAM_OBSERVABLE_AUDIT_PACKET.md",
    "N707_CURRENT_STRICT_T173_PREVIOUS_METHODOLOGY_SURVIVAL_AND_GLOBAL_GAP_BOUNDARY_THEOREM.md",
]

GREP_TERMS = [
    "ToE closure",
    "Operational ToE",
    "strict-core ToE",
    "QW-2191",
    "directed/sign-sensitive",
    "PO2",
    "EOM",
    "unitarity",
    "dressed amplitude",
    "projective",
    "next honest",
]


def load_generated(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode("utf-8")).hexdigest()


def file_sha(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def grep_evidence() -> list[dict[str, Any]]:
    rows = []
    for name in GREP_FILES:
        path = ROOT / name
        if not path.exists():
            rows.append({"path": name, "present": False, "matches": []})
            continue
        matches = []
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            lower = line.lower()
            if any(term.lower() in lower for term in GREP_TERMS):
                matches.append({"line": lineno, "text": line.strip()[:260]})
            if len(matches) >= 10:
                break
        rows.append({"path": name, "present": True, "sha256": file_sha(path), "matches": matches})
    return rows


def row(block_id: str, status: str, evidence: list[str], contribution: str, blocker: str) -> dict[str, Any]:
    return {
        "block_id": block_id,
        "status": status,
        "evidence": evidence,
        "toe_progress_contribution": contribution,
        "remaining_blocker": blocker,
    }


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)

    p1950 = load_generated("p1950_s900_strict_b1_renormalization_exact_cancellation_probe.json")
    p1954 = load_generated("p1954_s904_strict_dressed_amplitude_nonavailability_theorem.json")
    p1963 = load_generated("p1963_s913_strict_po3_double_run_machine_checker.json")
    p1965 = load_generated("p1965_s915_strict_delta_bg_yf_eom_normal_form_extraction_nonavailability.json")
    p1968 = load_generated("p1968_s918_strict_p1967_axes_residual_z2_downstream_observable_audit.json")
    p706 = load_generated("p706_current_release_7_strict_projective_operational_toe_os_closure_dashboard_probe_summary.json")
    p709 = load_generated("p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe_summary.json")

    status_rows = [
        row(
            "projective_operational_os_lane",
            "PASS_PROJECTIVE_OPERATIONAL_OS_READY",
            [
                f"P706 status={p706.get('status')}",
                f"P1968 release7_projective_os_safe_under_residual_z2={p1968.get('machine_checks', {}).get('release7_projective_os_safe_under_residual_z2')}",
                f"P709 sign_ok={p709.get('sign_ok')} patterns={p709.get('n_patterns_tested')}",
            ],
            "Projective/quadratic OS work can proceed without waiting for sign-sensitive QW-2191 closure.",
            "Does not imply strict-core/global ToE closure or physical-unit/SM identification.",
        ),
        row(
            "selector_axis_only_qw2191_lane",
            "PASS_AXIS_ONLY_O2_TO_Z2__GLOBAL_SIGN_OPEN",
            [
                f"P1968 global_sign_sensitive_selector_closure_pass={p1968.get('machine_checks', {}).get('global_sign_sensitive_selector_closure_pass')}",
                f"P1968 blocked_by_existing_audits={p1968.get('machine_checks', {}).get('global_sign_sensitive_selector_blocked_by_existing_audits')}",
            ],
            "QW-2191 no longer blocks projective/quadratic OS observables after the P1966-P1968 chain.",
            "Directed/sign-sensitive physical orientation and admissible S_sel_int remain open.",
        ),
        row(
            "po3_formal_domain",
            "PASS_MACHINE_CHECKED_FORMAL_NONEMPTY__COVARIANT_GLOBAL_CHECK_OPEN",
            [f"P1963 status={p1963.get('status')}", f"P1963 po3_restamp={p1963.get('po3_restamp', {}).get('after')}"] ,
            "The formal parameter/domain nonemptiness blocker is no longer merely hypothetical.",
            "Full covariant consistency/global background-independence transport remains tied to PO2/EOM export.",
        ),
        row(
            "po2_eom_normal_form",
            "OPEN_FULL_EOM_NORMAL_FORM_EXTRACTION_REQUIRED",
            [
                f"P1965 status={p1965.get('status')}",
                f"P1965 failed_or_partial_requirement_ids={p1965.get('failed_or_partial_requirement_ids')}",
            ],
            "P1964 algebraic sufficiency is retained, but P1965 prevents false promotion.",
            "Need explicit variational EOM normal-form extraction from full L_total; strongest blocker for background/PO2 closure.",
        ),
        row(
            "unitarity_dressed_cutkosky",
            "OPEN_DRESSED_AMPLITUDE_UNDERDETERMINED",
            [f"P1954 local_verdict={p1954.get('local_verdict')}", f"P1954 status={p1954.get('status')}"] ,
            "Seed unitarity work is scoped and audited.",
            "Need dressed amplitude, BRST projector, common-basis discontinuity/cut equality, and dressed pole residues.",
        ),
        row(
            "renormalization_b1",
            "LOCAL_PASS_DECLARED_B1_ANSATZ__GLOBAL_SCOPE_OPEN",
            [f"P1950 local_verdict={p1950.get('local_verdict')}", f"P1950 status={p1950.get('status')}"] ,
            "Declared B1 counterterm cancellation is machine-witnessed locally.",
            "Global/background-universal theorem and backend coefficient provenance remain open.",
        ),
    ]

    pass_like = [r for r in status_rows if r["status"].startswith("PASS") or "LOCAL_PASS" in r["status"]]
    open_like = [r for r in status_rows if r not in pass_like or "OPEN" in r["status"]]

    professor_decision = {
        "decision": "Do not spend the next strict-core step on another selector-axis expansion unless a directed observable is explicitly required.",
        "reason": (
            "P1966-P1968 moved QW-2191 from a projective/quadratic OS blocker to a residual sign-sensitive/global-selector blocker. "
            "The highest leverage ToE move is now the full L_total -> EOM -> DELTA_BG_Yf normal-form extraction, with dressed Cutkosky amplitude as the parallel hard block."
        ),
        "recommended_next_packet": "P1970_STRICT_HIGGS_YUKAWA_CURVATURE_VARIATIONAL_NORMAL_FORM_EXTRACTION_ATTEMPT",
        "acceptance_for_next_packet": [
            "freeze metric-density, covariant derivative, and integration-by-parts conventions",
            "derive Euler-Lagrange rows for the Higgs/Yukawa/nonminimal-curvature subsector",
            "project the result onto the P1930/P1964 DELTA_BG_Yf tensorial basis",
            "prove Omega_unexported = 0 or export the exact leftover term as a nonavailability witness",
        ],
    }

    out = {
        "packet_id": "P1969",
        "stage_id": "S919",
        "status": "STRICT_TOE_PROGRESS_DASHBOARD_AFTER_SELECTOR_AXIS_CHAIN__PROJECTIVE_OS_SAFE__GLOBAL_TOE_OPEN",
        "route": "strict_only_progress_synthesis_no_false_pass",
        "legacy_bridge_used": False,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "repo_grep_evidence": grep_evidence(),
        "input_hashes": {
            "p1950_sha256": digest(p1950),
            "p1954_sha256": digest(p1954),
            "p1963_sha256": digest(p1963),
            "p1965_sha256": digest(p1965),
            "p1968_sha256": digest(p1968),
            "p706_sha256": digest(p706),
            "p709_sha256": digest(p709),
        },
        "toe_progress_matrix": status_rows,
        "progress_counts": {
            "rows_total": len(status_rows),
            "pass_or_local_pass_rows": len(pass_like),
            "open_or_scope_limited_rows": len(open_like),
            "strict_global_toe_closed": False,
            "projective_operational_os_lane_safe": True,
        },
        "professorial_decision": professor_decision,
        "false_pass_guard": [
            "Projective Operational OS safety is not global strict-core ToE closure.",
            "Axis-only O(2)->Z2 is not directed/sign-sensitive physical orientation.",
            "P1964/P1965 prevent promoting conditional PO2 algebra into full EOM sufficiency.",
            "P1954 prevents promoting seed Cutkosky checks into dressed unitarity closure.",
        ],
        "lay_explanation": (
            "Postęp jest realny: dla wielkości projektorowych i kwadratowych problem wyboru osi nie blokuje już pracy operacyjnej. "
            "Ale pełna Teoria Wszystkiego wymaga jeszcze mechanicznego wyprowadzenia równań ruchu i pełnej amplitudy unitarności, a nie tylko bezpiecznych projekcyjnych obserwabli."
        ),
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
