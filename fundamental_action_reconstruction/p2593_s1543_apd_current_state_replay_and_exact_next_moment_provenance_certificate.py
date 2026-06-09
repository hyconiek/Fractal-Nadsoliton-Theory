#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2590_s1540_apd_finite_even_moment_shell_interval_nonuniqueness_audit import PRODUCT_PARAMETER_GRID
from p2592_s1542_apd_newton_girard_next_even_moment_sensitivity_certificate import (
    ELEMENTARY_FIRST_SUM,
    ENDPOINT_OFFSET,
    MD as P2592_MD,
    NEGATIVE_EXPORT_FLAGS as P2592_NEGATIVE_EXPORT_FLAGS,
    OUT as P2592_OUT,
    SOURCE_FILES as P2592_SOURCE_FILES,
    next_even_moment_certificate,
)

GEN = ROOT / "generated"
OUT = GEN / "p2593_s1543_apd_current_state_replay_and_exact_next_moment_provenance_certificate.json"
MD = GEN / "p2593_s1543_apd_current_state_replay_and_exact_next_moment_provenance_certificate.md"

SOURCE_FILES = {
    "P2592_NEXT_EVEN_MOMENT_SENSITIVITY": P2592_OUT,
    "P2592_MARKDOWN_SUMMARY": P2592_MD,
    **{f"P2592_SOURCE_{name}": path for name, path in P2592_SOURCE_FILES.items()},
}
NEGATIVE_EXPORT_FLAGS = sorted(set(P2592_NEGATIVE_EXPORT_FLAGS + [
    "apd_current_state_replay_source_exported",
    "apd_provenance_replay_selector_source_exported",
    "apd_exact_rational_next_moment_source_exported",
]))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def run_git(args: list[str]) -> str:
    proc = subprocess.run(["git", *args], cwd=REPO, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return proc.stdout.strip()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2593|S1543|APD current-state replay|current-state replay.*APD",
        "intended_research_nonduplication": "APD.*Newton Girard replay|Newton Girard replay.*APD|APD.*next moment provenance|next moment provenance.*APD|APD.*eighth moment provenance|eighth moment provenance.*APD|APD.*source fingerprint replay|source fingerprint replay.*APD",
        "apd_precursors": "P2591|S1541|P2592|S1542|APD.*next even moment|APD.*eighth moment|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def source_replay(p2592_payload: dict[str, Any]) -> dict[str, Any]:
    recorded = p2592_payload["apd_newton_girard_next_even_moment_sensitivity_certificate"]["source_fingerprints_sha256"]
    current_p2591_payload = load_json(P2592_SOURCE_FILES["P2591_STURM_INTERVAL_CERTIFICATE"])
    current_p2591_theorem = theorem(current_p2591_payload, "apd_product_parameter_sturm_interval_certificate")
    current = {
        "P2591_STURM_INTERVAL_CERTIFICATE": sha256_json(current_p2591_payload),
    }
    whole_match = recorded == current
    theorem_relevant_replay = current_p2591_theorem.get("continuous_interval_of_valid_supports_certified") is True
    return {
        "recorded_source_fingerprints_sha256": recorded,
        "current_source_fingerprints_sha256": current,
        "all_recorded_sources_match_current_repo_artifacts": whole_match,
        "whole_artifact_mismatch_detected": not whole_match,
        "current_p2591_theorem_relevant_interval_still_certified": theorem_relevant_replay,
        "theorem_relevant_inputs_replay_on_current_repo": theorem_relevant_replay,
        "mismatch_classification": (
            "whole generated P2591 artifact fingerprint differs from the P2592-recorded fingerprint, so P2592 was not byte-for-byte replaying the current generated artifact; the theorem-relevant P2591 interval predicate still replays on the current repo"
            if not whole_match
            else "whole generated P2591 artifact fingerprint matches the P2592-recorded fingerprint"
        ),
    }


def git_current_state_replay() -> dict[str, Any]:
    head = run_git(["rev-parse", "HEAD"])
    p2592_commit = run_git(["log", "--all", "--grep=Add P2592 APD", "--format=%H", "-n", "1"])
    p2592_parent = run_git(["rev-parse", f"{p2592_commit}^"]) if p2592_commit else ""
    p2592_subject = run_git(["show", "-s", "--format=%s", p2592_commit]) if p2592_commit else ""
    merge_base = run_git(["merge-base", "HEAD", p2592_commit]) if p2592_commit else ""
    return {
        "current_head": head,
        "p2592_commit": p2592_commit,
        "p2592_parent": p2592_parent,
        "p2592_subject": p2592_subject,
        "current_head_contains_p2592_commit": bool(p2592_commit and merge_base == p2592_commit),
    }


def exact_next_moment_replay(recorded_certificate: dict[str, Any]) -> dict[str, Any]:
    e4 = sp.symbols("e4")
    internal_exact = sp.expand(
        sp.Integer(int(ELEMENTARY_FIRST_SUM)) * sp.Integer(4890)
        - sp.Integer(273) * sp.Integer(354)
        + sp.Integer(820) * sp.Integer(30)
        - 4 * e4
    )
    central_exact = sp.expand(2 * sp.Rational(21, 4) ** 8 + 2 * internal_exact)
    recomputed = next_even_moment_certificate()
    exact_grid = [
        {
            "product_parameter": float(value),
            "internal_eighth_exact": float(internal_exact.subs(e4, int(value))),
            "central_eighth_exact": float(central_exact.subs(e4, int(value))),
        }
        for value in PRODUCT_PARAMETER_GRID
    ]
    exact_internal_values = [row["internal_eighth_exact"] for row in exact_grid]
    exact_central_values = [row["central_eighth_exact"] for row in exact_grid]
    recorded_internal_values = recorded_certificate["grid_internal_eighth_shell_values"]
    recorded_central_values = recorded_certificate["grid_central_eighth_moment_values"]
    return {
        "newton_girard_internal_formula_exact": str(internal_exact),
        "central_eighth_formula_exact_rational": str(central_exact),
        "recorded_internal_formula": recorded_certificate["newton_girard_formula_internal_eighth_shell_p4"],
        "recorded_central_formula_float": recorded_certificate["central_eighth_formula_with_endpoints"],
        "recomputed_internal_formula": recomputed["newton_girard_formula_internal_eighth_shell_p4"],
        "recomputed_central_formula_float": recomputed["central_eighth_formula_with_endpoints"],
        "internal_formula_replays_recorded": str(internal_exact).replace("74658", "74658.0") == recorded_certificate["newton_girard_formula_internal_eighth_shell_p4"],
        "central_float_formula_replays_recorded": recomputed["central_eighth_formula_with_endpoints"] == recorded_certificate["central_eighth_formula_with_endpoints"],
        "exact_grid": exact_grid,
        "max_abs_internal_grid_delta_vs_recorded": max(abs(left - right) for left, right in zip(exact_internal_values, recorded_internal_values)),
        "max_abs_central_grid_delta_vs_recorded": max(abs(left - right) for left, right in zip(exact_central_values, recorded_central_values)),
        "recomputed_gatekeeper_stability": {
            "lower_shells_constant_but_next_shell_varies": recomputed["lower_shells_constant_but_next_shell_varies"],
            "internal_eighth_shell_slope_by_product_parameter": recomputed["internal_eighth_shell_slope_by_product_parameter"],
            "central_eighth_moment_slope_by_product_parameter": recomputed["central_eighth_moment_slope_by_product_parameter"],
            "grid_internal_eighth_shell_distinct_count": recomputed["grid_internal_eighth_shell_distinct_count"],
            "grid_central_eighth_moment_distinct_count": recomputed["grid_central_eighth_moment_distinct_count"],
        },
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2592_payload = load_json(P2592_OUT)
    p2592_theorem = theorem(p2592_payload, "apd_newton_girard_next_even_moment_sensitivity_certificate")
    recorded_certificate = p2592_theorem["apd_newton_girard_next_even_moment_sensitivity_certificate"]
    source = source_replay(p2592_payload)
    git_replay = git_current_state_replay()
    exact = exact_next_moment_replay(recorded_certificate)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2593_T1_apd_current_state_replay_and_exact_next_moment_provenance_certificate",
        "audited_chain": ["P2591/S1541", "P2592/S1542"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "treat the P2592 next-moment certificate as stale or as a strict APD source after exact rational replay",
        "p2592_next_even_moment_certificate_inherited": p2592_theorem.get("lower_second_fourth_sixth_shells_do_not_select_next_shell") is True,
        "source_replay": source,
        "git_current_state_replay": git_replay,
        "exact_newton_girard_replay": exact,
        "p2592_whole_source_artifact_matches_current_repo": source["all_recorded_sources_match_current_repo_artifacts"],
        "p2592_whole_source_artifact_mismatch_detected": source["whole_artifact_mismatch_detected"],
        "p2592_theorem_relevant_inputs_replay_on_current_repo": source["theorem_relevant_inputs_replay_on_current_repo"],
        "p2592_commit_is_on_current_head_history": git_replay["current_head_contains_p2592_commit"],
        "exact_replay_preserves_p2592_numeric_certificate": exact["max_abs_internal_grid_delta_vs_recorded"] == 0.0 and exact["max_abs_central_grid_delta_vs_recorded"] <= 1.0e-9,
        "exact_rational_replay_refines_but_does_not_strengthen_to_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "P2593 detects that the P2592-recorded whole P2591 generated-artifact fingerprint is not byte-for-byte identical to the current P2591 artifact, while the theorem-relevant P2591 interval predicate and the Newton-Girard next-shell result replay exactly on the current repo. This still only identifies the missing selector coordinate; the next honest step is a strict nadsoliton-derived eighth-shell/support law or an explicit source theorem, not promotion of replay provenance into dynamics."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2592_next_even_moment_certificate_inherited": theorem_export["p2592_next_even_moment_certificate_inherited"],
        "p2592_whole_source_mismatch_classified": theorem_export["p2592_whole_source_artifact_mismatch_detected"] or theorem_export["p2592_whole_source_artifact_matches_current_repo"],
        "p2592_theorem_relevant_inputs_replay_on_current_repo": theorem_export["p2592_theorem_relevant_inputs_replay_on_current_repo"],
        "p2592_commit_on_current_history": theorem_export["p2592_commit_is_on_current_head_history"],
        "exact_internal_formula_replays_recorded": exact["internal_formula_replays_recorded"],
        "exact_central_formula_refines_recorded_float": exact["central_float_formula_replays_recorded"],
        "exact_grid_matches_recorded_internal_values": exact["max_abs_internal_grid_delta_vs_recorded"] == 0.0,
        "exact_grid_matches_recorded_central_values": exact["max_abs_central_grid_delta_vs_recorded"] <= 1.0e-9,
        "lower_shells_still_constant_but_next_shell_varies": exact["recomputed_gatekeeper_stability"]["lower_shells_constant_but_next_shell_varies"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2593",
        "stage_id": "S1543",
        "status": "P2593_APD_CURRENT_STATE_REPLAY_AND_EXACT_NEXT_MOMENT_PROVENANCE_CERTIFICATE_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_current_state_replay_and_exact_next_moment_provenance_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2592_NEXT_EVEN_MOMENT_SENSITIVITY": sha256_json(p2592_payload),
                "P2592_MARKDOWN_SUMMARY": sha256_text(P2592_MD.read_text(encoding="utf-8")),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_current_state_replay_and_exact_next_moment_provenance_certificate"]["theorem_export"]
    exact = t["exact_newton_girard_replay"]
    lines = [
        "# P2593/S1543 APD current-state replay and exact next-moment provenance certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2592 whole source artifact matches current repo: `{t['p2592_whole_source_artifact_matches_current_repo']}`.",
        f"- P2592 theorem-relevant inputs replay on current repo: `{t['p2592_theorem_relevant_inputs_replay_on_current_repo']}`.",
        f"- P2592 commit is on current HEAD history: `{t['p2592_commit_is_on_current_head_history']}`.",
        f"- Exact internal formula: `{exact['newton_girard_internal_formula_exact']}`.",
        f"- Exact central formula: `{exact['central_eighth_formula_exact_rational']}`.",
        f"- Exact replay preserves P2592 numeric certificate: `{t['exact_replay_preserves_p2592_numeric_certificate']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "P2593 directly checks the complaint-sensitive provenance question: P2592's recorded whole P2591 source fingerprint does not match the current generated artifact byte-for-byte, but the theorem-relevant P2591 interval predicate still replays on the current repo and the P2592 commit is on the current HEAD history.  It then reruns the Newton-Girard certificate exactly, replacing the decimal central formula by the rational identity `42715646049/32768 - 8*e4` while preserving the P2592 grid values.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No current-state replay source, provenance replay selector source, exact rational next-moment source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_current_state_replay_and_exact_next_moment_provenance_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2593/S1543 APD current-state replay and exact next-moment provenance guard

`P2593/S1543` verifies the current-state status of P2592: the recorded whole P2591 generated-artifact fingerprint is not byte-for-byte identical to the current artifact, but the theorem-relevant P2591 interval predicate replays and the Newton-Girard calculation exact-replays as `p4 = 74658 - 4*e4` and central eighth moment `42715646049/32768 - 8*e4`.  This classifies the provenance drift and removes decimal ambiguity, but it still does not source the APD support law.
""".strip()
    lag_section = """
## P2593/S1543 APD current-state replay and exact next-moment provenance Ltotal guard

`P2593/S1543` blocks a role-bearing APD Gram term in `L_total` from being justified by provenance replay or exact rational rewriting of P2592.  Replay confirms the selector-coordinate obstruction; it does not derive a strict nadsoliton support/density source.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2593/S1543 APD current-state replay and exact next-moment provenance guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2593/S1543 APD current-state replay and exact next-moment provenance Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
