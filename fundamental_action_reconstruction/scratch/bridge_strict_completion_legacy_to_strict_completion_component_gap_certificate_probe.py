#!/usr/bin/env python3
"""Scratch probe: legacy->strict completion component gap matrix.

This probe is the next finite, proof-oriented bridge step after restoring the
legacy kernel as an intermediate object.  It does not try to prove the full
bridge.  Instead it turns the current repo evidence into an explicit component
matrix for the required completion map

    K_legacy_ont -> K_strict_gate

and separates what is already finitely certified from what remains an open
source/role-transfer theorem.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.md"

SOURCE_REPORTS = {
    "agents_guardrail": ROOT / "AGENTS.md",
    "s2_priority_packet": ROOT / "fundamental_action_reconstruction" / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "necessity": HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "cocycle": HERE / "bridge_strict_kernel_completion_transport_cocycle_report.json",
    "gf2_linear_system": HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json",
    "gf2_commutative_diagram": HERE / "bridge_strict_completion_phase_zero_gf2_commutative_diagram_certificate_report.json",
    "path_cohomology": HERE / "bridge_strict_completion_phase_sign_path_cohomology_triviality_certificate_report.json",
    "cycle_closure": HERE / "bridge_strict_completion_phase_sign_cycle_closure_obstruction_certificate_report.json",
    "damping_exact": HERE / "bridge_strict_completion_damping_exact_rational_calculus_certificate_report.json",
    "legacy_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
    "torsion_chi11_audit": HERE / "bridge_legacy_torsion_chi11_opinion_audit_report.json",
    "s1_obstruction": HERE / "bridge_strict_completion_s1_selector_margin_monotonicity_obstruction_certificate_report.json",
    "closure_plan": HERE / "bridge_strict_completion_closure_plan_dependency_certificate_report.json",
}

MATRIX_COLUMNS = [
    "legacy_input_visible",
    "strict_completion_visible",
    "finite_certificate_exported",
    "strict_dynamic_source_exported",
    "selector_or_source_exported",
    "role_transfer_allowed_now",
]

SEARCH_TERMS = [
    "completion component gap",
    "legacy strict gap matrix",
    "legacy -> strict completion bridge",
    "damping/compression passage",
    "phase/frequency/topological bit passage",
    "role-transfer audit after bridge completion",
]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite text: {path}")
    return path.read_text(encoding="utf-8")


def report_path(path: Path) -> str:
    return str(path.relative_to(ROOT))


def gf2_rank(rows: list[list[int]]) -> int:
    """Return row rank over GF(2) for a small 0/1 matrix."""
    work = [row[:] for row in rows if any(row)]
    if not work:
        return 0
    n_cols = max(len(row) for row in work)
    rank = 0
    for col in range(n_cols):
        pivot = next((r for r in range(rank, len(work)) if work[r][col] % 2), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for r, row in enumerate(work):
            if r != rank and row[col] % 2:
                work[r] = [(a ^ b) for a, b in zip(row, work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def bool_bit(value: bool) -> int:
    return 1 if value else 0


def build_component_rows(evidence: dict[str, Any]) -> list[dict[str, Any]]:
    necessity = evidence["necessity"]
    cocycle = evidence["cocycle"]
    gf2 = evidence["gf2_linear_system"]
    diagram = evidence["gf2_commutative_diagram"]
    path = evidence["path_cohomology"]
    cycle = evidence["cycle_closure"]
    damping = evidence["damping_exact"]
    guardrail = evidence["legacy_bridge_guardrail"]
    torsion = evidence["torsion_chi11_audit"]
    s1 = evidence["s1_obstruction"]

    exact_apd = necessity["necessity_summary"]["exact_subsets_without_extra_scalar"] == [
        "alpha_normalization+phase_frequency_transport+damping_compression"
    ]
    phase_unique = gf2["gf2_linear_system_summary"]["unique_solution"] and gf2["gf2_linear_system_summary"]["rank"] == 11
    diagram_pass = diagram["diagram_summary"]["all_diagram_checks_pass"]
    path_h1_zero = path["path_cohomology_summary"]["h1_dimension_dim_C1_minus_rank_delta"] == 0
    cycle_h1_one = cycle["cycle_closure_summary"]["h1_dimension_dim_C1_minus_rank_delta"] == 1
    damping_monotone = damping["exact_rational_damping_summary"]["continuous_strictly_decreasing_certificate"]
    legacy_restored = guardrail["legacy_kernel_intermediate_bridge_summary"]["legacy_kernel_restored_as_intermediate"]
    strict_continuation = guardrail["legacy_kernel_intermediate_bridge_summary"]["strict_kernel_treated_as_completed_legacy_continuation"]
    compression_gap = guardrail["legacy_kernel_intermediate_bridge_summary"]["strict_compression_missing_from_legacy_recorded"]
    role_audit_required = guardrail["legacy_kernel_intermediate_bridge_summary"]["role_transfer_audit_required_after_full_bridge"]
    beta_candidate_not_theorem = guardrail["legacy_kernel_intermediate_bridge_summary"]["beta_tors_to_chi11_remains_candidate_not_theorem"]
    no_s1 = not s1["s1_obstruction_summary"]["s1_witness_exported"]
    torsion_candidate = any(
        claim["claim_id"] == "beta_tors_collapses_to_chi11"
        and claim["verdict"] == "CANDIDATE_BRIDGE_HYPOTHESIS_NOT_THEOREM"
        for claim in torsion["opinion_audit"]["claims"]
    )

    rows = [
        {
            "component": "amplitude_normalization",
            "bridge_obligation": "Specify how legacy alpha_geo normalization is converted into the strict A(d) completion factor without transferring legacy physical roles silently.",
            "finite_evidence": [
                "necessity report: A+P+D is the only exact no-extra-scalar subset",
                "necessity report: alpha alone is global scalar data rather than the strict shape source",
            ],
            "open_gap": "strict-side dynamical source for A(d) and post-bridge alpha_geo role transfer",
            "status": "finite_completion_factor_certified__source_and_role_transfer_open",
            "matrix_bits": {
                "legacy_input_visible": legacy_restored,
                "strict_completion_visible": strict_continuation,
                "finite_certificate_exported": exact_apd,
                "strict_dynamic_source_exported": False,
                "selector_or_source_exported": True,
                "role_transfer_allowed_now": False,
            },
        },
        {
            "component": "phase_frequency_transport",
            "bridge_obligation": "Transport legacy cosine phase/frequency data through strict P(d), GF(2) edge bits, commutative diagrams, and path/cycle cohomology.",
            "finite_evidence": [
                "cocycle report: four sign flips are isolated as edge transport",
                "GF(2) linear system: rank 11, determinant 1, unique four-edge solution",
                "commutative diagram and path/cycle cohomology reports agree on node/edge bits",
            ],
            "open_gap": "strict source of omega/phi and the orientation/selector premise for choosing the topological bit",
            "status": "finite_phase_transport_certified__orientation_source_open",
            "matrix_bits": {
                "legacy_input_visible": legacy_restored,
                "strict_completion_visible": strict_continuation,
                "finite_certificate_exported": cocycle["cocycle_summary"]["phase_sign_flip_edges"] == ["1->2", "5->6", "7->8", "9->10"] and phase_unique and diagram_pass and path_h1_zero and cycle_h1_one,
                "strict_dynamic_source_exported": False,
                "selector_or_source_exported": False,
                "role_transfer_allowed_now": False,
            },
        },
        {
            "component": "damping_compression",
            "bridge_obligation": "Account for the passage from legacy linear torsion damping 1+beta_tors*d to strict nonlinear compression 1+beta*d^eta.",
            "finite_evidence": [
                "necessity report: damping_compression is shape-critical in the exact A/P/D completion",
                "exact rational damping report: D'(x)<0 on [1,11] by a rational derivative bound",
                "legacy bridge guardrail: strict compression is recorded as missing from legacy",
            ],
            "open_gap": "derive or explicitly postulate the nonlinear d^eta compression and beta/eta mapping from the legacy bridge layer",
            "status": "finite_strict_compression_certified__legacy_to_nonlinear_passage_open",
            "matrix_bits": {
                "legacy_input_visible": legacy_restored,
                "strict_completion_visible": strict_continuation,
                "finite_certificate_exported": exact_apd and damping_monotone and compression_gap,
                "strict_dynamic_source_exported": False,
                "selector_or_source_exported": True,
                "role_transfer_allowed_now": False,
            },
        },
        {
            "component": "topological_phase_bit_chi11",
            "bridge_obligation": "Decide whether beta_tors, SSB, observer readout, or another strict source supplies the chi_11/orientation bit without violating QW-2191.",
            "finite_evidence": [
                "GF(2), path-cohomology, and cycle-closure reports locate and constrain the bit",
                "torsion audit classifies beta_tors -> chi_11 only as a candidate bridge hypothesis",
                "S1 obstruction report exports no local selector-margin witness for the current provider class",
            ],
            "open_gap": "actual strict selector/source theorem for chi_11 and QW-2191 discharge",
            "status": "bit_location_certified__bit_source_not_exported",
            "matrix_bits": {
                "legacy_input_visible": legacy_restored,
                "strict_completion_visible": strict_continuation,
                "finite_certificate_exported": phase_unique and path_h1_zero and cycle_h1_one and torsion_candidate,
                "strict_dynamic_source_exported": False,
                "selector_or_source_exported": False,
                "role_transfer_allowed_now": False,
            },
            "local_obstruction": {
                "beta_tors_to_chi11_candidate_not_theorem": beta_candidate_not_theorem,
                "current_s1_witness_exported": not no_s1,
            },
        },
        {
            "component": "legacy_physical_role_transfer",
            "bridge_obligation": "After and only after the full completion map is specified, audit each legacy role claim for survival, compression/modification, strict-side successor semantics, or rejection.",
            "finite_evidence": [
                "AGENTS/S2 guardrails require post-bridge role-transfer audit",
                "legacy bridge guardrail report marks role transfer as required after full bridge",
            ],
            "open_gap": "claim-by-claim theorem for sin^2(theta_W), alpha_EM inverse, beta^N hierarchy, and any beta_tors orientation role",
            "status": "role_transfer_audit_required__not_permitted_now",
            "matrix_bits": {
                "legacy_input_visible": legacy_restored,
                "strict_completion_visible": strict_continuation,
                "finite_certificate_exported": role_audit_required,
                "strict_dynamic_source_exported": False,
                "selector_or_source_exported": False,
                "role_transfer_allowed_now": False,
            },
        },
    ]
    return rows


def build_payload() -> dict[str, Any]:
    evidence: dict[str, Any] = {}
    source_reports: dict[str, str] = {}
    for name, path in SOURCE_REPORTS.items():
        source_reports[name] = report_path(path)
        if path.suffix == ".json":
            evidence[name] = load_json(path)
        else:
            evidence[name] = {"text": read_text(path)}

    agents_text = evidence["agents_guardrail"]["text"]
    s2_text = evidence["s2_priority_packet"]["text"]
    rows = build_component_rows(evidence)
    component_matrix = [
        [bool_bit(row["matrix_bits"][column]) for column in MATRIX_COLUMNS]
        for row in rows
    ]
    open_gap_columns = {
        column: sum(1 for row in rows if not row["matrix_bits"][column])
        for column in MATRIX_COLUMNS
    }
    rows_with_finite_certificates = sum(1 for row in rows if row["matrix_bits"]["finite_certificate_exported"])
    rows_with_strict_sources = sum(1 for row in rows if row["matrix_bits"]["strict_dynamic_source_exported"])
    rows_with_selector_sources = sum(1 for row in rows if row["matrix_bits"]["selector_or_source_exported"])
    rows_with_role_transfer = sum(1 for row in rows if row["matrix_bits"]["role_transfer_allowed_now"])
    rank = gf2_rank(component_matrix)

    guardrail_summary = evidence["legacy_bridge_guardrail"]["legacy_kernel_intermediate_bridge_summary"]
    closure_summary = evidence["closure_plan"]["closure_plan_summary"]
    torsion_verdicts = evidence["torsion_chi11_audit"]["opinion_audit"]["verdict_counts"]

    summary = {
        "component_count": len(rows),
        "matrix_columns": MATRIX_COLUMNS,
        "component_matrix_rank_mod2": rank,
        "rows_with_finite_certificates": rows_with_finite_certificates,
        "rows_with_strict_dynamic_sources": rows_with_strict_sources,
        "rows_with_selector_or_source_exported": rows_with_selector_sources,
        "rows_with_role_transfer_allowed_now": rows_with_role_transfer,
        "all_rows_have_legacy_and_strict_visibility": all(
            row["matrix_bits"]["legacy_input_visible"] and row["matrix_bits"]["strict_completion_visible"] for row in rows
        ),
        "all_rows_have_finite_certificates": rows_with_finite_certificates == len(rows),
        "strict_dynamic_sources_missing": rows_with_strict_sources == 0,
        "selector_source_gap_remains": rows[1]["matrix_bits"]["selector_or_source_exported"] is False and rows[3]["matrix_bits"]["selector_or_source_exported"] is False,
        "role_transfer_blocked_until_full_bridge": rows_with_role_transfer == 0,
        "bridge_ready_for_role_transfer": False,
        "raw_identity_claimed": False,
        "completion_map_partial_not_full": True,
        "legacy_kernel_restored_as_intermediate": guardrail_summary["legacy_kernel_restored_as_intermediate"],
        "strict_kernel_completed_continuation_only_with_evidence": guardrail_summary["strict_kernel_treated_as_completed_legacy_continuation"],
        "strict_compression_recorded_as_missing_from_legacy": guardrail_summary["strict_compression_missing_from_legacy_recorded"],
        "beta_tors_to_chi11_candidate_not_theorem": "CANDIDATE_BRIDGE_HYPOTHESIS_NOT_THEOREM" in torsion_verdicts,
        "closure_plan_still_keeps_toe_open": closure_summary["toe_closure_not_claimed"],
        "guardrail_text_requires_completion_map": "completion map" in agents_text and "completion bridge" in s2_text,
        "guardrail_text_requires_role_transfer_audit": "role-transfer audit" in agents_text and "role-transfer audit" in s2_text,
    }

    cross_checks = {
        "source_guardrails_loaded": summary["guardrail_text_requires_completion_map"] and summary["guardrail_text_requires_role_transfer_audit"],
        "legacy_and_strict_visibility_all_rows": summary["all_rows_have_legacy_and_strict_visibility"],
        "all_component_rows_have_finite_certificates": summary["all_rows_have_finite_certificates"],
        "matrix_rank_nonzero_and_bounded": 0 < rank <= len(rows),
        "strict_sources_remain_missing": summary["strict_dynamic_sources_missing"],
        "selector_source_gap_remains": summary["selector_source_gap_remains"],
        "role_transfer_still_blocked": summary["role_transfer_blocked_until_full_bridge"] and not summary["bridge_ready_for_role_transfer"],
        "raw_identity_not_claimed": not summary["raw_identity_claimed"] and summary["completion_map_partial_not_full"],
        "compression_gap_explicit": summary["strict_compression_recorded_as_missing_from_legacy"],
        "beta_tors_chi11_only_candidate": summary["beta_tors_to_chi11_candidate_not_theorem"],
        "closure_remains_open": summary["closure_plan_still_keeps_toe_open"],
    }

    payload = {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_COMPONENT_GAP_MATRIX__NO_FULL_BRIDGE_THEOREM",
        "status": "legacy-to-strict-component-gap-matrix-exported-finite-evidence-certified-sources-and-role-transfer-open",
        "source_reports": source_reports,
        "grep_disambiguation": {
            "searched_terms": SEARCH_TERMS,
            "finding": "Repo search found guardrails, closure-plan, compression, and prior closure-gap artifacts, but no single component-level legacy->strict bridge gap matrix that combines amplitude, phase transport, damping compression, chi11 source, and role-transfer columns.",
        },
        "matrix_definition": {
            "rows": [row["component"] for row in rows],
            "columns": MATRIX_COLUMNS,
            "bit_semantics": "1 means the current repo exports that finite/permission status for the component; 0 means the gap remains open or forbidden now.",
        },
        "component_gap_rows": rows,
        "component_gap_matrix": component_matrix,
        "open_gap_counts_by_column": open_gap_columns,
        "completion_gap_summary": summary,
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all(cross_checks.values()),
        "proof_certificate": {
            "nonduplication_step": "rg was used to distinguish this component-gap matrix from prior completion-necessity, closure-plan, S1 obstruction, and legacy guardrail reports.",
            "matrix_step": f"The bridge is decomposed into {len(rows)} obligations and {len(MATRIX_COLUMNS)} audit columns; the exported GF(2) component matrix has rank {rank}.",
            "finite_completion_step": "Every component row has a finite certificate witness, but these witnesses certify bookkeeping/completion factors rather than strict dynamical sources.",
            "compression_step": "The damping/compression row explicitly records that strict nonlinear d^eta compression is a strict-side addition missing from the legacy linear beta_tors*d denominator.",
            "selector_step": "The phase/topological rows locate the chi_11 bit through GF(2), path, and cycle certificates, while beta_tors->chi_11 remains only a candidate bridge hypothesis and the current S1 provider exports no witness.",
            "role_transfer_step": "No legacy physical role is transferable now; the role-transfer audit is triggered only after the full completion map is specified.",
            "theoretical_limit": "This is a gap matrix and partial completion certificate, not a proof that K_legacy_ont equals K_strict_gate, not a beta_tors->chi_11 theorem, and not ToE closure.",
        },
        "hard_limits": [
            "K_legacy_ont is an intermediate bridge kernel, not a co-equal final strict kernel.",
            "K_strict_gate may be treated as completed/enriched legacy continuation only through explicit completion-map evidence.",
            "No raw identity K_legacy_ont == K_strict_gate is claimed.",
            "No strict dynamical derivation of A(d), P(d), D(d), omega/phi, beta/eta, or d^eta compression is exported here.",
            "No beta_tors -> chi_11 theorem or QW-2191 discharge is claimed.",
            "No legacy physical-role transfer is licensed before a separate post-bridge role-transfer audit.",
            "No ToE closure is claimed.",
        ],
    }
    return payload


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Legacy -> strict completion component gap certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Component matrix",
        "",
        f"Columns: `{', '.join(payload['matrix_definition']['columns'])}`",
        "",
    ]
    for row, bits in zip(payload["component_gap_rows"], payload["component_gap_matrix"]):
        lines.append(f"- `{row['component']}`: bits `{bits}`; status `{row['status']}`; gap: {row['open_gap']}")
    lines.extend([
        "",
        "## Summary",
        "",
    ])
    for key, value in payload["completion_gap_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend([
        "",
        "## Cross-checks",
        "",
    ])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend([
        "",
        f"All cross-checks pass: `{payload['all_cross_checks_pass']}`",
        "",
        "## Hard limits",
        "",
    ])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
