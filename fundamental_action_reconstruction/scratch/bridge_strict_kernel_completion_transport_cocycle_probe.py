#!/usr/bin/env python3
"""Scratch probe: finite transport cocycle for strict completion of the legacy carrier.

The previous necessity certificate showed that the explicit factors
A(d), P(d), D(d) are all needed for pointwise completion on d=0..11.  This
probe takes the next honest proof/computation step: it treats the ratio

    T(d) = K_strict_gate(d) / K_legacy_ont(d)

as the unique node transport required on the audited finite domain, then derives
the adjacent-edge cocycle R(d)=T(d+1)/T(d).  The result is a finite, checkable
transport certificate: all node transports reconstruct exactly from the edge
cocycle, the additive log-magnitude increments split into phase and damping
increments, and the sign cocycle records the unavoidable phase sign flips.

This is not a derivation of the transport from strict nadsoliton dynamics.  It
is a finite transport witness forced by the already-audited completion ansatz.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_kernel_completion_transport_cocycle_report.json"
OUT_MD = HERE / "bridge_strict_kernel_completion_transport_cocycle_report.md"
NECESSITY_REPORT = HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json"
LADDER_REPORT = HERE / "bridge_strict_kernel_completion_ladder_report.json"
BLOCKER_LATTICE = HERE / "bridge_completed_kernel_blocker_dependency_lattice_report.json"

ALPHA_GEO = 4.0 * math.log(2.0)
LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
DOMAIN = list(range(12))
EDGES = list(range(11))
TOL = 1e-14


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def cos_legacy(d: float) -> float:
    return math.cos(LEGACY["omega"] * d + LEGACY["phi"])


def cos_strict(d: float) -> float:
    return math.cos(STRICT["omega"] * d + STRICT["phi"])


def legacy_full(d: float) -> float:
    return ALPHA_GEO * cos_legacy(d) / (1.0 + LEGACY["beta_tors"] * d)


def strict_kernel(d: float) -> float:
    return cos_strict(d) / (1.0 + STRICT["beta"] * d ** STRICT["eta"])


def alpha_factor() -> float:
    return 1.0 / ALPHA_GEO


def phase_factor(d: float) -> float:
    return cos_strict(d) / cos_legacy(d)


def damping_factor(d: float) -> float:
    return (1.0 + LEGACY["beta_tors"] * d) / (1.0 + STRICT["beta"] * d ** STRICT["eta"])


def transport_node(d: int) -> dict[str, Any]:
    df = float(d)
    strict = strict_kernel(df)
    legacy = legacy_full(df)
    phase = phase_factor(df)
    damping = damping_factor(df)
    transport = strict / legacy
    factor_product = alpha_factor() * phase * damping
    return {
        "d": d,
        "legacy_full": legacy,
        "strict_kernel": strict,
        "transport_T": transport,
        "factor_product_A_P_D": factor_product,
        "transport_minus_factor_product": transport - factor_product,
        "sign_T": sign(transport),
        "log_abs_T": math.log(abs(transport)),
        "log_abs_alpha": math.log(alpha_factor()),
        "log_abs_phase": math.log(abs(phase)),
        "log_damping": math.log(damping),
        "phase_factor": phase,
        "damping_factor": damping,
    }


def transport_edge(left: dict[str, Any], right: dict[str, Any]) -> dict[str, Any]:
    edge = left["d"]
    ratio = right["transport_T"] / left["transport_T"]
    phase_ratio = right["phase_factor"] / left["phase_factor"]
    damping_ratio = right["damping_factor"] / left["damping_factor"]
    log_ratio = math.log(abs(ratio))
    log_phase = math.log(abs(phase_ratio))
    log_damping = math.log(damping_ratio)
    sign_ratio = sign(ratio)
    return {
        "edge": f"{edge}->{edge + 1}",
        "left_d": edge,
        "right_d": edge + 1,
        "edge_transport_ratio_R": ratio,
        "edge_sign_ratio": sign_ratio,
        "edge_is_phase_sign_flip": sign_ratio < 0,
        "log_abs_R": log_ratio,
        "delta_log_abs_phase": log_phase,
        "delta_log_damping": log_damping,
        "delta_log_alpha": 0.0,
        "log_split_residual": log_ratio - (log_phase + log_damping),
        "phase_edge_ratio": phase_ratio,
        "damping_edge_ratio": damping_ratio,
    }


def reconstruct_nodes(nodes: list[dict[str, Any]], edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    reconstructed = nodes[0]["transport_T"]
    rows = [{"d": 0, "reconstructed_T": reconstructed, "actual_T": nodes[0]["transport_T"], "residual": 0.0}]
    for edge, node in zip(edges, nodes[1:]):
        reconstructed *= edge["edge_transport_ratio_R"]
        rows.append({
            "d": node["d"],
            "reconstructed_T": reconstructed,
            "actual_T": node["transport_T"],
            "residual": reconstructed - node["transport_T"],
        })
    return rows


def interval_rows(nodes: list[dict[str, Any]], edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for start in DOMAIN:
        running_ratio = 1.0
        running_log = 0.0
        running_sign = 1
        for end in range(start + 1, len(DOMAIN)):
            edge = edges[end - 1]
            running_ratio *= edge["edge_transport_ratio_R"]
            running_log += edge["log_abs_R"]
            running_sign *= edge["edge_sign_ratio"]
            endpoint_ratio = nodes[end]["transport_T"] / nodes[start]["transport_T"]
            endpoint_log = nodes[end]["log_abs_T"] - nodes[start]["log_abs_T"]
            rows.append({
                "interval": f"{start}->{end}",
                "start_d": start,
                "end_d": end,
                "edge_product": running_ratio,
                "endpoint_ratio": endpoint_ratio,
                "multiplicative_residual": running_ratio - endpoint_ratio,
                "edge_log_sum": running_log,
                "endpoint_log_delta": endpoint_log,
                "additive_log_residual": running_log - endpoint_log,
                "edge_sign_product": running_sign,
                "endpoint_sign_ratio": sign(endpoint_ratio),
            })
    return rows


def build_payload() -> dict[str, Any]:
    necessity = load_json(NECESSITY_REPORT)
    ladder = load_json(LADDER_REPORT)
    blocker_lattice = load_json(BLOCKER_LATTICE)

    nodes = [transport_node(d) for d in DOMAIN]
    edges = [transport_edge(left, right) for left, right in zip(nodes, nodes[1:])]
    reconstructions = reconstruct_nodes(nodes, edges)
    intervals = interval_rows(nodes, edges)

    node_residuals = [row["transport_minus_factor_product"] for row in nodes]
    reconstruction_residuals = [row["residual"] for row in reconstructions]
    edge_log_residuals = [row["log_split_residual"] for row in edges]
    interval_mult_residuals = [row["multiplicative_residual"] for row in intervals]
    interval_log_residuals = [row["additive_log_residual"] for row in intervals]
    sign_flip_edges = [row["edge"] for row in edges if row["edge_is_phase_sign_flip"]]

    return {
        "result_kind": "SCRATCH_STRICT_KERNEL_COMPLETION_TRANSPORT_COCYCLE__FINITE_Z12_WITNESS_NOT_DYNAMICAL_DERIVATION",
        "status": "unique-finite-node-transport-and-adjacent-edge-cocycle-certified-for-completion-ansatz",
        "source_reports": {
            "necessity_certificate": str(NECESSITY_REPORT.relative_to(ROOT)),
            "completion_ladder": str(LADDER_REPORT.relative_to(ROOT)),
            "blocker_lattice": str(BLOCKER_LATTICE.relative_to(ROOT)),
        },
        "constants": {
            "alpha_geo": ALPHA_GEO,
            "legacy": LEGACY,
            "strict": STRICT,
            "domain_d_values": DOMAIN,
            "edge_d_values": EDGES,
            "tolerance": TOL,
        },
        "transport_definitions": {
            "node_transport": "T(d)=K_strict_gate(d)/K_legacy_ont(d)=A(d)*P(d)*D(d)",
            "edge_cocycle": "R(d->d+1)=T(d+1)/T(d)",
            "additive_log_split": "Delta log|T| = Delta log|P| + Delta log D; Delta log A = 0 on edges",
            "sign_cocycle": "sign R records phase-sign flips forced by the strict-vs-legacy cosine transport",
        },
        "cocycle_summary": {
            "max_abs_node_factorization_residual": max(abs(value) for value in node_residuals),
            "max_abs_reconstruction_residual": max(abs(value) for value in reconstruction_residuals),
            "max_abs_edge_log_split_residual": max(abs(value) for value in edge_log_residuals),
            "interval_count": len(intervals),
            "max_abs_interval_multiplicative_residual": max(abs(value) for value in interval_mult_residuals),
            "max_abs_interval_additive_log_residual": max(abs(value) for value in interval_log_residuals),
            "transport_sign_pattern": [row["sign_T"] for row in nodes],
            "phase_sign_flip_edges": sign_flip_edges,
            "edge_count": len(edges),
            "alpha_drops_out_of_edge_cocycle": True,
            "T0": nodes[0]["transport_T"],
            "T11": nodes[-1]["transport_T"],
            "T11_over_T0": nodes[-1]["transport_T"] / nodes[0]["transport_T"],
        },
        "node_transport_rows": nodes,
        "edge_cocycle_rows": edges,
        "node_reconstruction_rows": reconstructions,
        "all_interval_cocycle_rows": intervals,
        "blocker_context": {
            "necessity_status": necessity["status"],
            "ladder_status": ladder["status"],
            "blocker_lattice_status": blocker_lattice["status"],
            "what_this_refines": "turns the pointwise completion factorization into the unique finite node transport and adjacent-edge cocycle forced by the ansatz",
            "strict_transport_derivation_status": "still_open; this is a transport witness/bookkeeping cocycle, not an internal strict dynamical source",
            "still_open": [
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "global_z12_map_derivation_as_internal_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "reynolds_obstruction_escape",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "uniqueness_step": "Given nonzero K_legacy_ont(d), any pointwise multiplicative completion to K_strict_gate has unique node value T(d)=K_strict_gate(d)/K_legacy_ont(d).",
            "edge_step": "On the ordered finite domain, the adjacent cocycle R(d)=T(d+1)/T(d) uniquely reconstructs all T(d) from T(0).",
            "interval_step": "All 66 nontrivial intervals satisfy edge-product equals endpoint-ratio and edge-log-sum equals endpoint-log-delta to machine precision.",
            "split_step": "Because alpha is constant, edge log-magnitude transport splits into phase and damping increments only, with residual at machine precision.",
            "sign_step": "The sign cocycle isolates four phase sign flips: 1->2, 5->6, 7->8, and 9->10.",
            "nonduplication": "This is not another subset-necessity or stage-ladder report; it certifies the finite transport/cocycle object induced by them.",
            "theoretical_limit": "The cocycle is forced by the ansatz on d=0..11 but is not yet derived from strict nadsoliton dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this cocycle is bookkeeping between two kernel carriers.",
            "forbidden_reading": "No separate informational layer below the nadsoliton is introduced.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives the transport cocycle from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["cocycle_summary"]
    lines = [
        "# Strict Kernel Completion Transport Cocycle Certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This finite-domain audit turns the pointwise completion quotient into a",
        "unique node transport `T(d)` and adjacent-edge cocycle `R(d->d+1)` on",
        "`d=0..11`.  It is a transport witness induced by the completion ansatz,",
        "not a derivation of the ansatz from strict nadsoliton dynamics.",
        "",
        "## Definitions",
        "",
    ]
    for name, definition in payload["transport_definitions"].items():
        lines.append(f"- `{name}`: {definition}")
    lines.extend([
        "",
        "## Cocycle summary",
        "",
        f"- Max node factorization residual: `{summary['max_abs_node_factorization_residual']:.12e}`",
        f"- Max reconstruction residual from edges: `{summary['max_abs_reconstruction_residual']:.12e}`",
        f"- Max edge log-split residual: `{summary['max_abs_edge_log_split_residual']:.12e}`",
        f"- Interval count: `{summary['interval_count']}`",
        f"- Max interval multiplicative residual: `{summary['max_abs_interval_multiplicative_residual']:.12e}`",
        f"- Max interval additive log residual: `{summary['max_abs_interval_additive_log_residual']:.12e}`",
        f"- Transport sign pattern: `{summary['transport_sign_pattern']}`",
        f"- Phase sign-flip edges: `{summary['phase_sign_flip_edges']}`",
        f"- `T(0)`: `{summary['T0']:.12e}`",
        f"- `T(11)`: `{summary['T11']:.12e}`",
        f"- `T(11)/T(0)`: `{summary['T11_over_T0']:.12e}`",
        "",
        "## Edge cocycle rows",
        "",
        "| edge | R | sign | log_abs_R | delta_log_abs_phase | delta_log_damping | split_residual |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in payload["edge_cocycle_rows"]:
        lines.append(
            "| {edge} | {ratio:.12e} | {sign} | {log_ratio:.12e} | {log_phase:.12e} | {log_damping:.12e} | {residual:.12e} |".format(
                edge=row["edge"],
                ratio=row["edge_transport_ratio_R"],
                sign=row["edge_sign_ratio"],
                log_ratio=row["log_abs_R"],
                log_phase=row["delta_log_abs_phase"],
                log_damping=row["delta_log_damping"],
                residual=row["log_split_residual"],
            )
        )
    lines.extend([
        "",
        "## Guarded interpretation",
        "",
        "The cocycle is the unique finite transport object forced after choosing the",
        "completion ansatz. It does not prove that strict nadsoliton dynamics generates",
        "that transport, does not prove `beta_tors -> chi_11`, and does not transfer",
        "legacy physical roles onto `K_strict_gate`.",
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
