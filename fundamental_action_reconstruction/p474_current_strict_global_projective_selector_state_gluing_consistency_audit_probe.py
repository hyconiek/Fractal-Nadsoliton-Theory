#!/usr/bin/env python3
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

STATE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
TRANSITION = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
ATLAS = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"

OUT = (
    GENERATED
    / "p474_current_strict_global_projective_selector_state_gluing_consistency_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p474_current_strict_global_projective_selector_state_gluing_consistency_audit_probe_summary.json"
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _p(path: Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def _mat_t(a: list[list[float]]) -> list[list[float]]:
    n = len(a)
    m = len(a[0]) if n else 0
    return [[a[i][j] for i in range(n)] for j in range(m)]


def _mat_mul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    n = len(a)
    k = len(a[0]) if n else 0
    k2 = len(b)
    m = len(b[0]) if k2 else 0
    if k != k2:
        raise ValueError(f"matrix shape mismatch: {n}x{k} times {k2}x{m}")
    out = [[0.0 for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for t in range(k):
            ait = a[i][t]
            if ait == 0.0:
                continue
            bt = b[t]
            for j in range(m):
                out[i][j] += ait * bt[j]
    return out


def _mat_sub_max_abs(a: list[list[float]], b: list[list[float]]) -> float:
    n = len(a)
    m = len(a[0]) if n else 0
    if n != len(b) or (m != (len(b[0]) if b else 0)):
        raise ValueError("matrix shape mismatch in residual computation")
    r = 0.0
    for i in range(n):
        for j in range(m):
            d = abs(a[i][j] - b[i][j])
            if d > r:
                r = d
    return r


def _mat_eye(n: int) -> list[list[float]]:
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def _outer_projector(u: list[float]) -> tuple[list[list[float]], float]:
    norm2 = sum(x * x for x in u)
    if norm2 <= 0.0:
        raise ValueError("nonpositive vector norm")
    inv = 1.0 / norm2
    n = len(u)
    p = [[u[i] * u[j] * inv for j in range(n)] for i in range(n)]
    return p, norm2


@dataclass(frozen=True)
class Edge:
    src: str
    dst: str
    key: str
    operator_ref: str


def _parse_edge_key(edge_key: str) -> tuple[str, str]:
    if "_to_" not in edge_key:
        raise ValueError(f"unexpected transition operator key: {edge_key}")
    src, dst = edge_key.split("_to_", 1)
    return src, dst


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = {"state": STATE, "transition": TRANSITION, "atlas": ATLAS}
    missing_required = [k for k, p in required.items() if not p.is_file()]
    if missing_required:
        payload = {
            "stage": "P474",
            "date": "2026-03-16",
            "status": "FAIL_MISSING_REQUIRED_INPUTS",
            "missing_required_inputs": missing_required,
            "required_paths": {k: _p(v) for k, v in required.items()},
            "no_false_pass": True,
            "hard_limits": [
                "no theorem-level pass",
                "no full-closure pass",
                "no sign-sensitive/directed selector state claim",
                "no claim of global QW-2191 discharge",
            ],
        }
        OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": payload["stage"],
                    "status": payload["status"],
                    "overall_pass": False,
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    state = _read_json(STATE)
    transition = _read_json(TRANSITION)
    atlas = _read_json(ATLAS)

    charts_state = state.get("charts", {})
    if not isinstance(charts_state, dict) or not charts_state:
        raise ValueError("unexpected selector state charts shape")

    chart_ids = sorted(charts_state.keys())

    # Load u_m vectors from the referenced local operator packets (A_m).
    local_operator_refs: dict[str, str] = {}
    u_by_chart: dict[str, list[float]] = {}
    u_norm2_by_chart: dict[str, float] = {}
    missing_local_operator_packets: list[str] = []

    for chart in chart_ids:
        local_op = (charts_state.get(chart) or {}).get("local_operator") or {}
        ref = local_op.get("ref")
        if not isinstance(ref, str):
            raise ValueError(f"missing local_operator.ref for chart {chart}")
        local_operator_refs[chart] = ref
        path = REPO / ref
        if not path.is_file():
            missing_local_operator_packets.append(ref)
            continue
        op_packet = _read_json(path)
        u_key = f"u_{chart.removeprefix('pair')}"
        u = (op_packet.get("data", {}) or {}).get(u_key)
        if not isinstance(u, list) or len(u) != 12:
            raise ValueError(f"missing 12-vector {u_key} in {ref}")
        p_u, norm2 = _outer_projector([float(x) for x in u])
        # Store u and norm; projector is recomputed per edge to keep the code explicit.
        u_by_chart[chart] = [float(x) for x in u]
        u_norm2_by_chart[chart] = float(norm2)
        if op_packet.get("no_false_pass") is not True:
            raise ValueError(f"expected no_false_pass=true in {ref}")

    # Transition edges from the exported global transition object.
    ops = transition.get("transition_operators", {})
    if not isinstance(ops, dict) or not ops:
        raise ValueError("unexpected selector transition operator map shape")

    edges: list[Edge] = []
    missing_operator_packets: list[str] = []
    for edge_key, spec in ops.items():
        if not isinstance(spec, dict):
            raise ValueError(f"unexpected operator spec for {edge_key}")
        src, dst = _parse_edge_key(edge_key)
        operator_ref = spec.get("operator_ref")
        if not isinstance(operator_ref, str):
            raise ValueError(f"missing operator_ref in transition operator spec {edge_key}")
        edges.append(Edge(src=src, dst=dst, key=edge_key, operator_ref=operator_ref))
        if not (REPO / operator_ref).is_file():
            missing_operator_packets.append(operator_ref)

    # Basic atlas coherence checks (chart list and overlap domain declarations).
    atlas_charts = atlas.get("charts", {})
    atlas_overlap_domains = atlas.get("overlap_domains", {})
    atlas_transition_map = atlas.get("transitions", {})
    atlas_ok = isinstance(atlas_charts, dict) and isinstance(atlas_overlap_domains, dict) and isinstance(atlas_transition_map, dict)

    # Numerical audits
    tol = 1e-12
    max_orth_residual = 0.0
    max_projector_transport_residual = 0.0
    max_projector_cocycle_residual = 0.0

    per_edge: dict[str, dict[str, Any]] = {}
    operator_matrix_cache: dict[str, list[list[float]]] = {}

    def load_operator_matrix(operator_ref: str) -> list[list[float]]:
        if operator_ref in operator_matrix_cache:
            return operator_matrix_cache[operator_ref]
        packet = _read_json(REPO / operator_ref)
        if packet.get("no_false_pass") is not True:
            raise ValueError(f"expected no_false_pass=true in {operator_ref}")
        outputs = packet.get("outputs", {})
        if not isinstance(outputs, dict) or "O12" not in outputs:
            # Most operators store the full matrix under a key matching the edge (e.g. O14/O25/O35...).
            # We accept the first 12x12 matrix-valued list in outputs as the carrier operator.
            candidate = None
            for v in outputs.values():
                if (
                    isinstance(v, list)
                    and len(v) == 12
                    and all(isinstance(row, list) and len(row) == 12 for row in v)
                ):
                    candidate = v
                    break
            if candidate is None:
                raise ValueError(f"could not locate 12x12 operator matrix in {operator_ref}")
            mat = [[float(x) for x in row] for row in candidate]
        else:
            mat = [[float(x) for x in row] for row in outputs["O12"]]
        operator_matrix_cache[operator_ref] = mat
        return mat

    # Transport check on the projector section encoded by the local u_m vectors.
    for e in edges:
        if e.src not in u_by_chart or e.dst not in u_by_chart:
            continue
        o = load_operator_matrix(e.operator_ref)
        ot = _mat_t(o)
        n = len(o)
        if n != 12 or len(o[0]) != 12:
            raise ValueError(f"unexpected operator matrix shape in {e.operator_ref}")

        # Orthogonality: O^T O ≈ I.
        eye = _mat_eye(12)
        orth = _mat_mul(ot, o)
        orth_res = _mat_sub_max_abs(orth, eye)
        max_orth_residual = max(max_orth_residual, orth_res)

        # Projector transport: P_dst ≈ O P_src O^T.
        p_src, _ = _outer_projector(u_by_chart[e.src])
        p_dst, _ = _outer_projector(u_by_chart[e.dst])
        transported = _mat_mul(_mat_mul(o, p_src), ot)
        proj_res = _mat_sub_max_abs(transported, p_dst)
        max_projector_transport_residual = max(max_projector_transport_residual, proj_res)

        overlap_key = f"{e.src}_INTERSECT_{e.dst}"
        per_edge[e.key] = {
            "src": e.src,
            "dst": e.dst,
            "operator_ref": e.operator_ref,
            "orthogonality_max_abs_residual": orth_res,
            "projector_transport_max_abs_residual": proj_res,
            "atlas_overlap_domain_declared": (atlas_ok and overlap_key in atlas_overlap_domains),
            "atlas_transition_declared": (atlas_ok and e.key in atlas_transition_map),
        }

    # Projector-level cocycle checks for all triples where all three directed edges exist.
    edge_map: dict[tuple[str, str], Edge] = {(e.src, e.dst): e for e in edges}
    for i in chart_ids:
        for j in chart_ids:
            if i == j:
                continue
            for k in chart_ids:
                if k in (i, j):
                    continue
                e_ij = edge_map.get((i, j))
                e_jk = edge_map.get((j, k))
                e_ik = edge_map.get((i, k))
                if not (e_ij and e_jk and e_ik):
                    continue
                if i not in u_by_chart:
                    continue
                oij = load_operator_matrix(e_ij.operator_ref)
                ojk = load_operator_matrix(e_jk.operator_ref)
                oik = load_operator_matrix(e_ik.operator_ref)

                p_i, _ = _outer_projector(u_by_chart[i])

                comp = _mat_mul(ojk, oij)
                comp_t = _mat_t(comp)
                direct_t = _mat_t(oik)

                via = _mat_mul(_mat_mul(comp, p_i), comp_t)
                direct = _mat_mul(_mat_mul(oik, p_i), direct_t)
                cocycle_res = _mat_sub_max_abs(via, direct)
                if cocycle_res > max_projector_cocycle_residual:
                    max_projector_cocycle_residual = cocycle_res

    passes = (
        not missing_local_operator_packets
        and not missing_operator_packets
        and max_orth_residual <= tol
        and max_projector_transport_residual <= tol
        and max_projector_cocycle_residual <= tol
    )
    status = "PASS_GLUING_CONSISTENT" if passes else "FAIL_GLUING_INCONSISTENT"

    artifact = {
        "stage": "P474",
        "date": "2026-03-16",
        "status": status,
        "goal": "audit_that_exported_global_projective_selector_state_object_is_projector_level_glued_consistently_by_exported_global_selector_transition_operators",
        "inputs": {
            "state": _p(STATE),
            "atlas": _p(ATLAS),
            "transition": _p(TRANSITION),
        },
        "checked": {
            "charts": chart_ids,
            "local_operator_refs": local_operator_refs,
            "missing_local_operator_packets": missing_local_operator_packets,
            "missing_operator_packets": missing_operator_packets,
            "u_norm2_by_chart": u_norm2_by_chart,
            "tolerance": tol,
        },
        "results": {
            "max_orthogonality_max_abs_residual": max_orth_residual,
            "max_projector_transport_max_abs_residual": max_projector_transport_residual,
            "max_projector_cocycle_max_abs_residual": max_projector_cocycle_residual,
            "per_edge": per_edge,
        },
        "no_false_pass": True,
        "hard_limits": [
            "no theorem-level pass",
            "no full-closure pass",
            "no sign-sensitive/directed selector state claim (projective/ray-level only)",
            "no claim of global QW-2191 discharge",
        ],
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "overall_pass": bool(passes),
        "max_projector_transport_max_abs_residual": max_projector_transport_residual,
        "max_projector_cocycle_max_abs_residual": max_projector_cocycle_residual,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

