#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2317_s1267_strict_fourier_degeneracy_existing_lane_audit_probe.json"
MD = GEN / "p2317_s1267_strict_fourier_degeneracy_existing_lane_audit_probe.md"

N = 12
OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8
TOL = 1e-10

SOURCE_FILES = {
    "P2315_KERNEL_SPECTRUM": GEN / "p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.json",
    "F454_SHANNON_MODE_ASSIGNMENT_NOTE": ROOT / "F454_CURRENT_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_MODE_INDEX_ASSIGNMENT_PACKET.md",
    "F459_DIAGONAL_HESSIAN_NOTE": ROOT / "F459_CURRENT_STRICT_DIAGONAL_LOCAL_PSI_HESSIAN_EIGENSYSTEM_VALUE_INSTANTIATION_PACKET.md",
    "N494_UNIQUENESS_UP_TO_CONJUGATION": ROOT / "N494_CURRENT_FIRST_STRICT_QW2190_DIAGONAL_LOCAL_MODE_INDEX_CANONICALIZATION_UNIQUENESS_UP_TO_CONJUGATION_THEOREM.md",
    "P454_O2_GAUGE_EQUIVALENCE_NOTE": ROOT / "P454_CURRENT_STRICT_QW2191_O2_ROTATION_GAUGE_EQUIVALENCE_AUDIT_PROBE.md",
    "P449_DIAGONAL_LOCAL_DEFECTS": GEN / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe.json",
    "SHANNON_MODE_ASSIGNMENT": GEN / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json",
    "DIAGONAL_MODE_ASSIGNMENT": GEN / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json",
    "PSI_HESSIAN_EIGENSYSTEM": GEN / "psi_hessian_diagonal_local_strict_derived_value_instantiated_v1.json",
    "R_ORD_Z12_REFERENCE": GEN / "r_ord_z12_v1_reference_distribution.json",
}

GREP_PATTERNS = (
    "Fourier pair",
    "pair plane",
    "pair_m",
    "O(2)",
    "QW-2191",
    "mode-index assignment",
    "element-order reference",
    "diagonal/local",
    "degenerate",
    "selector closure",
)


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def strict_kernel(distance: int) -> float:
    if distance == 0:
        return 0.0
    return math.cos(OMEGA * distance + PHI) / (1.0 + BETA * (distance ** ETA))


def cyclic_distance(i: int, j: int) -> int:
    raw = abs(i - j) % N
    return min(raw, N - raw)


def kernel_matrix() -> list[list[float]]:
    return [[strict_kernel(cyclic_distance(i, j)) if i != j else 0.0 for j in range(N)] for i in range(N)]


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def mat_vec(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [dot(row, vector) for row in matrix]


def real_fourier_pair(m: int) -> tuple[list[float], list[float]]:
    c = [math.sqrt(2.0 / N) * math.cos(2.0 * math.pi * m * x / N) for x in range(N)]
    s = [-math.sqrt(2.0 / N) * math.sin(2.0 * math.pi * m * x / N) for x in range(N)]
    return c, s


def pair_blocks_for_kernel(matrix: list[list[float]]) -> list[dict[str, Any]]:
    blocks: list[dict[str, Any]] = []
    for m in range(1, N // 2):
        c, s = real_fourier_pair(m)
        kc = mat_vec(matrix, c)
        ks = mat_vec(matrix, s)
        block = [[dot(c, kc), dot(c, ks)], [dot(s, kc), dot(s, ks)]]
        trace_half = 0.5 * (block[0][0] + block[1][1])
        scalar_residual = max(
            abs(block[0][0] - trace_half),
            abs(block[1][1] - trace_half),
            abs(block[0][1]),
            abs(block[1][0]),
        )
        split = abs(block[0][0] - block[1][1]) + abs(block[0][1]) + abs(block[1][0])
        blocks.append({
            "pair": f"pair{m}",
            "m": m,
            "partner_modes": [m, N - m],
            "real_fourier_block": block,
            "lambda_scalar": trace_half,
            "scalar_identity_residual": scalar_residual,
            "orientation_splitting_measure": split,
            "kernel_alone_degenerate": scalar_residual < TOL,
        })
    return blocks


def corpus_hits() -> list[dict[str, Any]]:
    paths = sorted(
        set(SOURCE_FILES.values())
        | set(ROOT.glob("F4*_CURRENT_STRICT*MODE_INDEX*.md"))
        | set(ROOT.glob("N49*_CURRENT*QW2191*.md"))
        | set(ROOT.glob("P45*_CURRENT_STRICT*QW2191*.md"))
        | set(GEN.glob("*mode_index_assignment*json"))
        | set(GEN.glob("p4*_current_strict*o2*json"))
    )
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not path.exists() or path.is_dir():
            continue
        text = read_text(path)
        lowered = text.lower()
        count = sum(lowered.count(pattern.lower()) for pattern in GREP_PATTERNS)
        if count == 0:
            continue
        first_line = 0
        first_excerpt = ""
        for idx, line in enumerate(text.splitlines(), start=1):
            if any(pattern.lower() in line.lower() for pattern in GREP_PATTERNS):
                first_line = idx
                first_excerpt = line.strip()[:240]
                break
        rows.append({
            "path": path.relative_to(REPO).as_posix(),
            "pattern_hit_count": count,
            "first_hit_line": first_line,
            "first_hit_excerpt": first_excerpt,
        })
    rows.sort(key=lambda row: (-int(row["pattern_hit_count"]), row["path"]))
    return rows


def diagonal_local_lane(p449: dict[str, Any], diag_assignment: dict[str, Any]) -> dict[str, Any]:
    pairs = p449.get("computed", {}).get("pairs", {})
    rows = []
    for name in sorted(pairs, key=lambda key: int(key.replace("pair", ""))):
        row = pairs[name]
        eig = row.get("eigenvalues", {})
        split = abs(float(eig.get("lambda_plus", 0.0)) - float(eig.get("lambda_minus", 0.0)))
        rows.append({
            "pair": name,
            "defect_abs": float(row.get("F", {}).get("abs", 0.0)),
            "cuts_O2_by_absF_nonzero": bool(row.get("cuts_O2_by_absF_nonzero")),
            "theta_star": float(row.get("theta_star", 0.0)),
            "eigenvalue_split": split,
        })
    return {
        "source": "P449 + mode_index_assignment_canonical_local_diagonal_strict_derived_v1",
        "scope": diag_assignment.get("scope"),
        "status": diag_assignment.get("status"),
        "all_pairs_cut": bool(p449.get("computed", {}).get("all_pairs_cut")),
        "pair_rows": rows,
        "min_defect_abs": min((row["defect_abs"] for row in rows), default=0.0),
        "min_eigenvalue_split": min((row["eigenvalue_split"] for row in rows), default=0.0),
        "lane_limitation": "Lane-scoped axis/Z2 mode-index export only; not a kernel-alone selector proof and not a Task-3 provider-to-policy-margin bridge.",
    }


def shannon_lane(shannon_assignment: dict[str, Any], r_ord: dict[str, Any]) -> dict[str, Any]:
    pairs = shannon_assignment.get("outputs", {}).get("pairs", {})
    rows = []
    for name in sorted(pairs, key=lambda key: int(key.replace("pair", ""))):
        row = pairs[name]
        eig = row.get("eigenvalues_on_diag_ord_restriction", {})
        split = abs(float(eig.get("lambda_plus", 0.0)) - float(eig.get("lambda_minus", 0.0)))
        rows.append({
            "pair": name,
            "defect_abs": float(row.get("F_2m_ord", {}).get("abs", 0.0)),
            "theta_star": float(row.get("theta_star", 0.0)),
            "objective_minimizer_vector": row.get("objective_minimizer_vector"),
            "eigenvalue_split": split,
        })
    return {
        "source": "F454/mode_index_assignment_shannon_element_order_reference_strict_core_v1 + r_ord_z12_v1",
        "scope": shannon_assignment.get("scope"),
        "status": shannon_assignment.get("status"),
        "r_ord_translation_invariance_note": r_ord.get("invariance_notes", []),
        "pair_rows": rows,
        "all_defects_nonzero": all(row["defect_abs"] > TOL for row in rows),
        "min_defect_abs": min((row["defect_abs"] for row in rows), default=0.0),
        "min_eigenvalue_split": min((row["eigenvalue_split"] for row in rows), default=0.0),
        "lane_limitation": "Strict Shannon element-order lane exports axis/Z2 mode-index data but explicitly does not claim strict-core selector closure, global QW-2191 discharge, or ToE closure.",
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {sid: load_json(path) for sid, path in SOURCE_FILES.items() if path.suffix == ".json"}
    source_hashes = {sid: sha256_file(path) for sid, path in SOURCE_FILES.items() if path.exists()}

    k_matrix = kernel_matrix()
    kernel_blocks = pair_blocks_for_kernel(k_matrix)
    p2315 = artifacts["P2315_KERNEL_SPECTRUM"]
    p449 = artifacts["P449_DIAGONAL_LOCAL_DEFECTS"]
    shannon_assignment = artifacts["SHANNON_MODE_ASSIGNMENT"]
    diag_assignment = artifacts["DIAGONAL_MODE_ASSIGNMENT"]
    r_ord = artifacts["R_ORD_Z12_REFERENCE"]

    existing_evidence = {
        "grep_patterns": list(GREP_PATTERNS),
        "hit_count": len(corpus_hits()),
        "top_hits": corpus_hits()[:30],
        "already_present_packets": [
            "P2315: kernel-alone Z12 strict spectrum and pair-plane degeneracy computation",
            "QW-2191/N494/P454: O(2) pair-plane obstruction and conjugation-equivalence hygiene",
            "P449/F453/F459: diagonal/local strict-derived multi-pair O(2)->Z2 lane",
            "F454/N514: Shannon element-order reference multi-pair O(2)->Z2 lane",
        ],
    }

    kernel_alone_audit = {
        "matrix_scope": "strict_gate_kernel_on_Z12_only",
        "parameters": {"n": N, "omega": OMEGA, "phi": PHI, "beta": BETA, "eta": ETA},
        "real_fourier_pair_blocks": kernel_blocks,
        "all_pair_blocks_scalar_identity": all(row["kernel_alone_degenerate"] for row in kernel_blocks),
        "max_scalar_identity_residual": max((row["scalar_identity_residual"] for row in kernel_blocks), default=0.0),
        "max_orientation_splitting_measure": max((row["orientation_splitting_measure"] for row in kernel_blocks), default=0.0),
        "p2315_loaded_status": p2315.get("status"),
        "p2315_pair_degeneracy_verified": p2315.get(
            "strict_schematic_lagrangian_eom_kernel_spectrum_probe", {}
        ).get("kernel_spectrum_audit", {}).get("pair_degeneracy_report", {}).get("all_pair_planes_degenerate"),
        "verdict": "Kernel-alone Fourier pair blocks are scalar multiples of the identity, so no orientation inside any pair plane is selected.",
    }

    diagonal = diagonal_local_lane(p449, diag_assignment)
    shannon = shannon_lane(shannon_assignment, r_ord)

    lane_comparison = {
        "kernel_alone": {
            "orientation_selected": False,
            "reason": "circulant/translation-invariant kernel gives scalar identity blocks on each real Fourier pair plane",
        },
        "diagonal_local_lane": {
            "orientation_selected_up_to_Z2_in_lane": bool(diagonal["all_pairs_cut"]),
            "min_defect_abs": diagonal["min_defect_abs"],
            "min_eigenvalue_split": diagonal["min_eigenvalue_split"],
            "not_a_task3_bridge": True,
        },
        "shannon_element_order_lane": {
            "orientation_selected_up_to_Z2_in_lane": bool(shannon["all_defects_nonzero"]),
            "min_defect_abs": shannon["min_defect_abs"],
            "min_eigenvalue_split": shannon["min_eigenvalue_split"],
            "not_a_task3_bridge": True,
        },
        "route_update": "Do not redo kernel spectrum as if new.  Existing FAR already has lane-scoped O(2)->Z2 cuts; the missing proof object is a bridge from one admissible lane source to Task-3 provider_lift_per_step / policy-margin semantics.",
    }

    theorem_export = {
        "theorem_name": "P2317 Fourier degeneracy and existing-lane audit",
        "claim": "The strict kernel alone cannot choose orientation in Fourier pair planes, and this was already computed in P2315/QW-2191-family evidence.  The repo also contains diagonal/local and Shannon element-order lane-scoped O(2)->Z2 computations, but these are not yet a Task-3 provider-to-policy-margin bridge.",
        "proof_parts": [
            "Direct recomputation: each real Fourier pair block of the strict Z12 kernel is a scalar identity block within numerical tolerance.",
            "P2315 confirms all five (k,12-k) pair planes are degenerate for the kernel-alone spectrum and after scalar 4 ln 2 scaling.",
            "P449/F453/F459 provide a diagonal/local strict-derived lane with nonzero multi-pair defects and residual Z2 only in that lane.",
            "F454/N514 provide a Shannon element-order lane with nonzero multi-pair defects and residual Z2 only in that lane.",
            "All lane exports preserve hard limits: no global QW-2191 discharge, no strict-core selector closure, no ToE closure, and no Task-3 G1/G3 update.",
        ],
        "next_required_object": "A theorem-grade bridge from diagonal/local or Shannon lane source data to Task-3 provider_lift_per_step / policy-margin lift, or a nonexistence theorem for that bridge class.",
    }

    probe = {
        "probe_id": "P2317_S1267_STRICT_FOURIER_DEGENERACY_EXISTING_LANE_AUDIT",
        "source_hashes": source_hashes,
        "existing_evidence_grep_audit": existing_evidence,
        "kernel_alone_pair_block_computation": kernel_alone_audit,
        "diagonal_local_lane_o2_to_z2_computation": diagonal,
        "shannon_element_order_lane_o2_to_z2_computation": shannon,
        "lane_comparison_and_route_decision": lane_comparison,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "far_grep_found_existing_fourier_selector_evidence": existing_evidence["hit_count"] >= 10,
        "p2315_loaded": p2315.get("packet_id") == "P2315",
        "p2315_pair_degeneracy_verified": kernel_alone_audit["p2315_pair_degeneracy_verified"] is True,
        "direct_kernel_pair_blocks_scalar_identity": kernel_alone_audit["all_pair_blocks_scalar_identity"],
        "direct_kernel_pair_block_residual_small": kernel_alone_audit["max_scalar_identity_residual"] < TOL,
        "diagonal_local_lane_all_pairs_cut": diagonal["all_pairs_cut"],
        "diagonal_local_lane_nonzero_min_defect": diagonal["min_defect_abs"] > TOL,
        "shannon_lane_all_defects_nonzero": shannon["all_defects_nonzero"],
        "shannon_lane_nonzero_min_defect": shannon["min_defect_abs"] > TOL,
        "existing_lane_exports_not_promoted_to_task3_bridge": True,
        "strict_g1_g3_not_updated": True,
        "no_selector_closure_claimed": True,
        "no_qw2191_global_discharge_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2317_s1267_v1",
        "packet_id": "P2317",
        "stage_id": "S1267",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_KERNEL_DEGENERACY_ALREADY_COMPUTED_EXISTING_LANES_NOT_TASK3_BRIDGE",
        "result_kind": "STRICT_FOURIER_PAIR_DEGENERACY_RECOMPUTATION_AND_EXISTING_SELECTOR_LANE_AUDIT_NO_G1_G3_UPDATE",
        "strict_fourier_degeneracy_existing_lane_audit_probe": probe,
        "recommended_next_honest_step": theorem_export["next_required_object"],
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")

    md_lines = [
        "# P2317 S1267: Fourier degeneracy and existing selector-lane audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## FAR grep result",
        f"- Existing Fourier/O(2)/selector evidence hits: `{existing_evidence['hit_count']}`.",
        "- Already present: `P2315`, `QW-2191/N494/P454`, `P449/F453/F459`, and `F454/N514` lanes.",
        "",
        "## Direct recomputation",
        f"- All strict-kernel real Fourier pair blocks scalar identity: `{kernel_alone_audit['all_pair_blocks_scalar_identity']}`.",
        f"- Max scalar-identity residual: `{kernel_alone_audit['max_scalar_identity_residual']:.3e}`.",
        "- Verdict: kernel alone does not choose orientation inside any `(k,12-k)` pair plane.",
        "",
        "## Existing lane computations",
        f"- Diagonal/local lane all pairs cut: `{diagonal['all_pairs_cut']}`; min defect abs: `{diagonal['min_defect_abs']:.6g}`.",
        f"- Shannon element-order lane all defects nonzero: `{shannon['all_defects_nonzero']}`; min defect abs: `{shannon['min_defect_abs']:.6g}`.",
        "- These are lane-scoped O(2)->Z2 axis exports, not a Task-3 provider-to-policy-margin bridge.",
        "",
        "## Route decision",
        f"{lane_comparison['route_update']}",
        "",
        "## Theorem fingerprint",
        f"`{probe['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
