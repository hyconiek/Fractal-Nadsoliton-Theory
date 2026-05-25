#!/usr/bin/env python3
"""P2045 S995: seed x norm joint robustness audit.

Computes worst-case residual-with-buffer table per (seed, norm_bucket) and
exports ranking-stability metrics (Spearman/Kendall style) across seeds.
Audit-level only; C3 remains open.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2045_s995_strict_same_scheme_seed_norm_joint_robustness_audit.json"
MD = GEN / "p2045_s995_strict_same_scheme_seed_norm_joint_robustness_audit.md"

SCHEMA_VERSION = "p2045_s995_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BASIS = ["R2_bar", "Ric2_bar", "Riem2_bar"]


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(v: Any) -> bool:
    return bool(v is True)


def rank_positions_desc(values: list[float]) -> list[int]:
    idx = sorted(range(len(values)), key=lambda i: (-values[i], i))
    pos = [0] * len(values)
    for r, i in enumerate(idx):
        pos[i] = r + 1
    return pos


def spearman_from_positions(p: list[int], q: list[int]) -> float:
    n = len(p)
    if n < 2:
        return 1.0
    d2 = sum((p[i] - q[i]) ** 2 for i in range(n))
    return 1.0 - 6.0 * d2 / (n * (n * n - 1))


def kendall_tau_bruteforce(p: list[int], q: list[int]) -> float:
    n = len(p)
    if n < 2:
        return 1.0
    conc = 0
    disc = 0
    for i in range(n):
        for j in range(i + 1, n):
            a = p[i] - p[j]
            b = q[i] - q[j]
            s = a * b
            if s > 0:
                conc += 1
            elif s < 0:
                disc += 1
    denom = n * (n - 1) / 2
    return (conc - disc) / denom if denom > 0 else 1.0


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2044 = load("p2044_s994_strict_same_scheme_seed_sensitivity_meta_envelope_audit.json")
    checks_2044 = p2044.get("gatekeeper_checks") or {}

    ready = (
        p2044.get("result_kind")
        == "PASS_SEED_SENSITIVITY_META_ENVELOPE_AUDIT_WITH_INTERSEED_CI__C3_STILL_OPEN"
        and as_bool(checks_2044.get("seed_count_ge_5"))
        and as_bool(checks_2044.get("meta_sup_finite"))
    )

    spec = p2044.get("seed_sensitivity_spec") or {}
    seeds = spec.get("seeds") or []
    norm_buckets = spec.get("norm_buckets") or [0.01, 0.015, 0.02, 0.025]
    key_order = [str(x) for x in norm_buckets]

    per_seed = p2044.get("per_seed_results") or []

    table = []
    ranks = {}
    for row in per_seed:
        s = int(row["seed"])
        bucket = row.get("bucket_sup_residual_with_buffer_linf") or {}
        vals = [float(bucket[k]) for k in key_order]
        ranks[s] = rank_positions_desc(vals)
        for k, v in zip(key_order, vals):
            table.append({"seed": s, "norm_bucket": float(k), "worst_case_residual_with_buffer_linf": v})

    # Pairwise ranking stability across seeds.
    pairwise = []
    sp_vals = []
    kd_vals = []
    for i in range(len(seeds)):
        for j in range(i + 1, len(seeds)):
            si = int(seeds[i])
            sj = int(seeds[j])
            if si not in ranks or sj not in ranks:
                continue
            spc = spearman_from_positions(ranks[si], ranks[sj])
            kdt = kendall_tau_bruteforce(ranks[si], ranks[sj])
            sp_vals.append(spc)
            kd_vals.append(kdt)
            pairwise.append({"seed_i": si, "seed_j": sj, "spearman": spc, "kendall_tau": kdt})

    mean_sp = float(np.mean(sp_vals)) if sp_vals else math.inf
    mean_kd = float(np.mean(kd_vals)) if kd_vals else math.inf
    min_sp = float(np.min(sp_vals)) if sp_vals else math.inf
    min_kd = float(np.min(kd_vals)) if kd_vals else math.inf

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2045",
        "stage_id": "S995",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_SEED_NORM_JOINT_ROBUSTNESS_AUDIT_WITH_RANK_STABILITY__C3_STILL_OPEN"
            if ready and math.isfinite(mean_sp) and math.isfinite(mean_kd)
            else "OPEN_SEED_NORM_JOINT_ROBUSTNESS_AUDIT_BLOCKED"
        ),
        "route": "strict_only_seed_norm_joint_robustness_audit",
        "depends_on": {
            "p2044_present": p2044.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2044_json_sha256": file_sha256(GEN / "p2044_s994_strict_same_scheme_seed_sensitivity_meta_envelope_audit.json"),
        },
        "audit_scope": {
            "controlled_background_pair": spec.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "seeds": seeds,
            "norm_buckets": norm_buckets,
        },
        "seed_norm_worst_case_table": table,
        "ranking_stability": {
            "pairwise": pairwise,
            "mean_spearman": mean_sp,
            "mean_kendall_tau": mean_kd,
            "min_spearman": min_sp,
            "min_kendall_tau": min_kd,
        },
        "c3_gate_update": {
            "C3_seed_norm_joint_robustness_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "theorem-grade operator identity across full background family",
                "cross-background finite-part transport identity theorem",
                "global finite-part lock/cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "table_nonempty": len(table) > 0,
            "pairwise_nonempty": len(pairwise) > 0,
            "mean_spearman_finite": math.isfinite(mean_sp),
            "mean_kendall_finite": math.isfinite(mean_kd),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2045 S995: seed x norm joint robustness audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Exported worst-case residual table per (seed, norm_bucket) and ranking",
        "stability metrics across seeds (Spearman + Kendall tau).",
        "",
        "## Gate update",
        "",
        "- `C3`: joint robustness audit computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
