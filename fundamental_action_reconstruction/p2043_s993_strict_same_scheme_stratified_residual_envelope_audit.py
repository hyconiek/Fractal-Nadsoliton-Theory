#!/usr/bin/env python3
"""P2043 S993: stratified residual envelope audit.

Extends P2042 by auditing envelope statistics stratified by norm bucket:
- per-bucket sup residual with uncertainty buffer,
- monotonicity check versus norm bucket,
- bootstrap CI for global sup_residual_with_buffer_linf.

This remains a witness/audit stage and does not discharge C3.
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
OUT = GEN / "p2043_s993_strict_same_scheme_stratified_residual_envelope_audit.json"
MD = GEN / "p2043_s993_strict_same_scheme_stratified_residual_envelope_audit.md"

SCHEMA_VERSION = "p2043_s993_v1"
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


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2042 = load("p2042_s992_strict_same_scheme_operator_level_residual_envelope_sweep.json")
    checks_2042 = p2042.get("gatekeeper_checks") or {}

    ready = (
        p2042.get("result_kind")
        == "PASS_OPERATOR_LEVEL_RESIDUAL_ENVELOPE_SWEEP_WITH_UNCERTAINTY_BUFFER__C3_STILL_OPEN"
        and as_bool(checks_2042.get("worst_case_finite"))
        and as_bool(checks_2042.get("sample_count_ge_40"))
    )

    samples = p2042.get("samples") or []
    if not samples:
        stratified = {}
        norms_sorted: list[float] = []
        sup_seq: list[float] = []
    else:
        stratified: dict[float, list[float]] = {}
        for row in samples:
            n = float(row["norm_bucket"])
            r = float(row["residual_with_uncertainty_buffer_linf"])
            stratified.setdefault(n, []).append(r)

        norms_sorted = sorted(stratified.keys())
        sup_seq = [max(stratified[n]) for n in norms_sorted]

    monotonic_non_decreasing = all(sup_seq[i] <= sup_seq[i + 1] + 1e-18 for i in range(len(sup_seq) - 1))

    per_bucket = []
    for n, s in zip(norms_sorted, sup_seq):
        per_bucket.append(
            {
                "norm_bucket": n,
                "sample_count": len(stratified[n]),
                "sup_residual_with_buffer_linf": s,
                "mean_residual_with_buffer_linf": float(np.mean(stratified[n])),
            }
        )

    global_sup = float(max(sup_seq)) if sup_seq else math.inf

    # Bootstrap CI for sup statistic over full sample list.
    full = np.array([float(r["residual_with_uncertainty_buffer_linf"]) for r in samples], dtype=float)
    rng = np.random.default_rng(2043)
    B = 600
    boot_sup = []
    if len(full) > 0:
        for _ in range(B):
            idx = rng.integers(0, len(full), size=len(full))
            boot_sup.append(float(np.max(full[idx])))
    boot_sup = np.array(boot_sup, dtype=float) if boot_sup else np.array([], dtype=float)

    if len(boot_sup) > 0:
        ci_low = float(np.quantile(boot_sup, 0.025))
        ci_high = float(np.quantile(boot_sup, 0.975))
    else:
        ci_low = math.inf
        ci_high = math.inf

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2043",
        "stage_id": "S993",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRATIFIED_RESIDUAL_ENVELOPE_AUDIT_WITH_BOOTSTRAP_CI__C3_STILL_OPEN"
            if ready and math.isfinite(global_sup)
            else "OPEN_STRATIFIED_RESIDUAL_ENVELOPE_AUDIT_BLOCKED"
        ),
        "route": "strict_only_stratified_envelope_audit",
        "depends_on": {
            "p2042_present": p2042.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2042_json_sha256": file_sha256(GEN / "p2042_s992_strict_same_scheme_operator_level_residual_envelope_sweep.json"),
        },
        "audit_scope": {
            "controlled_background_pair": (p2042.get("sweep_spec") or {}).get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "strata_key": "norm_bucket",
        },
        "stratified_envelope": {
            "per_bucket": per_bucket,
            "norm_bucket_order": norms_sorted,
            "sup_sequence": sup_seq,
            "monotonic_non_decreasing_vs_norm": monotonic_non_decreasing,
            "global_sup_residual_with_buffer_linf": global_sup,
        },
        "bootstrap_ci": {
            "bootstrap_seed": 2043,
            "replicates": B,
            "statistic": "sup_residual_with_buffer_linf",
            "ci_level": 0.95,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "ci_width": (ci_high - ci_low) if math.isfinite(ci_high) and math.isfinite(ci_low) else math.inf,
        },
        "c3_gate_update": {
            "C3_stratified_envelope_audit": "COMPUTED",
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
            "per_bucket_nonempty": all(row["sample_count"] > 0 for row in per_bucket),
            "bootstrap_ci_finite": math.isfinite(ci_low) and math.isfinite(ci_high),
            "global_sup_finite": math.isfinite(global_sup),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2043 S993: stratified residual envelope audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Exported per-norm-bucket sup envelope, monotonicity check vs norm, and",
        "bootstrap CI for global `sup_residual_with_buffer_linf`.",
        "",
        "## Gate update",
        "",
        "- `C3`: stratified audit computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
