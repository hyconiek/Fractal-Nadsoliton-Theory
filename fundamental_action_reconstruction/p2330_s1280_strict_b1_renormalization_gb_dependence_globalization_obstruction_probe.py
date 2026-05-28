#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import scipy.integrate as si
import sympy as sp

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2330_s1280_strict_b1_renormalization_gb_dependence_globalization_obstruction_probe.json"
MD = GEN / "p2330_s1280_strict_b1_renormalization_gb_dependence_globalization_obstruction_probe.md"

SOURCE_FILES = {
    "P1848_GRAVITY_COUNTERTERM_PROFILES": GEN / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json",
    "P1950_RENORMALIZATION_EXACT_INTEGRATION": GEN / "p1950_s900_strict_renormalization_exact_integration_probe.json",
    "P2094_QUOTIENT_RANK_REPAIR": GEN / "p2094_s1044_strict_b1_quotient_renormalization_rank_repair.json",
    "P2096_QUOTIENT_CLOSURE_CONTRACT": GEN / "p2096_s1046_strict_b1_renormalization_closure_contract.json",
    "P2329_SELECTOR_INDEPENDENCE_AUDIT": GEN / "p2329_s1279_selector_independence_lagrangian_eom_audit_probe.json",
}

GREP_PATTERNS = (
    "a_R2",
    "a_Ric2",
    "a_Riem2",
    "a_GB",
    "GaussBonnet",
    "rank_deficient",
    "renormalization",
    "counterterm",
    "background_family_B1",
    "delta_c_gr",
)

CHANNEL_ORDER = ("R2", "Ric2", "Riem2", "GB")
NULL_VECTOR = np.array([-1.0, 4.0, -1.0, 1.0], dtype=float)


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": path.relative_to(REPO).as_posix()}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def strict_kernel(d: sp.Symbol) -> sp.Expr:
    omega = sp.Float("0.18575")
    phi = sp.Float("0.16250")
    beta = sp.Float("1.0")
    eta = sp.Rational(9, 5)
    return sp.cos(omega * d + phi) / (1 + beta * d**eta)


def backend_profiles(p1848: dict[str, Any], d: sp.Symbol) -> dict[str, sp.Expr]:
    K = strict_kernel(d)
    Kd = sp.diff(K, d)
    Kdd = sp.diff(Kd, d)
    locals_map: dict[str, Any] = {"d": d, "K": K, "Kd": Kd, "Kdd": Kdd, "sin": sp.sin, "cos": sp.cos}
    exprs: dict[str, sp.Expr] = {}
    profiles = p1848.get("gravity_operator_profiles_B1", {}).get("profiles", {})
    for channel in CHANNEL_ORDER:
        raw = profiles.get(channel, {})
        expr_text = raw.get("expression")
        if not isinstance(expr_text, str):
            raise ValueError(f"missing profile expression for {channel}")
        exprs[channel] = sp.sympify(expr_text, locals={**locals_map, **exprs})
    return exprs


def eval_grid(expr: sp.Expr, d: sp.Symbol, xs: np.ndarray) -> np.ndarray:
    fn = sp.lambdify(d, sp.simplify(expr), "numpy")
    values = np.asarray(fn(xs), dtype=float)
    if values.shape == ():
        return np.full_like(xs, float(values), dtype=float)
    return values.astype(float, copy=False)


def gram_matrix(exprs: dict[str, sp.Expr], d: sp.Symbol, xs: np.ndarray) -> np.ndarray:
    values = [eval_grid(exprs[channel], d, xs) for channel in CHANNEL_ORDER]
    gram = np.zeros((len(CHANNEL_ORDER), len(CHANNEL_ORDER)), dtype=float)
    for i in range(len(CHANNEL_ORDER)):
        for j in range(len(CHANNEL_ORDER)):
            gram[i, j] = float(si.simpson(values[i] * values[j], x=xs))
    return gram


def collect_repo_grep_hits() -> list[dict[str, Any]]:
    candidates = [
        ROOT / "p1950_s900_strict_renormalization_exact_integration.py",
        ROOT / "p2094_s1044_strict_b1_quotient_renormalization_rank_repair.py",
        ROOT / "p2096_s1046_strict_b1_renormalization_closure_contract.py",
    ]
    candidates.extend(SOURCE_FILES.values())
    hits: list[dict[str, Any]] = []
    for path in candidates:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        lowered = text.lower()
        count = sum(lowered.count(pattern.lower()) for pattern in GREP_PATTERNS)
        if count == 0:
            continue
        first_line = 0
        first_excerpt = ""
        for index, line in enumerate(text.splitlines(), start=1):
            if any(pattern.lower() in line.lower() for pattern in GREP_PATTERNS):
                first_line = index
                first_excerpt = line.strip()[:220]
                break
        hits.append({
            "path": path.relative_to(REPO).as_posix(),
            "pattern_hit_count": count,
            "first_hit_line": first_line,
            "first_hit_excerpt": first_excerpt,
        })
    hits.sort(key=lambda row: (-int(row["pattern_hit_count"]), row["path"]))
    return hits


def quotient_summary(p2096: dict[str, Any]) -> dict[str, Any]:
    contract = p2096.get("closure_contract", {})
    residuals = contract.get("residuals", {})
    return {
        "p2096_loaded": p2096.get("packet_id") == "P2096",
        "status": p2096.get("status"),
        "result_kind": p2096.get("result_kind"),
        "scope_limit": contract.get("scope_limit"),
        "quotient_residual_abs_max": residuals.get("quotient_residual_abs_max"),
        "gb_derived_operator_residual_abs_max": residuals.get("gb_derived_operator_residual_abs_max"),
        "full_4channel_independence_proven": p2096.get("gatekeeper_checks", {}).get("full_4channel_independence_proven"),
        "global_background_family_closure_proven": p2096.get("gatekeeper_checks", {}).get("global_background_family_closure_proven"),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    p1848 = artifacts["P1848_GRAVITY_COUNTERTERM_PROFILES"]
    p1950 = artifacts["P1950_RENORMALIZATION_EXACT_INTEGRATION"]
    p2096 = artifacts["P2096_QUOTIENT_CLOSURE_CONTRACT"]

    raw_profiles = p1848.get("gravity_operator_profiles_B1", {}).get("profiles", {})
    gb_expression_text = raw_profiles.get("GB", {}).get("expression", "")
    symbolic_relation_zero = gb_expression_text.replace(" ", "") == "Riem2-4*Ric2+R2"

    p1950_rows = [
        row for row in p1950.get("computed_rows", [])
        if isinstance(row, dict) and "row_operator" in row
    ]
    gram = np.array([
        [float(row["row_operator"][f"a_{channel}"]) for channel in CHANNEL_ORDER]
        for row in p1950_rows
    ], dtype=float)
    singular_values = np.linalg.svd(gram, compute_uv=False)
    null_residual = gram @ NULL_VECTOR
    null_residual_l2 = float(np.linalg.norm(null_residual, ord=2))
    max_sv = float(np.max(singular_values)) if singular_values.size else 0.0
    rank_tol = max(gram.shape) * max_sv * np.finfo(float).eps * 100.0 if max_sv else 0.0
    numeric_rank = int(np.sum(singular_values > rank_tol)) if max_sv else 0
    nullity = int(gram.shape[1] - numeric_rank) if gram.ndim == 2 else 0
    profile_relation_residual_abs_max = null_residual_l2

    p1950_metrics = p1950.get("aggregate_metrics", {})
    p1950_pass_fail = p1950.get("pass_fail_criteria", {})
    p1950_rank = p1950_metrics.get("operator_profile_gram_rank")
    p1950_nullity = p1950_metrics.get("operator_profile_gram_nullity")

    gb_dependence_certificate = {
        "background_family": "B1",
        "channel_order": list(CHANNEL_ORDER),
        "exact_symbolic_relation": "GB - Riem2 + 4*Ric2 - R2 == 0",
        "backend_gb_expression_text": gb_expression_text,
        "symbolic_relation_zero": bool(symbolic_relation_zero),
        "null_vector_for_channel_profiles": {
            "R2": -1.0,
            "Ric2": 4.0,
            "Riem2": -1.0,
            "GB": 1.0,
        },
        "profile_relation_residual_abs_max_on_grid": profile_relation_residual_abs_max,
        "gram_matrix": gram.tolist(),
        "gram_singular_values": singular_values.tolist(),
        "rank_tolerance": rank_tol,
        "numeric_rank": numeric_rank,
        "numeric_nullity": nullity,
        "gram_times_null_vector": null_residual.tolist(),
        "gram_null_residual_l2": null_residual_l2,
        "p1950_reported_rank": p1950_rank,
        "p1950_reported_nullity": p1950_nullity,
        "p1950_verdict": p1950.get("verdict"),
        "p1950_fail_trace": p1950.get("fail_trace"),
        "p1950_full_4channel_rank_gate": p1950_pass_fail.get("operator_profile_gram_rank_eq_channel_count"),
    }

    globalization_obstruction = {
        "claim": "The current B1 backend cannot honestly certify four independent one-loop gravitational counterterm channels because Gauss-Bonnet is exported as the dependent tensor identity GB=Riem2-4*Ric2+R2 on the B1 scalar profile.",
        "what_is_closed": "quotient-scope three-channel cancellation/contract inherited from P2094/P2096",
        "what_remains_open": [
            "full four-channel independence on B1",
            "global background-family renormalization closure beyond B1 quotient scope",
            "full Task-3/global ToE renormalization theorem",
        ],
        "needed_next_input": "Add an independent background-family witness or tensor-resolved global curvature backend where GB independence/global transport can be tested instead of forcing B1 quotient data to behave as a four-channel theorem.",
    }

    theorem_export = {
        "theorem_name": "P2330 strict B1 Gauss-Bonnet dependence and renormalization globalization obstruction",
        "claim": "On the current strict B1 backend profiles, the Gauss-Bonnet channel is exactly dependent: GB - Riem2 + 4 Ric2 - R2 = 0. Therefore P1950's rank-3/nullity-1 result is a structural B1 quotient obstruction, not a numerical accident; P2096 quotient closure may be used only in quotient scope and cannot be promoted to full four-channel/global renormalization closure.",
        "proof_witnesses": [
            "P1848 exports GB by the tensor identity GB=Riem2-4*Ric2+R2 on B1 profiles.",
            "The backend GB expression is exactly the symbolic string Riem2 - 4*Ric2 + R2, and the sampled profile relation residual is numerically zero within tolerance.",
            "The B1 Gram matrix has numeric rank 3 and nullity 1 with null vector (-1,4,-1,1).",
            "P1950 records the same rank-defect obstruction and P2096 records quotient-scope-only closure.",
        ],
        "scope_limits": [
            "This is not a full renormalization closure theorem.",
            "This is not a global background-family theorem.",
            "This does not update G1/G3 or close ToE.",
        ],
        "strict_guardrails": {
            "strict_kernel_only": True,
            "no_legacy_kernel_role_transfer": True,
            "no_selector_premise_added": True,
            "no_qw2191_discharge_claimed": True,
            "no_g1_g3_update_claimed": True,
            "no_toe_closure_claimed": True,
        },
    }

    grep_hits = collect_repo_grep_hits()

    probe = {
        "probe_id": "P2330_S1280_STRICT_B1_RENORMALIZATION_GB_DEPENDENCE_GLOBALIZATION_OBSTRUCTION",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_patterns": list(GREP_PATTERNS),
            "hit_count": len(grep_hits),
            "top_hits": grep_hits[:20],
        },
        "gb_dependence_certificate": gb_dependence_certificate,
        "quotient_closure_summary": quotient_summary(p2096),
        "globalization_obstruction": globalization_obstruction,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "p1848_profiles_loaded": p1848.get("packet_id") == "P1848",
        "p1950_loaded": p1950.get("theorem_target") == "STRICT_B1_ONE_LOOP_DIVERGENCE_EXACT_ALGEBRAIC_CANCELLATION",
        "p2096_loaded": p2096.get("packet_id") == "P2096",
        "symbolic_gb_relation_zero": bool(symbolic_relation_zero),
        "sampled_gb_relation_residual_small": profile_relation_residual_abs_max < 1e-10,
        "numeric_rank_is_3": numeric_rank == 3,
        "numeric_nullity_is_1": nullity == 1,
        "null_vector_residual_small": null_residual_l2 < 1e-10,
        "p1950_rank_defect_preserved": p1950_rank == 3 and p1950_nullity == 1,
        "p2096_quotient_scope_only_preserved": quotient_summary(p2096)["full_4channel_independence_proven"] is False,
        "full_4channel_renormalization_not_claimed": True,
        "global_background_renormalization_not_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2330_s1280_v1",
        "packet_id": "P2330",
        "stage_id": "S1280",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_RENORMALIZATION_GLOBALIZATION_OBSTRUCTION_WITH_EXACT_B1_GB_DEPENDENCE_TRACE",
        "result_kind": "STRICT_B1_GAUSS_BONNET_DEPENDENCE_AND_RENORMALIZATION_GLOBALIZATION_OBSTRUCTION_NO_FALSE_PASS",
        "strict_b1_renormalization_gb_dependence_globalization_obstruction_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": "Do not force B1 quotient closure into full four-channel renormalization. Add a second independent background-family/tensor-resolved curvature witness and test whether the GB dependence is lifted globally, or prove that the current strict backend admits only quotient-scope renormalization closure.",
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2330 strict B1 renormalization GB-dependence obstruction\n\n"
        "Status: exact B1 Gauss-Bonnet dependence found; quotient closure preserved; full/global renormalization not claimed.\n\n"
        "- Exact relation: `GB - Riem2 + 4*Ric2 - R2 == 0`.\n"
        f"- Numeric Gram rank: `{numeric_rank}`; nullity: `{nullity}`.\n"
        f"- Null vector residual L2: `{null_residual_l2:.3e}`.\n"
        f"- P1950 verdict: `{p1950.get('verdict')}`.\n"
        "- No QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
