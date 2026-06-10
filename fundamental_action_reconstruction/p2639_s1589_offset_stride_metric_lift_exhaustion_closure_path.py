#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2639_s1589_offset_stride_metric_lift_exhaustion_closure_path.json"
MD = GEN / "p2639_s1589_offset_stride_metric_lift_exhaustion_closure_path.md"
REPO_ROOT = REPO

OMEGA = Fraction(1, 4)
PHI = Fraction(1, 6)
ETA = Fraction(9, 5)
BETA = Fraction(1, 1)
BETA_TORS = 0.01
LEGACY_NODES = [2, 5, 8, 11]
K0_RANGE = range(-2, 13)
STRIDE_RANGE = range(1, 7)
INTEGER_DOMAIN = list(range(1, 13))
UV_GRID = [Fraction(n, 12) for n in range(0, 13)]

SOURCE_FILES = {
    "P2637_PHASE_NODE_QUOTIENT_EXHAUSTION": GEN / "p2637_s1587_phase_node_quotient_exhaustion_toe_closure_path_audit.json",
    "P2638_METRIC_PUSHFORWARD_VIABILITY": GEN / "p2638_s1588_metric_pushforward_source_viability_and_neural_closure_audit.json",
    "P2636_BLOCKER_LATTICE": GEN / "p2636_s1586_current_toe_blocker_lattice_full_kernel_decision_audit.json",
    "DIAGRAMS_KERNEL_TRANSFORMATION": REPO_ROOT / "DIAGRAMS_KERNEL_TRANSFORMATION.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "offset_stride_metric_lift_source_exported",
    "phase_frequency_node_gauge_certificate_exported",
    "legacy_to_strict_bridge_completion_exported",
    "legacy_role_transfer_revalidated",
    "strict_kernel_full_kernel_claimed",
    "toe_closure_claimed",
    "positive_beta_renormalization_source_exported",
    "inverse_hierarchy_role_transfer_exported",
    "selector_source_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "blind_empirical_confirmation_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "*.lean",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: offset/stride/lift/source contents, not packet numbers.
    patterns = {
        "offset_stride_zero_lattice_content": (
            "zero lattice|integer node|node/gauge|offset|stride|lift|quotient|reindex|"
            "phase-frequency|gauge separation|cos\\(omega|cos\\(ω"
        ),
        "metric_coordinate_source_content": (
            "metric pushforward|distance-coordinate|domain-safe metric|coordinate warp|renormalized distance|"
            "topology.*selector|selector dynamics|fractal metric|source theorem"
        ),
        "inverse_hierarchy_attention_content": (
            "inverse hierarchy|distant octave|Wilson-loop|heavy-tailed attention|attention bias|"
            "effective support|fractal compression|nonlinear damping"
        ),
        "toe_closure_blocker_content": (
            "ToE closure|full kernel|bridge completion|role-transfer|QW-2191|selector source|"
            "role-bearing L_total|blind empirical|frozen-kernel"
        ),
        "neural_empirical_path_content": (
            "self-learning|samoucząca|neural universe|Energy-Based Model|positional encoding|"
            "CMB|LSS|GW/PTA|holdout|preregistration"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for offset/stride metric-lift closure search", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def frac_str(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def zero(index: int) -> Fraction:
    return Fraction(4, 3) + 4 * index


def phase(x: float) -> float:
    return math.cos(math.pi * (float(OMEGA) * x + float(PHI)))


def strict_abs(x: float) -> float:
    return abs(phase(x) / (1.0 + float(BETA) * (x ** float(ETA))))


def legacy_abs(d: int) -> float:
    return abs(phase(float(d)) / (1.0 + BETA_TORS * d))


def affine_for(k0: int, stride: int) -> tuple[Fraction, Fraction]:
    # Map legacy nodes 2,5,8,11 to zero indices k0, k0+stride, k0+2*stride, k0+3*stride.
    a = Fraction(4 * stride, 3)
    b = zero(k0) - a * LEGACY_NODES[0]
    return a, b


def evaluate_candidate(k0: int, stride: int) -> dict[str, Any]:
    a, b = affine_for(k0, stride)
    def r(d: Fraction) -> Fraction:
        return a * d + b

    node_rows = []
    for offset, node in enumerate(LEGACY_NODES):
        target_index = k0 + offset * stride
        rd = r(Fraction(node))
        node_rows.append({
            "legacy_node": node,
            "target_zero_index": target_index,
            "r_d": frac_str(rd),
            "target_zero": frac_str(zero(target_index)),
            "exact_match": rd == zero(target_index),
            "cos_residual": phase(float(rd)),
        })
    uv_positive = all(r(d) > 0 for d in UV_GRID if d > 0)
    r1 = r(Fraction(1))
    r7 = r(Fraction(7))
    strict_ratio = strict_abs(float(r7)) / strict_abs(float(r1)) if r1 > 0 and r7 > 0 else None
    legacy_ratio = legacy_abs(7) / legacy_abs(1)
    r_values = [r(Fraction(d)) for d in INTEGER_DOMAIN]
    weights = [strict_abs(float(rv)) ** 2 for rv in r_values if rv > 0]
    total = sum(weights)
    probs = [w / total for w in weights] if total else []
    entropy = -sum(p * math.log(p) for p in probs if p > 0)
    return {
        "k0": k0,
        "stride": stride,
        "map": f"r(d)=({frac_str(a)})*d+({frac_str(b)})",
        "a": frac_str(a),
        "b": frac_str(b),
        "node_rows": node_rows,
        "exact_node_lift": all(row["exact_match"] for row in node_rows),
        "uv_domain_positive_on_grid_0_1_step_1_12": uv_positive,
        "r_1": frac_str(r1),
        "r_7": frac_str(r7),
        "legacy_abs_k7_over_k1": legacy_ratio,
        "strict_lifted_abs_k7_over_k1": strict_ratio,
        "inverse_hierarchy_ratio_above_one": strict_ratio is not None and strict_ratio > 1.0,
        "effective_support_points_d1_to_d12": math.exp(entropy) if probs else 0.0,
        "attention_argmax_d": INTEGER_DOMAIN[probs.index(max(probs))] if probs else None,
        "nonlocality_cost_r1": float(r1) if r1 > 0 else None,
    }


def exhaustion() -> dict[str, Any]:
    rows = [evaluate_candidate(k0, stride) for stride in STRIDE_RANGE for k0 in K0_RANGE]
    exact = [row for row in rows if row["exact_node_lift"]]
    uv_safe = [row for row in exact if row["uv_domain_positive_on_grid_0_1_step_1_12"]]
    role_like = [row for row in uv_safe if row["inverse_hierarchy_ratio_above_one"]]
    best_role = max(role_like, key=lambda row: row["strict_lifted_abs_k7_over_k1"], default=None)
    nearest_local_uv_safe = min(uv_safe, key=lambda row: (row["stride"], row["k0"]), default=None)
    return {
        "audited_k0_range": [min(K0_RANGE), max(K0_RANGE)],
        "audited_stride_range": [min(STRIDE_RANGE), max(STRIDE_RANGE)],
        "total_candidates": len(rows),
        "exact_lift_count": len(exact),
        "uv_safe_exact_lift_count": len(uv_safe),
        "uv_safe_and_inverse_hierarchy_count": len(role_like),
        "nearest_local_uv_safe_exact_lift": nearest_local_uv_safe,
        "best_inverse_hierarchy_uv_safe_lift": best_role,
        "role_like_candidates": role_like[:12],
        "classification": "OFFSET_STRIDE_LIFTS_EXIST_BUT_REQUIRE_UNSOURCED_QUOTIENT_OFFSET_OR_STRIDE_AND_DO_NOT_CLOSE_BRIDGE",
    }


def closure(exh: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "some_domain_safe_exact_lift_exists": exh["uv_safe_exact_lift_count"] > 0,
        "some_domain_safe_inverse_hierarchy_lift_exists": exh["uv_safe_and_inverse_hierarchy_count"] > 0,
        "canonical_offset_stride_source_theorem_exists": False,
        "nonlocal_zero_lift_physically_licensed": False,
        "beta_selector_ltotal_empirical_blockers_closed": False,
    }
    return {
        "gates": gates,
        "promote_any_lift_to_bridge_completion": all(gates.values()),
        "full_kernel_now": False,
        "professorial_verdict": (
            "P2639 corrects the overly local reading of the first affine repair: offset/stride lifts can be exact and UV-positive, and some can even make the K7/K1 proxy exceed one. "
            "However those successes occur only after choosing an unsourced zero-lattice offset/stride, often a nonlocal lift into later strict zeros.  This is a sharper candidate class, not a ToE closure theorem."
        ),
        "next_honest_step": (
            "Search topology/selector artifacts for a canonical source of the offset k0 and stride m.  If no such source exists, demote node/gauge transfer to an approximate or rejected legacy role; if it exists, rerun beta, inverse-hierarchy, QW-2191, L_total, and blind frozen-kernel gates under that sourced lift."
        ),
    }


def professorial_path() -> list[dict[str, str]]:
    return [
        {
            "rank": "1",
            "task": "canonical offset/stride source theorem",
            "why": "The finite search found mathematical lifts; ToE closure needs the nadsoliton dynamics to choose one without retuning.",
        },
        {
            "rank": "2",
            "task": "role-transfer rerun under the sourced lift",
            "why": "Node repair, inverse hierarchy, beta normalization, and alpha_geo/beta_tors roles must survive together, not one at a time.",
        },
        {
            "rank": "3",
            "task": "neural/frozen-kernel empirical holdout",
            "why": "A self-learning-universe reading becomes physics only when the frozen attention kernel predicts CMB/LSS, GW/PTA, or cross-sector data better than baselines without retuning.",
        },
    ]


def write_markdown(payload: dict[str, Any]) -> None:
    exh = payload["offset_stride_exhaustion"]
    cl = payload["closure_decision"]
    lines = [
        "# P2639/S1589 offset-stride metric-lift exhaustion closure path",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps zero-lattice offset/stride, metric-coordinate source, inverse-hierarchy attention, ToE blockers, and neural empirical path content before adding the finite search.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Exhaustion result",
        "",
        f"Audited candidates: `{exh['total_candidates']}` over k0=`{exh['audited_k0_range']}` and stride=`{exh['audited_stride_range']}`.",
        f"UV-safe exact node lifts: `{exh['uv_safe_exact_lift_count']}`.",
        f"UV-safe exact lifts with |K(7)|/|K(1)| > 1: `{exh['uv_safe_and_inverse_hierarchy_count']}`.",
        f"Nearest UV-safe exact lift: `{exh['nearest_local_uv_safe_exact_lift']['map']}` with k0=`{exh['nearest_local_uv_safe_exact_lift']['k0']}`, stride=`{exh['nearest_local_uv_safe_exact_lift']['stride']}`.",
    ])
    best = exh["best_inverse_hierarchy_uv_safe_lift"]
    if best:
        lines.append(f"Best inverse-hierarchy UV-safe lift in the audited box: `{best['map']}`, k0=`{best['k0']}`, stride=`{best['stride']}`, lifted ratio=`{best['strict_lifted_abs_k7_over_k1']:.10f}`.")
    lines.extend([
        "",
        "## Closure decision",
        "",
        cl["professorial_verdict"],
        "",
        f"Promote any lift to bridge completion? `{cl['promote_any_lift_to_bridge_completion']}`.",
        f"Full kernel now? `{cl['full_kernel_now']}`.",
        "",
        "## Professorial closure path",
        "",
    ])
    for row in payload["professorial_closure_path"]:
        lines.append(f"{row['rank']}. **{row['task']}** — {row['why']}")
    lines.extend([
        "",
        "## Recommended next honest step",
        "",
        cl["next_honest_step"],
        "",
        "No ToE closure, full-kernel finality, bridge completion, role-transfer, selector-source discharge, positive beta source, blind empirical confirmation, or role-bearing `L_total` is claimed.",
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    exh = exhaustion()
    cl = closure(exh)
    payload: dict[str, Any] = {
        "status": "P2639_OFFSET_STRIDE_METRIC_LIFT_EXHAUSTION_NO_PROMOTION",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_fingerprints": {name: {"path": rel(path), "exists": path.exists(), "sha256": sha256_file(path)} for name, path in SOURCE_FILES.items()},
        "source_payloads_present": {name: not load_json(path).get("missing", False) for name, path in SOURCE_FILES.items() if path.suffix == ".json"},
        "offset_stride_exhaustion": exh,
        "closure_decision": cl,
        "professorial_closure_path": professorial_path(),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2639/S1589 offset-stride metric-lift guard",
        "\n## P2639/S1589 offset-stride metric-lift guard\n\n"
        "`P2639/S1589` widens the P2637/P2638 affine repair from the first zero-lattice block to offset/stride lifts.  The finite search finds UV-safe exact node lifts and even some lifts with `|K(7)|/|K(1)| > 1`, but only by choosing an unsourced zero-lattice offset/stride.  Therefore these lifts are sharper bridge candidates, not bridge completion: a topology/selector theorem must canonically choose the offset and stride before node/gauge role-transfer, beta normalization, inverse hierarchy, `QW-2191`, role-bearing `L_total`, or ToE closure can reopen.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2639/S1589 offset-stride lift Ltotal guard",
        "\n## P2639/S1589 offset-stride lift Ltotal guard\n\n"
        "`P2639/S1589` does not re-enable `L_total`: offset/stride metric lifts are admissible mathematical candidates only until a strict topology/selector source chooses one and the full role-transfer matrix is rerun under that sourced lift.\n",
    )
    return payload


if __name__ == "__main__":
    main()
