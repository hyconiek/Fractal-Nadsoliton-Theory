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
OUT = GEN / "p2638_s1588_metric_pushforward_source_viability_and_neural_closure_audit.json"
MD = GEN / "p2638_s1588_metric_pushforward_source_viability_and_neural_closure_audit.md"
REPO_ROOT = REPO

OMEGA = Fraction(1, 4)  # pi units
PHI = Fraction(1, 6)  # pi units
ETA = Fraction(9, 5)
BETA = Fraction(1, 1)
ALPHA_GEO = 4.0 * math.log(2.0)
BETA_TORS = 0.01
PUSH_A = Fraction(4, 3)
PUSH_B = Fraction(-4, 3)
LEGACY_NODES = [2, 5, 8, 11]
INTEGER_DOMAIN = list(range(1, 13))
UV_DOMAIN = [Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1)]

SOURCE_FILES = {
    "P2637_PHASE_NODE_QUOTIENT_EXHAUSTION": GEN / "p2637_s1587_phase_node_quotient_exhaustion_toe_closure_path_audit.json",
    "P2636_BLOCKER_LATTICE": GEN / "p2636_s1586_current_toe_blocker_lattice_full_kernel_decision_audit.json",
    "P2635_TOE_NEURAL_UNIVERSE_SIGNATURE": GEN / "p2635_s1585_toe_neural_universe_empirical_signature_audit.json",
    "P2634_STABILITY_VS_ROLE_COMPLETENESS": GEN / "p2634_s1584_strict_stability_evidence_vs_role_completeness_audit.json",
    "DIAGRAMS_KERNEL_TRANSFORMATION": REPO_ROOT / "DIAGRAMS_KERNEL_TRANSFORMATION.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "metric_pushforward_source_theorem_exported",
    "phase_frequency_node_gauge_certificate_exported",
    "legacy_to_strict_bridge_completion_exported",
    "legacy_role_transfer_revalidated",
    "strict_kernel_full_kernel_claimed",
    "toe_closure_claimed",
    "self_learning_universe_proven",
    "positive_beta_renormalization_source_exported",
    "inverse_hierarchy_role_transfer_exported",
    "selector_source_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "blind_empirical_confirmation_claimed",
]


def sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "*.lean",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO_ROOT,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:90]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: search the intended research content, not ticket IDs.
    patterns = {
        "metric_pushforward_source_content": (
            "metric pushforward|distance-coordinate|renormalized distance|coordinate warp|metric warp|"
            "distance redefinition|topology.*selector|selector dynamics|learned.*distance|fractal metric"
        ),
        "phase_node_role_content": (
            "node/gauge|integer node|zero lattice|phase-frequency|phase node|gauge separation|"
            "cos\\(omega|cos\\(ω|SU\\(3\\).*SU\\(2\\).*U\\(1\\)"
        ),
        "attention_damping_viability_content": (
            "heavy-tailed attention|attention bias|nonlinear damping|fractal compression|inverse hierarchy|"
            "Wilson-loop|distant octave|effective support|information throughput"
        ),
        "neural_universe_empirical_closure_content": (
            "self-learning|samoucząca|neural universe|Energy-Based Model|Graph Neural Network|"
            "positional encoding|CMB|LSS|GW/PTA|frozen-kernel|holdout|preregistration"
        ),
        "toe_blocker_closure_content": (
            "ToE closure|full kernel|blocker lattice|role-transfer|QW-2191|selector source|"
            "role-bearing L_total|positive_beta_renormalization_source"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for metric-pushforward source viability, neural interpretation, and ToE closure blockers",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def source_fingerprints() -> dict[str, Any]:
    return {name: {"path": rel(path), "exists": path.exists(), "sha256": sha256_file(path)} for name, path in SOURCE_FILES.items()}


def frac_str(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def r_push(d: Fraction) -> Fraction:
    return PUSH_A * d + PUSH_B


def phase(d: float) -> float:
    return math.cos(math.pi * (float(OMEGA) * d + float(PHI)))


def strict_kernel(d: float) -> float:
    return phase(d) / (1.0 + float(BETA) * (d ** float(ETA)))


def pushed_strict_kernel(d: Fraction) -> float | None:
    r = r_push(d)
    if r < 0:
        return None
    rf = float(r)
    return phase(rf) / (1.0 + float(BETA) * (rf ** float(ETA)))


def legacy_kernel_amplitude_normalized(d: float) -> float:
    return phase(d) / (1.0 + BETA_TORS * d)


def phase_node_repair_check() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for node in LEGACY_NODES:
        d = Fraction(node)
        r = r_push(d)
        rows.append(
            {
                "legacy_node_d": node,
                "pushed_strict_distance_r_d": frac_str(r),
                "cos_phase_after_pushforward": phase(float(r)),
                "is_exact_zero_with_tolerance_1e_12": abs(phase(float(r))) < 1.0e-12,
            }
        )
    return {
        "candidate_pushforward": f"r(d)=({frac_str(PUSH_A)})*d+({frac_str(PUSH_B)})",
        "node_rows": rows,
        "all_legacy_nodes_become_strict_phase_zeros": all(row["is_exact_zero_with_tolerance_1e_12"] for row in rows),
        "source_status": "constructive_math_repair_only_not_topology_selector_source",
    }


def domain_admissibility_check() -> dict[str, Any]:
    uv_rows: list[dict[str, Any]] = []
    for d in UV_DOMAIN:
        r = r_push(d)
        uv_rows.append(
            {
                "d": frac_str(d),
                "r_d": frac_str(r),
                "real_fractional_power_damping_defined": r >= 0,
                "comment": "negative r(d) makes r(d)^(9/5) non-real on the principal real branch" if r < 0 else "real branch ok",
            }
        )
    integer_rows: list[dict[str, Any]] = []
    for d_int in INTEGER_DOMAIN:
        d = Fraction(d_int)
        r = r_push(d)
        integer_rows.append(
            {
                "d": d_int,
                "r_d": frac_str(r),
                "collapses_positive_distance_to_zero": d_int > 0 and r == 0,
                "real_fractional_power_damping_defined": r >= 0,
            }
        )
    return {
        "uv_rows": uv_rows,
        "integer_rows": integer_rows,
        "is_real_on_open_uv_interval_0_to_1": all(row["real_fractional_power_damping_defined"] for row in uv_rows if row["d"] != "1"),
        "collapses_d_equals_1_to_strict_zero_distance": r_push(Fraction(1)) == 0,
        "requires_domain_cut_or_shifted_uv_origin": True,
        "admissibility_verdict": "fails_as_global_positive_distance_metric; admissible only as a chart-local integer-lattice map with extra UV/domain theorem",
    }


def damping_and_attention_distortion() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for d_int in INTEGER_DOMAIN:
        d = Fraction(d_int)
        pushed = pushed_strict_kernel(d)
        identity = strict_kernel(float(d_int))
        legacy = legacy_kernel_amplitude_normalized(float(d_int))
        rows.append(
            {
                "d": d_int,
                "r_d": frac_str(r_push(d)),
                "legacy_amplitude_normalized_abs_kernel": abs(legacy),
                "strict_identity_abs_kernel": abs(identity),
                "strict_pushed_abs_kernel": None if pushed is None else abs(pushed),
            }
        )
    row1 = rows[0]
    row7 = next(row for row in rows if row["d"] == 7)
    row10 = next(row for row in rows if row["d"] == 10)
    pushed_ratio_7_1 = row7["strict_pushed_abs_kernel"] / row1["strict_pushed_abs_kernel"]
    identity_ratio_7_1 = row7["strict_identity_abs_kernel"] / row1["strict_identity_abs_kernel"]
    legacy_ratio_7_1 = row7["legacy_amplitude_normalized_abs_kernel"] / row1["legacy_amplitude_normalized_abs_kernel"]
    pushed_ratio_10_1 = row10["strict_pushed_abs_kernel"] / row1["strict_pushed_abs_kernel"]

    weights = [row["strict_pushed_abs_kernel"] ** 2 for row in rows if row["strict_pushed_abs_kernel"] is not None]
    total = sum(weights)
    probs = [w / total for w in weights]
    entropy = -sum(p * math.log(p) for p in probs if p > 0.0)
    effective_support = math.exp(entropy)
    return {
        "domain": INTEGER_DOMAIN,
        "rows": rows,
        "ratio_abs_k7_over_abs_k1": {
            "legacy_amplitude_normalized": legacy_ratio_7_1,
            "strict_identity": identity_ratio_7_1,
            "strict_after_metric_pushforward": pushed_ratio_7_1,
            "pushforward_preserves_legacy_inverse_hierarchy_above_one": pushed_ratio_7_1 > 1.0,
        },
        "ratio_abs_k10_over_abs_k1_after_metric_pushforward": pushed_ratio_10_1,
        "normalized_attention_proxy_after_pushforward": {
            "entropy_nats": entropy,
            "effective_support_points": effective_support,
            "max_probability": max(probs),
            "argmax_d": rows[probs.index(max(probs))]["d"],
        },
        "distortion_verdict": "node repair is exact, but pushed strict damping becomes near-neighbor dominated and still does not recover the legacy inverse-hierarchy role on the same integer domain",
    }


def closure_decision(phase_repair: dict[str, Any], domain: dict[str, Any], damping: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "phase_nodes_repaired_exactly": phase_repair["all_legacy_nodes_become_strict_phase_zeros"],
        "global_positive_distance_metric_on_uv_domain": domain["is_real_on_open_uv_interval_0_to_1"] and not domain["collapses_d_equals_1_to_strict_zero_distance"],
        "inverse_hierarchy_role_preserved_after_pushforward": damping["ratio_abs_k7_over_abs_k1"]["pushforward_preserves_legacy_inverse_hierarchy_above_one"],
        "independent_topology_selector_source_present": False,
        "beta_uv_and_ltotal_sources_present": False,
        "blind_empirical_holdout_confirmed": False,
    }
    return {
        "gates": gates,
        "metric_pushforward_promoted_to_bridge_theorem": all(gates.values()),
        "full_kernel_now": False,
        "professorial_verdict": (
            "The affine pushforward is a useful exact mathematical repair for the node/gauge list, but it is not yet a physical completion map. "
            "It fails as a global positive UV distance metric, collapses d=1 to strict distance zero, and does not preserve the legacy inverse-hierarchy role on the same domain. "
            "Thus it narrows the proof obligation instead of closing the ToE."
        ),
    }


def professorial_neural_toe_path(decision: dict[str, Any]) -> dict[str, Any]:
    return {
        "strongest_toe_symptoms_now": [
            {
                "symptom": "single stable strict kernel",
                "professorial_reading": "A single robust kernel with repeated stability evidence is ToE-like because it behaves like a universal architecture rather than a sector-by-sector fit.",
                "how_to_check_with_modern_physics": "freeze omega, phi, beta, eta and test no-retune predictions in CMB/LSS correlation tails, GW/PTA stochastic backgrounds, and cross-sector spectra.",
            },
            {
                "symptom": "neural positional plus heavy-tailed attention structure",
                "professorial_reading": "cos(omega*d+phi) is a positional/resonance channel and 1/(1+beta*d^eta) is a power-law attention bias, so the kernel resembles an energy-based graph neural layer on a fractal graph.",
                "how_to_check_with_modern_physics": "compare frozen heavy-tailed attention to exponential and flexible spline baselines under preregistered likelihood/Bayes-factor criteria.",
            },
            {
                "symptom": "constructive but source-guarded node repair",
                "professorial_reading": "The exact affine repair shows the mismatch is mathematically structured, not random; however a self-learning universe may learn a coordinate only if the theory derives the coordinate from its variational/topological dynamics.",
                "how_to_check_with_modern_physics": "derive r(d)=4/3*(d-1) from topology/selector equations, then test whether the induced phase nodes align with observed gauge/generation separation without retuning.",
            },
        ],
        "closure_path_after_p2638": [
            "Either derive a topology/selector source for a domain-safe metric coordinate, or demote the legacy integer node/gauge role.",
            "Derive a dimensionless beta/UV normalization invariant instead of using beta=1 by convention.",
            "Recompute role-transfer under the sourced metric, especially inverse hierarchy and Wilson-loop distant-octave claims.",
            "Resolve QW-2191 by a real strict-core selector source or keep the closure explicitly non-strict/axiom-augmented.",
            "Only after bridge and role-transfer pass, build role-bearing L_total and preregister frozen-kernel empirical tests.",
        ],
        "next_honest_step": (
            "Do not add another global symptom score.  The next proof-grade computation should search the repo's topology/selector certificates for a domain-safe metric-coordinate law and test whether any law derives an affine map equivalent to r(d)=4/3*(d-1) while keeping r(d)>0 on the physical UV domain and preserving the inverse-hierarchy role."
        ),
        "full_kernel_answer": "No: strict remains a strong ToE-like working kernel, not a full kernel, until the metric source, beta source, inverse-hierarchy role-transfer, selector, L_total, and blind empirical blockers close.",
        "negative_exports_respected": all(value is False for value in decision.get("negative_export_flags", {}).values()) if "negative_export_flags" in decision else True,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    phase = payload["phase_node_repair_check"]
    domain = payload["domain_admissibility_check"]
    damping = payload["damping_and_attention_distortion"]
    decision = payload["closure_decision"]
    neural = payload["professorial_neural_toe_path"]
    lines = [
        "# P2638/S1588 metric-pushforward source viability and neural closure audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps metric-pushforward/source, phase-node, attention/damping, neural-universe empirical, and ToE-blocker content before adding the viability test.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Exact repair versus physical viability",
        "",
        f"Candidate: `{phase['candidate_pushforward']}`.",
        f"All legacy nodes become strict phase zeros? `{phase['all_legacy_nodes_become_strict_phase_zeros']}`.",
        f"Domain verdict: `{domain['admissibility_verdict']}`.",
        f"Legacy |K(7)|/|K(1)|: `{damping['ratio_abs_k7_over_abs_k1']['legacy_amplitude_normalized']:.10f}`.",
        f"Strict identity |K(7)|/|K(1)|: `{damping['ratio_abs_k7_over_abs_k1']['strict_identity']:.10f}`.",
        f"Strict pushed |K(7)|/|K(1)|: `{damping['ratio_abs_k7_over_abs_k1']['strict_after_metric_pushforward']:.10f}`.",
        f"Pushed attention effective support on d=1..12: `{damping['normalized_attention_proxy_after_pushforward']['effective_support_points']:.6f}` points, argmax d=`{damping['normalized_attention_proxy_after_pushforward']['argmax_d']}`.",
        "",
        "## Closure decision",
        "",
        decision["professorial_verdict"],
        "",
        f"Promote metric pushforward to bridge theorem? `{decision['metric_pushforward_promoted_to_bridge_theorem']}`.",
        f"Full kernel now? `{decision['full_kernel_now']}`.",
        "",
        "## Professorial neural/ToE reading",
        "",
    ])
    for row in neural["strongest_toe_symptoms_now"]:
        lines.append(f"- **{row['symptom']}**: {row['professorial_reading']} Check: {row['how_to_check_with_modern_physics']}")
    lines.extend([
        "",
        "## Closure path after P2638",
        "",
    ])
    for idx, item in enumerate(neural["closure_path_after_p2638"], start=1):
        lines.append(f"{idx}. {item}")
    lines.extend([
        "",
        "## Recommended next honest step",
        "",
        neural["next_honest_step"],
        "",
        neural["full_kernel_answer"],
        "",
        "No ToE closure, full-kernel finality, bridge completion, role-transfer, selector-source discharge, positive beta source, blind empirical confirmation, or role-bearing `L_total` is claimed.",
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    phase = phase_node_repair_check()
    domain = domain_admissibility_check()
    damping = damping_and_attention_distortion()
    decision = closure_decision(phase, domain, damping)
    payload: dict[str, Any] = {
        "status": "P2638_METRIC_PUSHFORWARD_SOURCE_VIABILITY_AND_NEURAL_CLOSURE_AUDIT_NO_PROMOTION",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_fingerprints": source_fingerprints(),
        "source_payloads_present": {name: load_json(path).get("missing", False) is False for name, path in SOURCE_FILES.items() if path.suffix == ".json"},
        "phase_node_repair_check": phase,
        "domain_admissibility_check": domain,
        "damping_and_attention_distortion": damping,
        "closure_decision": decision,
        "professorial_neural_toe_path": professorial_neural_toe_path(decision),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)

    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2638/S1588 metric-pushforward source viability guard",
        "\n## P2638/S1588 metric-pushforward source viability guard\n\n"
        "`P2638/S1588` audits the P2637 affine node repair as a physical metric candidate.  The map `r(d)=4/3*(d-1)` exactly sends the legacy node list `2,5,8,11` to strict phase zeros, but it fails as a global positive UV distance metric because it is negative for `0<d<1` and collapses `d=1` to strict distance zero.  On `d=1..12` it also does not restore the legacy inverse-hierarchy role: the pushed strict `|K(7)|/|K(1)|` remains below one.  Therefore the repair is a source-guarded mathematical candidate, not bridge completion, role-transfer, full-kernel finality, or ToE closure.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2638/S1588 metric-pushforward Ltotal viability guard",
        "\n## P2638/S1588 metric-pushforward Ltotal viability guard\n\n"
        "`P2638/S1588` prevents promoting the affine node repair into `L_total`: a role-bearing dynamics would first need a topology/selector-derived, domain-safe metric coordinate that keeps the strict damping real on the UV domain and preserves any transferred inverse-hierarchy role.  Neural self-learning language remains a test-generating analogy until that source and blind frozen-kernel empirical checks exist.\n",
    )
    return payload


if __name__ == "__main__":
    main()
