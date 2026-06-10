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
OUT = GEN / "p2637_s1587_phase_node_quotient_exhaustion_toe_closure_path_audit.json"
MD = GEN / "p2637_s1587_phase_node_quotient_exhaustion_toe_closure_path_audit.md"
REPO_ROOT = REPO

OMEGA = Fraction(1, 4)  # in pi units: omega = pi/4
PHI = Fraction(1, 6)    # in pi units: phi = pi/6
DECLARED_LEGACY_NODES = [Fraction(2), Fraction(5), Fraction(8), Fraction(11)]
ZERO_INDICES = [0, 1, 2, 3]

SOURCE_FILES = {
    "P2633_DIAGRAM_GROUNDED_RETENTION": GEN / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.json",
    "P2634_STABILITY_VS_ROLE_COMPLETENESS": GEN / "p2634_s1584_strict_stability_evidence_vs_role_completeness_audit.json",
    "P2635_TOE_NEURAL_UNIVERSE_SIGNATURE": GEN / "p2635_s1585_toe_neural_universe_empirical_signature_audit.json",
    "P2636_BLOCKER_LATTICE": GEN / "p2636_s1586_current_toe_blocker_lattice_full_kernel_decision_audit.json",
    "DIAGRAMS_KERNEL_TRANSFORMATION": REPO_ROOT / "DIAGRAMS_KERNEL_TRANSFORMATION.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "toe_closure_claimed",
    "strict_kernel_full_kernel_claimed",
    "phase_frequency_node_gauge_certificate_exported",
    "metric_warp_source_theorem_exported",
    "legacy_to_strict_bridge_completion_exported",
    "legacy_role_transfer_revalidated",
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
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: search research content before adding this certificate.
    patterns = {
        "phase_node_zero_lattice_content": (
            "node/gauge|integer node|zero lattice|cos\\(omega|cos\\(ω|omega=π/4|omega = pi/4|"
            "phase-frequency|phase frequency|gauge separation|SU\\(3\\).*SU\\(2\\).*U\\(1\\)"
        ),
        "quotient_reindexing_metric_warp_content": (
            "quotient|reindex|re-index|coordinate warp|metric warp|distance redefinition|"
            "renormalized distance|pushforward|pullback|selector map"
        ),
        "neural_positional_attention_content": (
            "positional encoding|attention bias|heavy-tailed attention|Transformer|ALiBi|"
            "Energy-Based Model|self-learning|samoucząca|Graph Neural Network"
        ),
        "toe_symptom_and_blocker_content": (
            "ToE-like|ToE symptoms|full kernel|complete kernel|blocker lattice|"
            "bridge completion|role-transfer|QW-2191|role-bearing L_total|blind empirical"
        ),
        "modern_physics_closure_tests_content": (
            "CMB|LSS|large scale structure|GW/PTA|gravitational wave|frozen-kernel|holdout|"
            "preregistration|phase node|Wilson-loop|inverse hierarchy"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for phase/node quotient exhaustion and ToE closure path; not ticket-number search",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def source_fingerprints() -> dict[str, Any]:
    return {name: {"path": rel(path), "exists": path.exists(), "sha256": sha256_file(path)} for name, path in SOURCE_FILES.items()}


def exact_zero(k: int) -> Fraction:
    # omega*d + phi = 1/2 + k in pi units.
    return (Fraction(1, 2) + Fraction(k) - PHI) / OMEGA


def cos_residual_at_strict_input(d: Fraction) -> float:
    return math.cos(math.pi * float(OMEGA * d + PHI))


def frac_str(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def regression_affine(xs: list[Fraction], ys: list[Fraction]) -> tuple[Fraction, Fraction]:
    n = Fraction(len(xs))
    xbar = sum(xs, Fraction(0)) / n
    ybar = sum(ys, Fraction(0)) / n
    denom = sum((x - xbar) ** 2 for x in xs)
    if denom == 0:
        raise ValueError("degenerate affine fit")
    a = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / denom
    b = ybar - a * xbar
    return a, b


def max_abs_fraction(values: list[Fraction]) -> Fraction:
    return max((abs(v) for v in values), default=Fraction(0))


def phase_node_quotient_exhaustion() -> dict[str, Any]:
    zeros = [exact_zero(k) for k in ZERO_INDICES]
    declared = DECLARED_LEGACY_NODES

    identity_rows = [
        {
            "declared_node": frac_str(d),
            "strict_phase_cos_value": cos_residual_at_strict_input(d),
            "is_exact_formula_zero": abs(cos_residual_at_strict_input(d)) < 1.0e-12,
        }
        for d in declared
    ]

    identity_distance_residuals = [d - z for d, z in zip(declared, zeros)]
    pure_translation_first = zeros[0] - declared[0]
    pure_translation_residuals = [(d + pure_translation_first) - z for d, z in zip(declared, zeros)]
    pure_translation_lsq = sum((z - d) for d, z in zip(declared, zeros)) / Fraction(len(declared))
    pure_translation_lsq_residuals = [(d + pure_translation_lsq) - z for d, z in zip(declared, zeros)]

    pure_scale_a = sum(d * z for d, z in zip(declared, zeros)) / sum(d * d for d in declared)
    pure_scale_residuals = [pure_scale_a * d - z for d, z in zip(declared, zeros)]

    affine_a, affine_b = regression_affine(declared, zeros)
    affine_residuals = [affine_a * d + affine_b - z for d, z in zip(declared, zeros)]
    omega_eff = OMEGA * affine_a
    phi_eff = PHI + OMEGA * affine_b

    declared_spacing = [declared[i + 1] - declared[i] for i in range(len(declared) - 1)]
    zero_spacing = [zeros[i + 1] - zeros[i] for i in range(len(zeros) - 1)]

    classes = [
        {
            "map_class": "identity_same_distance_coordinate",
            "definition": "r(d)=d",
            "max_distance_residual_to_exact_zero_lattice": frac_str(max_abs_fraction(identity_distance_residuals)),
            "passes_exact_node_certificate": all(row["is_exact_formula_zero"] for row in identity_rows),
            "interpretation": "Fails: declared integer nodes are not zeros of cos(pi*d/4+pi/6).",
        },
        {
            "map_class": "pure_translation_fixed_by_first_node",
            "definition": f"r(d)=d+({frac_str(pure_translation_first)})",
            "max_distance_residual_to_exact_zero_lattice": frac_str(max_abs_fraction(pure_translation_residuals)),
            "passes_exact_node_certificate": max_abs_fraction(pure_translation_residuals) == 0,
            "interpretation": "Fails because the declared node spacing is 3 while the strict zero-lattice spacing is 4.",
        },
        {
            "map_class": "least_squares_pure_translation",
            "definition": f"r(d)=d+({frac_str(pure_translation_lsq)})",
            "max_distance_residual_to_exact_zero_lattice": frac_str(max_abs_fraction(pure_translation_lsq_residuals)),
            "passes_exact_node_certificate": max_abs_fraction(pure_translation_lsq_residuals) == 0,
            "interpretation": "Fails; translation alone cannot change the lattice period.",
        },
        {
            "map_class": "least_squares_pure_scale_about_origin",
            "definition": f"r(d)=({frac_str(pure_scale_a)})*d",
            "max_distance_residual_to_exact_zero_lattice": frac_str(max_abs_fraction(pure_scale_residuals)),
            "passes_exact_node_certificate": max_abs_fraction(pure_scale_residuals) == 0,
            "interpretation": "Fails; scaling about the origin cannot simultaneously align phase and period.",
        },
        {
            "map_class": "monotone_affine_metric_pushforward",
            "definition": f"r(d)=({frac_str(affine_a)})*d+({frac_str(affine_b)})",
            "max_distance_residual_to_exact_zero_lattice": frac_str(max_abs_fraction(affine_residuals)),
            "passes_exact_node_certificate": max_abs_fraction(affine_residuals) == 0 and affine_a > 0,
            "induced_phase_in_legacy_coordinate": f"cos(({frac_str(omega_eff)})*pi*d + ({frac_str(phi_eff)})*pi)",
            "interpretation": (
                "Constructive but non-silent: exact alignment exists only after a nontrivial distance-coordinate pushforward. "
                "It changes the effective phase parameters in the legacy coordinate, so it needs an independent metric-warp/selector source theorem."
            ),
        },
    ]

    exact_constructive = next(row for row in classes if row["map_class"] == "monotone_affine_metric_pushforward")
    silent_classes = [row for row in classes if row["map_class"] != "monotone_affine_metric_pushforward"]
    return {
        "strict_phase_parameters_pi_units": {"omega": frac_str(OMEGA), "phi": frac_str(PHI)},
        "declared_legacy_nodes": [frac_str(d) for d in declared],
        "exact_zero_lattice_formula": "d_k = 4/3 + 4*k for k in Z",
        "audited_zero_lattice": [frac_str(z) for z in zeros],
        "declared_spacing": [frac_str(s) for s in declared_spacing],
        "strict_zero_spacing": [frac_str(s) for s in zero_spacing],
        "identity_phase_rows": identity_rows,
        "exhausted_map_classes": classes,
        "silent_identity_or_simple_reindexing_passes": any(row["passes_exact_node_certificate"] for row in silent_classes),
        "constructive_affine_metric_pushforward_exists": exact_constructive["passes_exact_node_certificate"],
        "constructive_map_requires_new_source_theorem": True,
        "metric_pushforward_obligation": {
            "needed_theorem": "derive r(d)=4/3*(d-1) from nadsoliton topology/selector dynamics before using the legacy integer node/gauge role",
            "why": "Without this theorem the affine map is an exact mathematical repair but can still be a post-hoc coordinate warp.",
            "neural_network_reading": "A Transformer-like positional channel may use a learned/renormalized distance coordinate, but a ToE kernel must derive that coordinate rather than fit it.",
        },
    }


def professorial_toe_closure_path(phase_cert: dict[str, Any], src: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2636 = src["P2636_BLOCKER_LATTICE"]
    p2635 = src["P2635_TOE_NEURAL_UNIVERSE_SIGNATURE"]
    stability_positive = p2636.get("full_kernel_decision", {}).get("gates", {}).get("strict_kernel_stability_positive", False)
    toe_symptoms = p2636.get("full_kernel_decision", {}).get("gates", {}).get("toe_symptoms_present", False)
    score = p2635.get("weighted_toe_signature_score", {}).get("toe_likeness_score_0_to_1_not_probability", None)
    blockers = p2636.get("current_blocker_lattice", [])

    path = [
        {
            "rank": 1,
            "closure_task": "phase-frequency/node metric-pushforward source theorem",
            "current_status": "constructive affine map found, but source theorem missing",
            "computable_exit_condition": "derive r(d)=4/3*(d-1) or reject the legacy integer node/gauge role",
        },
        {
            "rank": 2,
            "closure_task": "positive beta / UV normalization source",
            "current_status": "not discharged by neural information-flux or Z_beta normalization audits",
            "computable_exit_condition": "derive a dimensionless beta invariant from nadsoliton dynamics before setting beta=1",
        },
        {
            "rank": 3,
            "closure_task": "inverse-hierarchy role-transfer theorem",
            "current_status": "strict compression is stable but does not preserve legacy distant-octave ratios on the same d-grid",
            "computable_exit_condition": "state the distance/domain measure under which Wilson-loop distant-octave coupling is preserved, or downgrade the role",
        },
        {
            "rank": 4,
            "closure_task": "strict-core selector/source obstruction QW-2191",
            "current_status": "still a real selector obstruction",
            "computable_exit_condition": "export a non-axiomatic symmetry-breaking/selector premise or keep closure non-strict",
        },
        {
            "rank": 5,
            "closure_task": "role-bearing L_total dynamics",
            "current_status": "strict kernel has variational/stability evidence, not a full role-bearing L_total source",
            "computable_exit_condition": "derive Euler-Lagrange/RG stationary dynamics that also carries transferred roles",
        },
        {
            "rank": 6,
            "closure_task": "blind modern-physics frozen-kernel confirmation",
            "current_status": "test matrix exists but no external blind confirmation is claimed",
            "computable_exit_condition": "pre-register frozen-kernel CMB/LSS or GW/PTA no-retune predictions against exponential/spline baselines",
        },
    ]

    return {
        "where_the_theory_currently_looks_most_toe_like": [
            "one stable strict kernel organizing multiple sectors",
            "sinusoidal positional/resonance encoding plus heavy-tailed attention bias, which is structurally close to a graph neural/energy-based information architecture",
            "reproducibility and stability evidence for the strict damping exponent and selected kernel",
            "a constructive metric-pushforward candidate now exists for the node/gauge mismatch, though it is not yet sourced",
        ],
        "neural_universe_interpretation": (
            "If the universe is modeled as a self-learning information network, the strict kernel looks like a frozen learned layer: "
            "cos(omega*d+phi) supplies positional resonance, 1/(1+beta*d^eta) supplies heavy-tailed attention memory, and variational/RG stationarity plays the role of energy-based learning. "
            "The observed symptom would not be metaphor alone if the same frozen kernel predicts independent CMB/LSS, GW/PTA, and cross-sector spectra without retuning."
        ),
        "strict_kernel_full_kernel_decision": {
            "toe_symptoms_present": toe_symptoms,
            "strict_stability_positive": stability_positive,
            "toe_likeness_score_0_to_1_not_probability": score,
            "open_blocker_count_from_p2636": len([row for row in blockers if row.get("is_closed") is False]),
            "is_full_kernel_now": False,
            "classification": "TOE_LIKE_STABLE_WORKING_KERNEL_WITH_A_CONSTRUCTIVE_NODE_REPAIR_CANDIDATE_NOT_FULL_KERNEL",
        },
        "professorial_closure_path": path,
        "next_honest_step": (
            "Promote neither the affine repair nor the strict kernel to finality yet.  The next proof-grade task is to derive or falsify the metric-pushforward source r(d)=4/3*(d-1) from nadsoliton topology/selector dynamics. "
            "If it is derived, it becomes the first stability-to-role interface theorem; if it fails, remove the legacy integer node/gauge role from the transfer list and continue with beta/source and inverse-hierarchy blockers."
        ),
    }


def acceptance(phase_cert: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "identity_same_distance_nodes_exact": all(row["is_exact_formula_zero"] for row in phase_cert["identity_phase_rows"]),
        "silent_reindexing_or_simple_map_suffices": phase_cert["silent_identity_or_simple_reindexing_passes"],
        "constructive_affine_pushforward_exists": phase_cert["constructive_affine_metric_pushforward_exists"],
        "metric_pushforward_source_theorem_present": False,
        "legacy_node_gauge_role_transfer_allowed_now": False,
    }
    return {
        "gates": gates,
        "phase_frequency_node_gauge_certificate_exported": all(gates.values()),
        "verdict": (
            "P2637 improves the situation from a bare mismatch to a constructive exact repair candidate: r(d)=4/3*(d-1) maps the legacy integer nodes to the strict cosine zero lattice. "
            "But because that repair is a nontrivial metric pushforward, not identity retention, it does not by itself export the node/gauge role-transfer certificate or ToE closure."
        ),
    }


def write_markdown(payload: dict[str, Any]) -> None:
    phase = payload["phase_node_quotient_exhaustion_certificate"]
    acc = payload["source_acceptance"]
    closure = payload["professorial_toe_closure_path"]
    lines = [
        "# P2637/S1587 phase-node quotient exhaustion and ToE closure path audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps phase/node zero-lattice content, quotient/reindexing/metric-warp content, neural positional-attention content, ToE blocker content, and modern-physics closure tests before adding the finite certificate.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Exact phase/node quotient exhaustion",
        "",
        f"Strict phase parameters: `omega={phase['strict_phase_parameters_pi_units']['omega']}*pi`, `phi={phase['strict_phase_parameters_pi_units']['phi']}*pi`.",
        f"Exact zero lattice: `{phase['exact_zero_lattice_formula']}`.",
        f"Legacy declared nodes: `{phase['declared_legacy_nodes']}`.",
        f"Strict zero lattice sample: `{phase['audited_zero_lattice']}`.",
        f"Silent identity/simple reindexing passes? `{phase['silent_identity_or_simple_reindexing_passes']}`.",
        f"Constructive affine metric pushforward exists? `{phase['constructive_affine_metric_pushforward_exists']}`.",
        "",
        "| map class | definition | passes exact node certificate? | interpretation |",
        "| --- | --- | --- | --- |",
    ])
    for row in phase["exhausted_map_classes"]:
        lines.append(f"| {row['map_class']} | `{row['definition']}` | {row['passes_exact_node_certificate']} | {row['interpretation']} |")
    lines.extend([
        "",
        "## Acceptance",
        "",
        acc["verdict"],
        "",
        "## Professorial ToE closure path",
        "",
        "Where the theory currently looks most ToE-like:",
    ])
    for item in closure["where_the_theory_currently_looks_most_toe_like"]:
        lines.append(f"- {item}")
    lines.extend([
        "",
        closure["neural_universe_interpretation"],
        "",
        f"Full-kernel now? `{closure['strict_kernel_full_kernel_decision']['is_full_kernel_now']}`.",
        f"Classification: `{closure['strict_kernel_full_kernel_decision']['classification']}`.",
        "",
        "| rank | closure task | current status | computable exit condition |",
        "| ---: | --- | --- | --- |",
    ])
    for row in closure["professorial_closure_path"]:
        lines.append(f"| {row['rank']} | {row['closure_task']} | {row['current_status']} | {row['computable_exit_condition']} |")
    lines.extend([
        "",
        "## Recommended next honest step",
        "",
        closure["next_honest_step"],
        "",
        "No ToE closure, full-kernel finality, bridge completion, role-transfer, selector-source discharge, or role-bearing `L_total` is claimed.",
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    src = {name: load_json(path) for name, path in SOURCE_FILES.items() if path.suffix == ".json"}
    phase_cert = phase_node_quotient_exhaustion()
    closure_path = professorial_toe_closure_path(phase_cert, src)
    payload: dict[str, Any] = {
        "status": "P2637_PHASE_NODE_QUOTIENT_EXHAUSTION_TOE_CLOSURE_PATH_AUDIT_NO_SILENT_ROLE_TRANSFER",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_fingerprints": source_fingerprints(),
        "phase_node_quotient_exhaustion_certificate": phase_cert,
        "source_acceptance": acceptance(phase_cert),
        "professorial_toe_closure_path": closure_path,
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)

    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2637/S1587 phase-node metric-pushforward guard",
        "\n## P2637/S1587 phase-node metric-pushforward guard\n\n"
        "`P2637/S1587` exhausts simple phase/node repairs for the legacy integer node/gauge list versus the strict cosine zero lattice.  Identity, pure translation, and pure scale/reindexing do not certify the legacy nodes.  A constructive affine metric pushforward `r(d)=4/3*(d-1)` maps `2,5,8,11` to the strict zero lattice exactly, but this is a non-silent distance-coordinate theorem obligation, not a completed role-transfer certificate.  The strict kernel remains ToE-like and stable, but not a full kernel until this metric-pushforward source and the remaining beta, inverse-hierarchy, selector, `L_total`, and blind empirical blockers are discharged.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2637/S1587 phase-node Ltotal closure-path guard",
        "\n## P2637/S1587 phase-node Ltotal closure-path guard\n\n"
        "`P2637/S1587` supplies a constructive node/gauge repair candidate through a metric pushforward, but it does not re-enable role-bearing `L_total`: the pushforward must be derived from nadsoliton topology/selector dynamics before the legacy node/gauge role can transfer to strict dynamics.  Neural-universe language is therefore a test generator, not a closure theorem.\n",
    )
    return payload


if __name__ == "__main__":
    main()
