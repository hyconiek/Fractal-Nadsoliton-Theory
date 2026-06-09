#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from statistics import mean
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2632_s1582_neural_legacy_strict_characteristic_retention_audit.json"
MD = GEN / "p2632_s1582_neural_legacy_strict_characteristic_retention_audit.md"

REPO_ROOT = REPO
QW1729_JSON = REPO_ROOT / "report_qw1729_nadsoliton_kernel_characteristics_map.json"
P2625_JSON = GEN / "p2625_s1575_nonlinear_damping_completion_source_classification.json"
P2629_JSON = GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.json"
P2630_JSON = GEN / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.json"
P2631_JSON = GEN / "p2631_s1581_neural_information_flux_beta_criticality_audit.json"

ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA = math.pi / 4.0
PHI = math.pi / 6.0
BETA_TORS = 0.01
BETA_STRICT = 1.0
ETA_STRICT = 9.0 / 5.0
SAMPLE_D = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0]

NEGATIVE_EXPORT_FLAGS = [
    "all_legacy_characteristics_preserved_claimed",
    "strict_kernel_finality_claimed",
    "positive_beta_renormalization_source_exported",
    "nonlinear_damping_completion_source_exported",
    "role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "qw2191_discharged",
    "toe_closure_claimed",
]


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
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.lean", "-g", "*.json",
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
    # Content-first anti-duplication: characteristic-preservation content, not ticket names/numbers.
    patterns = {
        "nadsoliton_legacy_characteristics_content": (
            "4-bit informational capacity|octave resonance|hexagonal lattice|inter-layer torsion|"
            "topological tunneling|combined effective form|nadsoliton.*kernel.*characteristics"
        ),
        "legacy_strict_bridge_split_content": (
            "K_legacy_ont|K_strict_gate|legacy kernel|strict kernel|completion bridge|"
            "strict-side characteristics|intermediate bridge kernel|completed.*strict"
        ),
        "neural_attention_prism_content": (
            "positional encoding|heavy-tailed attention|attention bias|Energy-Based|Boltzmann|"
            "edge of chaos|information throughput|information flux"
        ),
        "damping_completion_obstruction_content": (
            "beta_tors|positive_beta_renormalization_source|Z_beta|normalization-gauge|"
            r"beta_micro/beta_strict|nonlinear damping|d\^eta|d\^\(9/5\)|9/5"
        ),
        "closure_guard_content": (
            "role-transfer|role-bearing L_total|QW-2191|ToE closure|selector source|"
            "not.*full bridge|no.*closure"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for neural legacy-to-strict characteristic retention", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def legacy_denominator(d: float) -> float:
    return 1.0 + BETA_TORS * d


def strict_denominator(d: float) -> float:
    return 1.0 + BETA_STRICT * (d ** ETA_STRICT)


def phase(d: float) -> float:
    return math.cos(OMEGA * d + PHI)


def legacy_kernel(d: float) -> float:
    return ALPHA_GEO * phase(d) / legacy_denominator(d)


def strict_kernel(d: float) -> float:
    return phase(d) / strict_denominator(d)


def effective_log_slope_denominator(beta: float, power: float, d: float) -> float:
    # d/d(log d) log(1 + beta*d^power)
    return (power * beta * (d ** power)) / (1.0 + beta * (d ** power))


def finite_sample_geometry() -> dict[str, Any]:
    rows = []
    for d in SAMPLE_D:
        legacy_denom = legacy_denominator(d)
        strict_denom = strict_denominator(d)
        legacy_attention = 1.0 / legacy_denom
        strict_attention = 1.0 / strict_denom
        rows.append({
            "d": d,
            "phase_legacy": phase(d),
            "phase_strict": phase(d),
            "phase_difference": 0.0,
            "legacy_denominator": legacy_denom,
            "strict_denominator": strict_denom,
            "legacy_attention_bias": legacy_attention,
            "strict_attention_bias": strict_attention,
            "strict_attention_over_legacy_attention": strict_attention / legacy_attention,
            "legacy_effective_log_slope": effective_log_slope_denominator(BETA_TORS, 1.0, d),
            "strict_effective_log_slope": effective_log_slope_denominator(BETA_STRICT, ETA_STRICT, d),
            "legacy_kernel_raw": legacy_kernel(d),
            "strict_kernel_raw": strict_kernel(d),
            "strict_over_amplitude_normalized_legacy": strict_kernel(d) / (legacy_kernel(d) / ALPHA_GEO) if phase(d) != 0 else None,
        })
    return {
        "sample_domain": SAMPLE_D,
        "rows": rows,
        "phase_linf_difference": max(abs(row["phase_difference"]) for row in rows),
        "mean_strict_attention_over_legacy_attention": mean(row["strict_attention_over_legacy_attention"] for row in rows),
        "max_strict_denominator_minus_legacy_denominator": max(row["strict_denominator"] - row["legacy_denominator"] for row in rows),
        "strict_attention_at_d10_over_legacy": next(row["strict_attention_over_legacy_attention"] for row in rows if row["d"] == 10.0),
        "strict_slope_at_d10_minus_legacy_slope": next(row["strict_effective_log_slope"] - row["legacy_effective_log_slope"] for row in rows if row["d"] == 10.0),
        "interpretation": (
            "The resonance/positional phase is syntactically retained, but the denominator is not a small perturbation of the legacy torsion denominator.  "
            "On the audited window the strict attention becomes much more compressed at long range, especially near d=10."
        ),
    }


def characteristic_retention_matrix(qw1729: dict[str, Any], p2625: dict[str, Any], p2629: dict[str, Any], p2630: dict[str, Any], p2631: dict[str, Any]) -> list[dict[str, Any]]:
    p2625_text = json.dumps(p2625, sort_keys=True)
    p2629_accepts = bool(p2629.get("exact_source_gate", {}).get("accepts_positive_beta_renormalization_source", False))
    p2630_accepts = bool(p2630.get("bridge_source_truth_table", {}).get("current_accepts_bridge_positive_zbeta_source", False))
    p2631_accepts = bool(p2631.get("source_acceptance", {}).get("accepts_information_conservation_beta_identity", False))
    return [
        {
            "legacy_characteristic": "4-bit informational capacity / info-geometry amplitude alpha_geo",
            "neural_readout": "global gain / feature normalization of the coupling layer",
            "strict_retention_verdict": "modified_not_silent_preservation",
            "evidence": "strict kernel uses unit numerator amplitude; alpha_geo belongs to legacy/canonical role layer unless an explicit normalization/role-transfer theorem is supplied",
            "computable_status": "amplitude-normalized comparisons are allowed; raw amplitude role is not transferred",
            "preserved_as_final": False,
        },
        {
            "legacy_characteristic": "octave resonance frequency omega=pi/4",
            "neural_readout": "sinusoidal positional/resonance encoding frequency",
            "strict_retention_verdict": "formally_retained_but_source_still_guarded",
            "evidence": "legacy and strict phases share omega in the working formulas; first-principles selector/source closure remains guarded",
            "computable_status": "phase difference is zero on finite samples by formula identity",
            "preserved_as_final": False,
        },
        {
            "legacy_characteristic": "hexagonal/vacuum phase offset phi=pi/6",
            "neural_readout": "positional encoding phase offset / symmetry bias",
            "strict_retention_verdict": "formally_retained_but_source_still_guarded",
            "evidence": "legacy and strict phases share phi in the working formulas; QW-2191-style selector/source closure is not discharged",
            "computable_status": "phase difference is zero on finite samples by formula identity",
            "preserved_as_final": False,
        },
        {
            "legacy_characteristic": "inter-layer torsion damping beta_tors=0.01",
            "neural_readout": "legacy temperature/length-scale prior for attention damping",
            "strict_retention_verdict": "not_preserved_as_bare_parameter",
            "evidence": "strict damping uses beta=1 and eta=9/5; P2629/P2630 leave positive beta renormalization and legacy/UV normalization source negative",
            "computable_status": f"positive_beta_source_accepts={p2629_accepts}; bridge_source_accepts={p2630_accepts}",
            "preserved_as_final": False,
        },
        {
            "legacy_characteristic": "hyperbolic denominator 1/(1+beta_tors*d)",
            "neural_readout": "attention bias / memory tail",
            "strict_retention_verdict": "enriched_strict_successor_not_exact_retention",
            "evidence": "strict denominator is 1/(1+beta*d^eta); P2625 identifies exponent pushforward as partial but coefficient source remains missing",
            "computable_status": "conditional_completion_schema_present=" + str("fractal_pushforward_plus_independent_Z_beta" in p2625_text),
            "preserved_as_final": False,
        },
        {
            "legacy_characteristic": "joint effective nadsoliton coupling form",
            "neural_readout": "energy-based fractal attention kernel with oscillatory positional channel and heavy-tailed damping channel",
            "strict_retention_verdict": "working_completion_candidate_incomplete_as_toe_kernel",
            "evidence": "strict kernel preserves useful structure and adds nonlinear compression, but P2631 does not prove beta=1 as a neural critical point and P2629/P2630 keep bridge source open",
            "computable_status": f"p2631_beta_identity_accepts={p2631_accepts}",
            "preserved_as_final": False,
        },
    ]


def professorial_readout(matrix: list[dict[str, Any]], geometry: dict[str, Any], qw1729: dict[str, Any]) -> dict[str, Any]:
    retained_formally = [row for row in matrix if "formally_retained" in row["strict_retention_verdict"]]
    rejected_or_modified = [row for row in matrix if row["strict_retention_verdict"] not in {"formally_retained_but_source_still_guarded"}]
    closure_score = float(qw1729.get("kernel_closure_score", 0.0))
    return {
        "professorial_summary": (
            "A theoretical physicist who also understands modern neural networks would likely read the strict kernel as a sophisticated working successor of the legacy coupling kernel, not as a finished ToE kernel.  "
            "The cosine channel looks like a retained sinusoidal positional/resonance encoding; the denominator looks like a heavy-tailed attention bias on a fractal measure.  "
            "However, preserving the language of the architecture is weaker than proving that all nadsoliton characteristics have survived as theorem-level physical roles."
        ),
        "does_strict_preserve_all_legacy_characteristics": False,
        "why_not": [
            "alpha_geo is not silently preserved as the strict numerator amplitude or physical role",
            "omega and phi are formula-retained but their selector/source status remains guarded",
            "beta_tors=0.01 is not the strict beta=1 without the missing positive renormalization / UV normalization source",
            "linear legacy damping is replaced by nonlinear strict compression, which is a strict-side addition rather than exact retention",
            "neural edge-of-chaos/information-flux language did not produce a noncircular beta=1 identity in P2631",
        ],
        "formal_retention_count": len(retained_formally),
        "modified_or_not_final_count": len(rejected_or_modified),
        "legacy_qw1729_kernel_closure_score": closure_score,
        "long_range_compression_witness": {
            "strict_attention_at_d10_over_legacy": geometry["strict_attention_at_d10_over_legacy"],
            "strict_slope_at_d10_minus_legacy_slope": geometry["strict_slope_at_d10_minus_legacy_slope"],
        },
        "status": "STRICT_KERNEL_IS_ROBUST_WORKING_SUCCESSOR_BUT_INCOMPLETE_FOR_TOE_FINALITY",
    }


def source_acceptance(readout: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "all_legacy_characteristics_preserved_as_theorems": False,
        "strict_side_additions_have_completion_map_sources": False,
        "beta_tors_to_beta_1_bridge_source_closed": False,
        "phase_topological_selector_source_closed": False,
        "neural_criticality_beta_identity_closed": False,
        "no_role_transfer_or_toe_overclaim": True,
    }
    return {
        "gates": gates,
        "accepts_strict_kernel_as_complete_toe_kernel": all(gates.values()),
        "failed_gates": [name for name, value in gates.items() if not value],
        "status": "NO_FINALITY_PROMOTION_PROFESSORIAL_RETENTION_AUDIT",
        "reason": (
            "The strict kernel is better viewed as a mathematically sharper, neural-attention-like working completion candidate built from the legacy kernel, "
            "not as a final kernel that has already preserved every intended nadsoliton characteristic with theorem-level role transfer."
        ),
    }


def recommendation() -> str:
    return (
        "Next honest step: build a characteristic-by-characteristic completion certificate instead of a global identity claim.  "
        "For each legacy nadsoliton role, prove one typed statement: amplitude normalization/role map, phase-frequency selector source, "
        "damping exponent pushforward, positive beta/UV normalization source, and neural/information stability functional.  "
        "The most urgent computable substep is a dimensionless phase-stability or information-throughput functional whose extremum is derived before inserting beta=1; if that fails, formalize the strict kernel as an approximate finite-domain effective successor with explicit residuals."
    )


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2632/S1582 neural legacy-to-strict characteristic retention audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps research content about nadsoliton characteristics, bridge completion, neural attention language, damping obstructions, and closure guards before adding the audit.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Professorial readout",
        "",
        payload["professorial_readout"]["professorial_summary"],
        "",
        f"Verdict: `{payload['professorial_readout']['status']}`.",
        "",
        "The strict kernel does **not** yet preserve all legacy/nadsoliton characteristics as final theorem-level roles.",
        "",
        "## Finite computational witness",
        "",
        payload["finite_sample_geometry"]["interpretation"],
        "",
        "| d | strict attention / legacy attention | legacy slope | strict slope |",
        "| ---: | ---: | ---: | ---: |",
    ])
    for row in payload["finite_sample_geometry"]["rows"]:
        lines.append(
            f"| {row['d']:.2f} | {row['strict_attention_over_legacy_attention']:.10f} | "
            f"{row['legacy_effective_log_slope']:.10f} | {row['strict_effective_log_slope']:.10f} |"
        )
    lines.extend([
        "",
        "## Characteristic retention matrix",
        "",
        "| legacy characteristic | neural readout | retention verdict |",
        "| --- | --- | --- |",
    ])
    for row in payload["characteristic_retention_matrix"]:
        lines.append(f"| {row['legacy_characteristic']} | {row['neural_readout']} | `{row['strict_retention_verdict']}` |")
    lines.extend([
        "",
        "## Verdict",
        "",
        payload["source_acceptance"]["reason"],
        "",
        "No role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure is reopened.",
        "",
        "## Recommended next honest step",
        "",
        payload["recommended_next_honest_step"],
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    qw1729 = load_json(QW1729_JSON)
    p2625 = load_json(P2625_JSON)
    p2629 = load_json(P2629_JSON)
    p2630 = load_json(P2630_JSON)
    p2631 = load_json(P2631_JSON)
    geometry = finite_sample_geometry()
    matrix = characteristic_retention_matrix(qw1729, p2625, p2629, p2630, p2631)
    readout = professorial_readout(matrix, geometry, qw1729)
    acceptance = source_acceptance(readout)
    payload: dict[str, Any] = {
        "status": "P2632_NEURAL_LEGACY_STRICT_CHARACTERISTIC_RETENTION_AUDIT_NO_FINALITY_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "sources": {
            "QW1729_JSON": rel(QW1729_JSON),
            "P2625_JSON": rel(P2625_JSON),
            "P2629_JSON": rel(P2629_JSON),
            "P2630_JSON": rel(P2630_JSON),
            "P2631_JSON": rel(P2631_JSON),
        },
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "finite_sample_geometry": geometry,
        "characteristic_retention_matrix": matrix,
        "professorial_readout": readout,
        "source_acceptance": acceptance,
        "recommended_next_honest_step": recommendation(),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2632/S1582 neural legacy-to-strict characteristic retention audit guard",
        "\n## P2632/S1582 neural legacy-to-strict characteristic retention audit guard\n\n"
        "`P2632/S1582` gives the professorial neural/physics readout of the legacy-to-strict passage.  The strict kernel is a strong working successor: the cosine channel retains the sinusoidal positional/resonance encoding, and the strict denominator supplies a heavy-tailed fractal attention/compression channel.  But the finite audit shows that this is not exact preservation of all legacy nadsoliton characteristics: `alpha_geo` is not silently retained as strict amplitude, `omega/phi` remain selector/source-guarded, `beta_tors=0.01` is not `beta=1` without the missing UV/positive renormalization source, and the nonlinear denominator is a strict-side addition rather than the old linear torsion denominator.  Therefore the strict kernel remains an incomplete/working ToE kernel candidate; no role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure is reopened.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2632/S1582 neural legacy-strict retention Ltotal guard",
        "\n## P2632/S1582 neural legacy-strict retention Ltotal guard\n\n"
        "`P2632/S1582` keeps `L_total` role claims closed: the neural-attention reading makes the strict kernel intelligible as an enriched successor of the legacy kernel, but it does not by itself prove characteristic-by-characteristic role transfer or final ToE kernel status.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "strict_preserves_all": result["professorial_readout"]["does_strict_preserve_all_legacy_characteristics"],
        "d10_attention_ratio": result["finite_sample_geometry"]["strict_attention_at_d10_over_legacy"],
        "accepts_final_kernel": result["source_acceptance"]["accepts_strict_kernel_as_complete_toe_kernel"],
        "recommended_next": result["recommended_next_honest_step"],
        "out": rel(OUT),
        "md": rel(MD),
    }, indent=2, sort_keys=True, ensure_ascii=False))
