#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.json"
MD = GEN / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.md"

REPO_ROOT = REPO
EXTERNAL_DIAGRAM = Path("/media/magazyn/home/teresa/Dokumenty/HYC/TOE/edison/DIAGRAMS_KERNEL_TRANSFORMATION.md")
LOCAL_DIAGRAM = REPO_ROOT / "DIAGRAMS_KERNEL_TRANSFORMATION.md"
P2625_JSON = GEN / "p2625_s1575_nonlinear_damping_completion_source_classification.json"
P2631_JSON = GEN / "p2631_s1581_neural_information_flux_beta_criticality_audit.json"
P2632_JSON = GEN / "p2632_s1582_neural_legacy_strict_characteristic_retention_audit.json"

ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA = math.pi / 4.0
PHI = math.pi / 6.0
BETA_TORS = 0.01
STRICT_BETA = 1.0
STRICT_ETA = 9.0 / 5.0
DISTANCE_GRID = [1.0, 2.0, 3.0, 5.0, 7.0, 8.0, 10.0, 11.0, 12.0]
DECLARED_INTEGER_NODES = [2.0, 5.0, 8.0, 11.0]

NEGATIVE_EXPORT_FLAGS = [
    "all_diagram_characteristics_preserved_claimed",
    "legacy_inverse_hierarchy_preserved_by_strict_claimed",
    "diagram_integer_node_pattern_formula_exact_claimed",
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


def read_diagram() -> tuple[Path, str, bool]:
    if EXTERNAL_DIAGRAM.exists():
        return EXTERNAL_DIAGRAM, EXTERNAL_DIAGRAM.read_text(encoding="utf-8"), True
    return LOCAL_DIAGRAM, LOCAL_DIAGRAM.read_text(encoding="utf-8"), False


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
    # Content-first anti-duplication: diagram mechanisms and strict residual characteristics, not ticket numbers.
    patterns = {
        "legacy_four_mechanism_content": (
            "K_geo|K_res|K_tors|K_topo|viscosity|lepko|resonance|torsion|"
            "topological paths|topological tunneling|winding"
        ),
        "legacy_effective_kernel_content": (
            r"alpha_geo|beta_tors|cos\(omega|cos\(ω|1/\(1\+β|1/\(1\+beta|inverse hierarchy|hyperbolic"
        ),
        "diagram_node_and_wilson_content": (
            "Wilson loops|13\\.6|d=7|nodes at d=2|period: 3|SU\\(3\\)|SU\\(2\\)|U\\(1\\)"
        ),
        "strict_fractal_compression_content": (
            r"fractal compression|nonlinear compression|strict damping|d\^eta|d\^\(9/5\)|D_f=eta|heavy-tailed attention|information flux"
        ),
        "bridge_guard_content": (
            "positive_beta_renormalization_source|normalization-gauge|role-bearing L_total|role-transfer|QW-2191|ToE closure|selector source"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for diagram-grounded legacy characteristic preservation", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def phase(d: float) -> float:
    return math.cos(OMEGA * d + PHI)


def legacy_attention(d: float) -> float:
    return 1.0 / (1.0 + BETA_TORS * d)


def strict_attention(d: float) -> float:
    return 1.0 / (1.0 + STRICT_BETA * d ** STRICT_ETA)


def legacy_kernel_amplitude_normalized(d: float) -> float:
    return phase(d) * legacy_attention(d)


def strict_kernel(d: float) -> float:
    return phase(d) * strict_attention(d)


def extract_diagram_claims(text: str) -> dict[str, Any]:
    claims = {
        "has_four_mechanism_product": "K_total = K_geo" in text and "K_res" in text and "K_tors" in text and "K_topo" in text,
        "has_effective_legacy_kernel": "K(d) = α" in text or "K(d) = alpha" in text or "K(d) = α_geo" in text,
        "claims_95_percent_accuracy": "95% accuracy" in text,
        "claims_100x_faster": "100× faster" in text or "100x faster" in text,
        "claims_wilson_amplification_13_6": "13.6×" in text or "13.6x" in text,
        "claims_inverse_hierarchy": "inverse hierarchy" in text.lower(),
        "claims_hyperbolic_not_exponential": "Hyperbolic, not exponential" in text or "hyperbolic damping" in text.lower(),
        "claims_integer_node_pattern": "nodes at d=2,5,8,11" in text or "nodes at d=2" in text,
        "claims_fractal_path_count": "Path count" in text and "d^1.6" in text,
        "claims_four_minimal_parameters_absorb_physics": "4 minimal parameters absorb" in text,
    }
    snippets = {}
    for key in ["K_total", "K(d) =", "13.6", "nodes at d=2", "Hyperbolic"]:
        match = re.search(re.escape(key), text)
        if match:
            start = max(match.start() - 220, 0)
            end = min(match.end() + 360, len(text))
            snippets[key] = text[start:end]
    return {"claims": claims, "claim_count": sum(bool(v) for v in claims.values()), "snippets": snippets}


def true_phase_zeros(n: int = 5) -> list[float]:
    return [((math.pi / 2.0 + k * math.pi) - PHI) / OMEGA for k in range(n)]


def finite_witness() -> dict[str, Any]:
    rows = []
    for d in DISTANCE_GRID:
        legacy = legacy_kernel_amplitude_normalized(d)
        strict = strict_kernel(d)
        rows.append({
            "d": d,
            "cos_phase": phase(d),
            "legacy_attention_linear_beta_tors": legacy_attention(d),
            "strict_attention_power_beta_1_eta_9_5": strict_attention(d),
            "strict_attention_over_legacy_attention": strict_attention(d) / legacy_attention(d),
            "legacy_amplitude_normalized_abs_kernel": abs(legacy),
            "strict_abs_kernel": abs(strict),
            "strict_abs_over_legacy_abs": abs(strict) / abs(legacy) if abs(legacy) > 0 else None,
        })
    d1 = next(row for row in rows if row["d"] == 1.0)
    d7 = next(row for row in rows if row["d"] == 7.0)
    node_rows = []
    for d in DECLARED_INTEGER_NODES:
        node_rows.append({
            "declared_node_d": d,
            "cos_phase_value": phase(d),
            "abs_cos_residual": abs(phase(d)),
            "is_formula_zero_with_tolerance_1e_12": abs(phase(d)) < 1.0e-12,
        })
    zeros = true_phase_zeros(5)
    return {
        "distance_grid": DISTANCE_GRID,
        "rows": rows,
        "inverse_hierarchy_ratio_abs_k7_over_abs_k1": {
            "legacy_amplitude_normalized": d7["legacy_amplitude_normalized_abs_kernel"] / d1["legacy_amplitude_normalized_abs_kernel"],
            "strict": d7["strict_abs_kernel"] / d1["strict_abs_kernel"],
            "strict_preserves_legacy_ratio_above_one": (d7["strict_abs_kernel"] / d1["strict_abs_kernel"]) > 1.0,
        },
        "declared_integer_node_audit": {
            "declared_nodes": DECLARED_INTEGER_NODES,
            "rows": node_rows,
            "all_declared_integer_nodes_are_formula_zeros": all(row["is_formula_zero_with_tolerance_1e_12"] for row in node_rows),
            "true_continuous_zero_lattice_first_values": zeros,
            "note": "With omega=pi/4 and phi=pi/6, exact continuous zeros are d=4/3+4n, not the declared integer pattern 2,5,8,11. Strict has the same phase, so it shares this formula-level issue rather than repairing it.",
        },
        "d10_long_range_compression": {
            "strict_attention_over_legacy_attention": next(row["strict_attention_over_legacy_attention"] for row in rows if row["d"] == 10.0),
            "strict_abs_over_legacy_abs": next(row["strict_abs_over_legacy_abs"] for row in rows if row["d"] == 10.0),
        },
    }


def preservation_matrix(p2625: dict[str, Any], p2631: dict[str, Any], witness: dict[str, Any]) -> list[dict[str, Any]]:
    p2625_text = json.dumps(p2625, sort_keys=True)
    p2631_accepts = bool(p2631.get("source_acceptance", {}).get("accepts_information_conservation_beta_identity", False))
    inverse_preserved = witness["inverse_hierarchy_ratio_abs_k7_over_abs_k1"]["strict_preserves_legacy_ratio_above_one"]
    integer_nodes_exact = witness["declared_integer_node_audit"]["all_declared_integer_nodes_are_formula_zeros"]
    return [
        {
            "diagram_characteristic": "four-mechanism product K_geo*K_res*(1+0.2K_tors)*K_topo",
            "strict_reading": "strict kernel compresses this to a two-channel architecture: phase/resonance numerator plus fractal attention denominator",
            "preservation_verdict": "compressed_architectural_successor_not_full_mechanism_retention",
            "computable_evidence": "strict formula no longer exposes four independent mechanism factors",
            "preserved_as_final_toe_role": False,
        },
        {
            "diagram_characteristic": "cos(omega*d+phi) from torsion oscillation plus resonance synchronization",
            "strict_reading": "same sinusoidal positional/resonance encoding is retained",
            "preservation_verdict": "formula_retained_source_guarded",
            "computable_evidence": "phase channel identical on all sampled d, but selector/source theorem is still separate",
            "preserved_as_final_toe_role": False,
        },
        {
            "diagram_characteristic": "integer node pattern d=2,5,8,11 used for gauge/generation interpretation",
            "strict_reading": "strict inherits the same phase formula, but the declared integer nodes are not exact zeros of that formula",
            "preservation_verdict": "diagram_claim_not_formula_exact",
            "computable_evidence": f"all_declared_integer_nodes_are_formula_zeros={integer_nodes_exact}",
            "preserved_as_final_toe_role": False,
        },
        {
            "diagram_characteristic": "inverse hierarchy / distant octaves remain strongly coupled",
            "strict_reading": "strict beta=1, eta=9/5 produces much stronger long-range compression on the same d grid",
            "preservation_verdict": "not_preserved_numerically_without_reinterpretation",
            "computable_evidence": f"strict_abs_K7_over_K1={witness['inverse_hierarchy_ratio_abs_k7_over_abs_k1']['strict']:.6f}; legacy_abs_K7_over_K1={witness['inverse_hierarchy_ratio_abs_k7_over_abs_k1']['legacy_amplitude_normalized']:.6f}",
            "preserved_as_final_toe_role": False,
        },
        {
            "diagram_characteristic": "hyperbolic not exponential damping from fractal topological path summation",
            "strict_reading": "strict retains a heavy-tailed non-exponential denominator but changes linear beta_tors*d to beta*d^(9/5)",
            "preservation_verdict": "modified_fractal_compression_successor",
            "computable_evidence": "P2625 conditional schema present=" + str("fractal_pushforward_plus_independent_Z_beta" in p2625_text),
            "preserved_as_final_toe_role": False,
        },
        {
            "diagram_characteristic": "beta_tors controls optimal inverse-hierarchy range 0.01..0.08",
            "strict_reading": "strict uses beta=1; P2631 did not prove beta=1 as a target-independent neural critical point",
            "preservation_verdict": "not_preserved_as_legacy_temperature_parameter",
            "computable_evidence": f"p2631_information_conservation_beta_identity_accepts={p2631_accepts}",
            "preserved_as_final_toe_role": False,
        },
        {
            "diagram_characteristic": "alpha_geo, omega, phi, beta absorb all physics in the effective legacy form",
            "strict_reading": "strict removes alpha_geo from the raw numerator and changes beta semantics; absorption claim needs role-transfer audit",
            "preservation_verdict": "parameter_absorption_not_silently_transferred",
            "computable_evidence": "strict raw numerator amplitude is 1 while legacy amplitude is alpha_geo=4*ln(2)",
            "preserved_as_final_toe_role": False,
        },
    ]


def acceptance(matrix: list[dict[str, Any]], witness: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "diagram_four_mechanisms_retained_as_typed_sources": False,
        "integer_node_pattern_formula_exact_or_repaired": witness["declared_integer_node_audit"]["all_declared_integer_nodes_are_formula_zeros"],
        "legacy_inverse_hierarchy_numerically_retained": witness["inverse_hierarchy_ratio_abs_k7_over_abs_k1"]["strict_preserves_legacy_ratio_above_one"],
        "beta_tors_to_strict_beta_source_closed": False,
        "alpha_geo_role_transfer_closed": False,
        "no_role_transfer_or_toe_overclaim": True,
    }
    return {
        "gates": gates,
        "accepts_all_diagram_characteristics_preserved_by_strict": all(gates.values()),
        "failed_gates": [name for name, value in gates.items() if not value],
        "status": "DIAGRAM_GROUNDED_PRESERVATION_REJECTS_STRICT_FINALITY",
        "reason": (
            "The diagram-grounded audit strengthens the previous verdict: strict is an enriched successor, but it does not preserve every legacy diagram characteristic as a theorem-level ToE role.  "
            "In particular, the same phase formula does not exactly realize the diagram's integer node list, and strict long-range damping suppresses the d=7/d=1 inverse-hierarchy ratio on the audited amplitude-normalized grid."
        ),
    }


def recommendation() -> str:
    return (
        "Next honest step: split the bridge into two falsifiable certificates before any role audit: "
        "(1) a phase-frequency certificate reconciling the diagram's integer node/gauge interpretation with the exact cos(omega*d+phi) zero lattice, "
        "and (2) an inverse-hierarchy preservation or reinterpretation theorem showing which distance measure/domain lets strict heavy-tailed compression preserve the legacy Wilson-loop/distant-octave role.  "
        "Only after those pass should alpha_geo/beta_tors role transfer be retried; otherwise mark strict as an approximate finite-domain successor with explicit residuals."
    )


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2633/S1583 diagram-grounded strict kernel characteristic preservation audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps diagram mechanisms, node/inverse-hierarchy content, strict fractal compression content, and closure guards before adding the audit.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    src = payload["diagram_source"]
    lines.extend([
        "",
        "## Diagram source",
        "",
        f"Read diagram file: `{src['path']}`; external path available: `{src['external_path_available']}`.",
        "",
        "## Computational witnesses",
        "",
        "### Inverse hierarchy check",
        "",
        f"Legacy amplitude-normalized `|K(7)|/|K(1)|`: `{payload['finite_witness']['inverse_hierarchy_ratio_abs_k7_over_abs_k1']['legacy_amplitude_normalized']:.10f}`.",
        f"Strict `|K(7)|/|K(1)|`: `{payload['finite_witness']['inverse_hierarchy_ratio_abs_k7_over_abs_k1']['strict']:.10f}`.",
        "",
        "### Declared integer node check",
        "",
        payload["finite_witness"]["declared_integer_node_audit"]["note"],
        "",
        "| declared node d | cos residual | formula zero? |",
        "| ---: | ---: | --- |",
    ])
    for row in payload["finite_witness"]["declared_integer_node_audit"]["rows"]:
        lines.append(f"| {row['declared_node_d']:.0f} | {row['abs_cos_residual']:.10f} | {row['is_formula_zero_with_tolerance_1e_12']} |")
    lines.extend([
        "",
        "## Preservation matrix",
        "",
        "| diagram characteristic | strict reading | verdict |",
        "| --- | --- | --- |",
    ])
    for row in payload["preservation_matrix"]:
        lines.append(f"| {row['diagram_characteristic']} | {row['strict_reading']} | `{row['preservation_verdict']}` |")
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
    diagram_path, diagram_text, external_available = read_diagram()
    p2625 = load_json(P2625_JSON)
    p2631 = load_json(P2631_JSON)
    p2632 = load_json(P2632_JSON)
    diagram_claims = extract_diagram_claims(diagram_text)
    witness = finite_witness()
    matrix = preservation_matrix(p2625, p2631, witness)
    source_acceptance = acceptance(matrix, witness)
    payload: dict[str, Any] = {
        "status": "P2633_DIAGRAM_GROUNDED_STRICT_KERNEL_CHARACTERISTIC_PRESERVATION_AUDIT_NO_FINALITY_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "diagram_source": {
            "requested_external_path": str(EXTERNAL_DIAGRAM),
            "external_path_available": external_available,
            "path": str(diagram_path) if diagram_path.is_absolute() else rel(diagram_path),
            "sha256": hashlib.sha256(diagram_text.encode("utf-8")).hexdigest(),
        },
        "sources": {
            "P2625_JSON": rel(P2625_JSON),
            "P2631_JSON": rel(P2631_JSON),
            "P2632_JSON": rel(P2632_JSON),
            "p2632_status": p2632.get("status"),
        },
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "diagram_claims": diagram_claims,
        "finite_witness": witness,
        "preservation_matrix": matrix,
        "source_acceptance": source_acceptance,
        "recommended_next_honest_step": recommendation(),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2633/S1583 diagram-grounded strict characteristic preservation audit guard",
        "\n## P2633/S1583 diagram-grounded strict characteristic preservation audit guard\n\n"
        "`P2633/S1583` grounds the retention audit in `DIAGRAMS_KERNEL_TRANSFORMATION.md`: legacy `K_total` encodes four mechanisms, inverse hierarchy/Wilson-loop distant-octave persistence, hyperbolic-not-exponential damping, integer node/gauge interpretation, and four-parameter absorption.  The strict kernel retains the cosine phase channel and gives a nonlinear heavy-tailed fractal compression successor, but the finite audit rejects full preservation: the declared integer nodes `d=2,5,8,11` are not exact zeros of `cos(pi*d/4+pi/6)`, strict `|K(7)|/|K(1)|` on the amplitude-normalized grid is below one while legacy is above one, and `beta_tors/alpha_geo` roles remain untransferred.  Thus strict remains an enriched working successor rather than a final ToE kernel; no role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure is reopened.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2633/S1583 diagram-grounded retention Ltotal guard",
        "\n## P2633/S1583 diagram-grounded retention Ltotal guard\n\n"
        "`P2633/S1583` keeps `L_total` role claims closed after reading the legacy diagram characteristics: the strict kernel's neural/fractal-compression architecture is a working successor, not a proof that the diagram's inverse hierarchy, node/gauge, `alpha_geo`, or `beta_tors` roles have transferred to strict dynamics.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "diagram_path": result["diagram_source"]["path"],
        "external_path_available": result["diagram_source"]["external_path_available"],
        "legacy_k7_over_k1": result["finite_witness"]["inverse_hierarchy_ratio_abs_k7_over_abs_k1"]["legacy_amplitude_normalized"],
        "strict_k7_over_k1": result["finite_witness"]["inverse_hierarchy_ratio_abs_k7_over_abs_k1"]["strict"],
        "integer_nodes_exact": result["finite_witness"]["declared_integer_node_audit"]["all_declared_integer_nodes_are_formula_zeros"],
        "accepts_all_preserved": result["source_acceptance"]["accepts_all_diagram_characteristics_preserved_by_strict"],
        "recommended_next": result["recommended_next_honest_step"],
        "out": rel(OUT),
        "md": rel(MD),
    }, indent=2, sort_keys=True, ensure_ascii=False))
