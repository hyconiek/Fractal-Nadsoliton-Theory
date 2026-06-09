#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import math
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2625_s1575_nonlinear_damping_completion_source_classification.json"
MD = GEN / "p2625_s1575_nonlinear_damping_completion_source_classification.md"

SOURCE_FILES = {
    "K1_KERNEL_SPLIT_NOTE": ROOT / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md",
    "S2_STRATEGIC_PRIORITY": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "P2618_COMPLETION_OBSTRUCTION": GEN / "p2618_s1568_analytic_legacy_to_strict_completion_obstruction.json",
    "P2620_TWO_OBSTRUCTION_CUT": GEN / "p2620_s1570_p2618_p2619_bridge_two_obstruction_cut.json",
    "P2624_BLOCKER_RECOMMENDATION": GEN / "p2624_s1574_current_blockers_and_next_step_recommendation.json",
}

BETA_TORS = Fraction(1, 100)
BETA_STRICT = Fraction(1, 1)
ETA_STRICT = Fraction(9, 5)
Z_REQUIRED = BETA_STRICT / BETA_TORS

NEGATIVE_EXPORT_FLAGS = [
    "unconditional_nonlinear_damping_completion_source_exported",
    "scalar_beta_tors_to_beta_theorem_exported",
    "full_legacy_to_strict_bridge_revalidated",
    "orientation_odd_selector_source_exported",
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
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.lean",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:240]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: the patterns are research concepts, not packet identifiers.
    patterns = {
        "nonlinear_compression_law_content": (
            r"nonlinear damping|nonlinear compression|d\^eta|d\^eta|power-law damping|"
            "fractal projection|fractal measure|measure pushforward|scale-covariant|codimension"
        ),
        "legacy_linear_torsion_content": (
            r"beta_tors|torsion damping|linear torsion|1 \+ beta_tors|legacy kernel|"
            "K_legacy_ont|linear denominator"
        ),
        "strict_beta_eta_content": (
            "Z_beta|delta_eta|eta=1.8|eta=9/5|beta=1.0|strict damping|"
            "K_strict_gate|strict denominator"
        ),
        "hydrodynamic_source_content": (
            "hydrodynamic|compressible|compressibility|attenuation|renormalization source|"
            "micro-supported|source theorem|positive character"
        ),
        "guard_and_nonfit_content": (
            "completion map|bridge-source cut|no fit|no numerical fit|a priori|obstruction|"
            "role-bearing L_total|QW-2191|ToE closure"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic anti-duplication audit for nonlinear damping completion", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def monomial_identity_obstruction(left_power: Fraction, left_coeff: Fraction, right_power: Fraction, right_coeff: Fraction) -> dict[str, Any]:
    exact = left_power == right_power and left_coeff == right_coeff
    return {
        "left": {"coefficient": str(left_coeff), "power": str(left_power)},
        "right": {"coefficient": str(right_coeff), "power": str(right_power)},
        "identity_for_all_positive_d": exact,
        "proof": (
            "For positive coefficients, equality 1+a*d^p = 1+b*d^q for every d>0 implies "
            "a*d^p=b*d^q. Dividing by d^q gives a*d^(p-q)=b for every d, hence p=q; then a=b."
        ),
        "obstruction_if_not_exact": [] if exact else [
            "power_mismatch" if left_power != right_power else "coefficient_mismatch",
            "coefficient_mismatch" if left_coeff != right_coeff else "power_mismatch",
        ],
    }


def value(power: Fraction, coeff: Fraction, d: Fraction) -> Fraction:
    if power.denominator == 1:
        return Fraction(1, 1) + coeff * (d ** power.numerator)
    # The finite certificate uses float only for residual magnitudes; the theorem above is exact.
    return Fraction.from_float(1.0 + float(coeff) * (float(d) ** (power.numerator / power.denominator))).limit_denominator(10**12)


def candidate_models() -> list[dict[str, Any]]:
    sample_ds = [Fraction(1, 1), Fraction(2, 1), Fraction(5, 1), Fraction(10, 1)]
    models: list[dict[str, Any]] = []

    def add_model(name: str, power: Fraction, coeff: Fraction, required_atoms: list[str], verdict: str, note: str) -> None:
        identity = monomial_identity_obstruction(power, coeff, ETA_STRICT, BETA_STRICT)
        residuals = []
        for d in sample_ds:
            lhs = 1.0 + float(coeff) * (float(d) ** (power.numerator / power.denominator))
            rhs = 1.0 + float(BETA_STRICT) * (float(d) ** (ETA_STRICT.numerator / ETA_STRICT.denominator))
            residuals.append({"d": str(d), "candidate_minus_strict": lhs - rhs})
        models.append({
            "name": name,
            "candidate_denominator": f"1 + ({coeff}) * d^({power})",
            "required_atoms": required_atoms,
            "identity_check": identity,
            "sample_residuals": residuals,
            "max_abs_sample_residual": max(abs(r["candidate_minus_strict"]) for r in residuals),
            "accepted_as_unconditional_completion_source": False,
            "accepted_as_conditional_completion_schema": verdict == "conditional_accept",
            "verdict": verdict,
            "note": note,
        })

    add_model(
        "legacy_scalar_linear_torsion",
        Fraction(1, 1),
        BETA_TORS,
        [],
        "reject",
        "Legacy linear torsion keeps power 1 and coefficient beta_tors; it fails both strict exponent and strict beta.",
    )
    add_model(
        "scalar_beta_renormalization_without_projection",
        Fraction(1, 1),
        BETA_STRICT,
        ["positive_beta_renormalization_source"],
        "reject",
        "Changing only beta still leaves the power linear, so it is not nonlinear strict damping.",
    )
    add_model(
        "fractal_measure_pushforward_only",
        ETA_STRICT,
        BETA_TORS,
        ["fractal_measure_projection_source", "scale_domain_typing"],
        "partial",
        "The pushforward d -> d^(9/5) supplies the nonlinear exponent but leaves beta=beta_tors rather than beta=1.",
    )
    add_model(
        "fractal_pushforward_plus_independent_Z_beta",
        ETA_STRICT,
        BETA_TORS * Z_REQUIRED,
        ["fractal_measure_projection_source", "positive_beta_renormalization_source", "scale_domain_typing"],
        "conditional_accept",
        "This is exactly strict damping if an independent positive source exports Z_beta=beta/beta_tors=100; without that source it is a parametrization, not a proof.",
    )
    add_model(
        "scale_dependent_beta_eff_smuggling",
        ETA_STRICT,
        BETA_STRICT,
        ["scale_dependent_beta_eff_source"],
        "reject_as_untyped_shortcut",
        "Writing beta_eff(d)*d = beta*d^(9/5) is algebraically exact but just hides the strict law in beta_eff unless beta_eff has an independent field equation/source.",
    )
    return models


def source_lattice() -> dict[str, Any]:
    atoms = ["fractal_measure_projection_source", "positive_beta_renormalization_source", "scale_domain_typing"]
    rows = []
    for values in itertools.product([False, True], repeat=len(atoms)):
        assignment = dict(zip(atoms, values))
        accepts = all(assignment.values())
        rows.append({
            "assignment": assignment,
            "nonlinear_damping_completion_source_conditionally_repaired": accepts,
            "missing_atoms": [k for k, v in assignment.items() if not v],
        })
    current_assignment = {
        "fractal_measure_projection_source": True,
        "positive_beta_renormalization_source": False,
        "scale_domain_typing": True,
    }
    current = next(row for row in rows if row["assignment"] == current_assignment)
    return {
        "atoms": atoms,
        "rows": rows,
        "row_count": len(rows),
        "accepting_row_count": sum(1 for row in rows if row["nonlinear_damping_completion_source_conditionally_repaired"]),
        "current_assignment": current_assignment,
        "current_unconditional_repair": False,
        "current_missing_atoms": current["missing_atoms"],
        "minimal_new_required_atom": "positive_beta_renormalization_source",
    }


def p2620_update(nonlinear_repaired: bool) -> dict[str, Any]:
    rows = []
    for orientation in [False, True]:
        assignment = {
            "nonlinear_damping_completion_source": nonlinear_repaired,
            "orientation_odd_selector_source": orientation,
        }
        rows.append({
            "assignment": assignment,
            "bridge_source_cut_repaired": assignment["nonlinear_damping_completion_source"] and assignment["orientation_odd_selector_source"],
        })
    return {
        "nonlinear_atom_repaired_by_p2625_unconditionally": nonlinear_repaired,
        "rows_over_remaining_orientation_atom": rows,
        "accepting_row_count": sum(1 for row in rows if row["bridge_source_cut_repaired"]),
        "note": "P2625 does not change the orientation atom; even a future nonlinear repair would still need orientation_odd_selector_source.",
    }


def theorem_export(models: list[dict[str, Any]], lattice: dict[str, Any]) -> dict[str, Any]:
    return {
        "theorem_name": "P2625_T1_nonlinear_damping_completion_source_classification",
        "exact_positive_result": (
            "A typed fractal measure pushforward q(d)=d^(9/5) plus an independent positive beta-renormalization source "
            "Z_beta=beta/beta_tors transforms 1+beta_tors*q(d) into 1+beta*d^(9/5)."
        ),
        "exact_negative_result": (
            "Neither scalar rescaling of the legacy linear denominator nor the exponent eta=9/5 alone is a completion-source theorem. "
            "For 1+a*d^p = 1+b*d^q on all d>0 with a,b>0, exact equality forces p=q and a=b."
        ),
        "required_Z_beta": str(Z_REQUIRED),
        "required_Z_beta_float": float(Z_REQUIRED),
        "current_verdict": "CONDITIONAL_SCHEMA_ONLY_NO_UNCONDITIONAL_NONLINEAR_DAMPING_COMPLETION_SOURCE",
        "why_current_not_closed": lattice["current_missing_atoms"],
        "candidate_model_count": len(models),
        "conditional_accepting_model_names": [m["name"] for m in models if m["accepted_as_conditional_completion_schema"]],
        "rejected_shortcut_names": [m["name"] for m in models if m["verdict"].startswith("reject")],
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items() if path.suffix == ".json"}
    models = candidate_models()
    lattice = source_lattice()
    theorem = theorem_export(models, lattice)
    payload = {
        "packet_id": "P2625",
        "slice_id": "S1575",
        "status": "P2625_NONLINEAR_DAMPING_COMPLETION_CLASSIFICATION_CONDITIONAL_ZBETA_SOURCE_REQUIRED_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "inherited_artifact_status": {name: art.get("status", art.get("packet_id", "UNKNOWN")) for name, art in artifacts.items()},
        "constants": {
            "beta_tors": str(BETA_TORS),
            "beta_strict": str(BETA_STRICT),
            "eta_strict": str(ETA_STRICT),
            "required_Z_beta": str(Z_REQUIRED),
        },
        "candidate_completion_models": models,
        "nonlinear_completion_source_lattice": lattice,
        "theorem_export": theorem,
        "p2620_update": p2620_update(lattice["current_unconditional_repair"]),
        "recommended_next_honest_step": {
            "target": "derive or obstruct positive_beta_renormalization_source / Z_beta from micro-supported strict dynamics",
            "reason": "P2625 isolates the exact missing damping atom after exponent pushforward: beta_tors cannot become beta=1 without an independent positive source Z_beta=100.",
            "avoid": "do not return to selector/chirality loops until the independent damping coefficient-source question is settled or explicitly postponed",
        },
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["certificate_hash"] = sha256_json({k: v for k, v in payload.items() if k != "certificate_hash"})
    return payload


def write_markdown(payload: dict[str, Any]) -> None:
    theorem = payload["theorem_export"]
    lattice = payload["nonlinear_completion_source_lattice"]
    lines = [
        "# P2625/S1575 Nonlinear damping completion-source classification",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication grep audit",
        "",
        f"Mode: `{payload['semantic_rg_antiduplication_audit']['mode']}`.",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits; samples retained in JSON certificate.")
    lines.extend([
        "",
        "## Theorem export",
        "",
        f"Positive conditional result: {theorem['exact_positive_result']}",
        "",
        f"Negative result: {theorem['exact_negative_result']}",
        "",
        f"Required coefficient source: `Z_beta={theorem['required_Z_beta']}`.",
        f"Current verdict: `{theorem['current_verdict']}`.",
        f"Current missing atoms: `{theorem['why_current_not_closed']}`.",
        "",
        "## Candidate model verdicts",
        "",
    ])
    for model in payload["candidate_completion_models"]:
        lines.append(f"- `{model['name']}`: `{model['verdict']}`; max sample residual `{model['max_abs_sample_residual']}`; {model['note']}")
    lines.extend([
        "",
        "## Source lattice",
        "",
        f"Atoms: `{lattice['atoms']}`.",
        f"Rows: `{lattice['row_count']}`; accepting rows: `{lattice['accepting_row_count']}`.",
        f"Current assignment: `{lattice['current_assignment']}`.",
        f"Minimal new required atom: `{lattice['minimal_new_required_atom']}`.",
        "",
        "## P2620 update",
        "",
        f"Unconditional nonlinear atom repaired by P2625: `{payload['p2620_update']['nonlinear_atom_repaired_by_p2625_unconditionally']}`.",
        f"Accepting bridge rows after P2625 alone: `{payload['p2620_update']['accepting_row_count']}`.",
        "",
        "## Recommended next honest step",
        "",
        f"Target: {payload['recommended_next_honest_step']['target']}.",
        f"Reason: {payload['recommended_next_honest_step']['reason']}",
        f"Avoid: {payload['recommended_next_honest_step']['avoid']}.",
        "",
        "## Negative export flags",
        "",
    ])
    for flag, value in payload["negative_export_flags"].items():
        lines.append(f"- `{flag}`: `{value}`")
    lines.append("")
    MD.write_text("\n".join(lines), encoding="utf-8")


def update_docs() -> None:
    equation_note = """
## P2625/S1575 nonlinear damping completion-source classification guard

`P2625/S1575` supplies a sharper damping-side proof boundary.  A fractal measure pushforward `q(d)=d^(9/5)` plus an independent positive coefficient source `Z_beta=beta/beta_tors=100` would exactly transform the legacy linear torsion denominator into `1+beta*d^(9/5)`, but the current repository does not yet export that coefficient source.  Scalar rescaling, `eta=9/5` alone, or a scale-dependent `beta_eff(d)` shortcut are rejected as completion proofs; P2620 therefore remains unrepaired without the independent `positive_beta_renormalization_source` and the separate orientation source.
"""
    ltotal_note = """
## P2625/S1575 nonlinear damping completion Ltotal guard

`P2625/S1575` keeps role-bearing `L_total` closed.  It classifies the nonlinear damping completion as conditional on a typed fractal pushforward and an independent positive `Z_beta` source; because that coefficient-source theorem is not exported, there is still no full bridge repair, no role-transfer rerun, no `QW-2191` discharge, and no ToE closure.
"""
    append_once(ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md", "P2625/S1575 nonlinear damping completion-source classification guard", equation_note)
    append_once(ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md", "P2625/S1575 nonlinear damping completion Ltotal guard", ltotal_note)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    update_docs()
    print(json.dumps({
        "packet_id": payload["packet_id"],
        "status": payload["status"],
        "required_Z_beta": payload["constants"]["required_Z_beta"],
        "current_missing_atoms": payload["nonlinear_completion_source_lattice"]["current_missing_atoms"],
        "p2620_accepting_rows_after_p2625": payload["p2620_update"]["accepting_row_count"],
        "certificate_hash": payload["certificate_hash"],
    }, indent=2, sort_keys=True, ensure_ascii=False))


if __name__ == "__main__":
    main()
