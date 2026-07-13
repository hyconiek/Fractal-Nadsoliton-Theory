"""P3153/S2103: axiom-branch GR/EH nonproxy coupling audit.

P3152 closed the current SM source-selection stack as conditional and
recommended pivoting to a metric/EH interface.  This packet constructs the next
missing object, `G_EH^ax`: a finite nonproxy coupling audit that asks whether the
axiom branch now supplies an Einstein-Hilbert metric source, Newton/action unit,
and stress-energy coupling.  It computes exact FRW Einstein-tensor residuals for
several candidate scale factors and checks current repo-backed source rows.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3153_s2103_axiom_gr_eh_nonproxy_coupling_audit.json"
MD = GEN / "p3153_s2103_axiom_gr_eh_nonproxy_coupling_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3152": GEN / "p3152_s2102_y_sm_charge_unit_normalization_obstruction.json",
    "P3151": GEN / "p3151_s2101_rsm_representation_content_source_selection_audit.json",
    "P3149": GEN / "p3149_s2099_brs_ltotal_interface_invariance_audit.json",
    "P3146": GEN / "p3146_s2096_axiom_unit_action_measure_bridge.json",
    "P3145": GEN / "p3145_s2095_strict_kernel_reverse_sm_gr_layout.json",
    "P3094": GEN / "p3094_s2044_stress_energy_metric_response_obstruction_audit.json",
    "P2686": GEN / "p2686_s1636_shared_background_nonproxy_component_residual_table.json",
}

t = sp.symbols("t", positive=True)
H = sp.symbols("H", nonzero=True)


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def frw_einstein_covariant(a: sp.Expr) -> dict[str, sp.Expr]:
    adot = sp.diff(a, t)
    addot = sp.diff(a, t, 2)
    return {
        "G_00": sp.simplify(3 * (adot / a) ** 2),
        "G_ii_covariant_each": sp.simplify(-(2 * a * addot + adot**2)),
        "Ricci_scalar": sp.simplify(6 * (addot / a + (adot / a) ** 2)),
    }


def is_zero_expr(expr: sp.Expr) -> bool:
    return bool(sp.simplify(expr) == 0)


def frw_rows() -> list[dict[str, Any]]:
    candidates = [
        ("Minkowski_static", sp.Integer(1), "zero-curvature vacuum baseline only"),
        ("radiation_like_power_p_1_2", sp.sqrt(t), "common cosmology receiver shape; needs stress-energy source"),
        ("matter_like_power_p_2_3", t ** sp.Rational(2, 3), "common cosmology receiver shape; needs stress-energy source"),
        ("linear_coasting_power_p_1", t, "curved receiver with nonzero spatial residual"),
        ("deSitter_exp_Ht", sp.exp(H * t), "cosmological-constant-like receiver; needs Lambda/source unit"),
    ]
    rows = []
    for name, a, note in candidates:
        tensor = frw_einstein_covariant(a)
        nonzero = {k: not is_zero_expr(v) for k, v in tensor.items()}
        rows.append({
            "candidate_metric": name,
            "scale_factor_a_t": str(a),
            "G_00": str(tensor["G_00"]),
            "G_ii_covariant_each": str(tensor["G_ii_covariant_each"]),
            "Ricci_scalar": str(tensor["Ricci_scalar"]),
            "nonzero_components": {k: bool(v) for k, v in nonzero.items()},
            "vacuum_EH_residual_zero": not nonzero["G_00"] and not nonzero["G_ii_covariant_each"],
            "source_needed_for_nonflat_solution": nonzero["G_00"] or nonzero["G_ii_covariant_each"],
            "note": note,
        })
    return rows


def source_interface_rows() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "P3149 local BRST matter/Higgs interface",
            "metric_bundle_exported": False,
            "stress_energy_tensor_exported": False,
            "newton_or_action_unit_exported": False,
            "nonproxy_variation_exported": False,
            "noncircular_strict_source": False,
            "verdict": "local gauge invariance receiver; no metric variation, expectation state, or G_N/hbar unit",
        },
        {
            "candidate": "P3146 Lambda_unit^ax length/time/action postulates",
            "metric_bundle_exported": False,
            "stress_energy_tensor_exported": False,
            "newton_or_action_unit_exported": True,
            "nonproxy_variation_exported": False,
            "noncircular_strict_source": False,
            "verdict": "conditional units only; no strict metric source or EH coupling theorem",
        },
        {
            "candidate": "P3094 stress-energy metric-response obstruction lane",
            "metric_bundle_exported": False,
            "stress_energy_tensor_exported": False,
            "newton_or_action_unit_exported": False,
            "nonproxy_variation_exported": False,
            "noncircular_strict_source": True,
            "verdict": "repo-backed obstruction: metric response remains missing rather than exported",
        },
        {
            "candidate": "P2686 shared-background EA/EH/ELg residual table",
            "metric_bundle_exported": True,
            "stress_energy_tensor_exported": False,
            "newton_or_action_unit_exported": False,
            "nonproxy_variation_exported": True,
            "noncircular_strict_source": True,
            "verdict": "nonproxy residual evidence exists, but EH/ELg and Bianchi-I rows remain open/nonzero",
        },
        {
            "candidate": "P3145 reverse SM/GR layout",
            "metric_bundle_exported": False,
            "stress_energy_tensor_exported": False,
            "newton_or_action_unit_exported": False,
            "nonproxy_variation_exported": False,
            "noncircular_strict_source": False,
            "verdict": "receiver layout only; not a source/coupling theorem",
        },
    ]


def build_payload() -> dict[str, Any]:
    frw = frw_rows()
    interfaces = source_interface_rows()
    accepted = [
        r for r in interfaces
        if r["metric_bundle_exported"] and r["stress_energy_tensor_exported"] and r["newton_or_action_unit_exported"] and r["nonproxy_variation_exported"] and r["noncircular_strict_source"]
    ]
    nonflat_source_needed = [r for r in frw if r["candidate_metric"] != "Minkowski_static" and r["source_needed_for_nonflat_solution"]]
    counts = {
        "frw_metric_candidates": len(frw),
        "vacuum_zero_rows": sum(r["vacuum_EH_residual_zero"] for r in frw),
        "nonflat_rows_needing_source": len(nonflat_source_needed),
        "source_interface_rows": len(interfaces),
        "accepted_full_eh_coupling_rows": len(accepted),
    }
    return {
        "status": "P3153_AXIOM_GR_EH_NONPROXY_COUPLING_AUDIT_RESIDUAL_SOURCE_GAP",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "constructed_object": {
            "name": "G_EH^ax FRW/EH nonproxy coupling audit",
            "classification": "finite_metric_residual_and_source_interface_obstruction",
            "scope": "flat FRW Einstein-tensor residuals plus repo-backed source/interface rows after P3152",
        },
        "frw_einstein_residual_rows": frw,
        "source_interface_rows": interfaces,
        "finite_theorem": {
            "name": "P3153_T1_axiom_branch_eh_coupling_source_gap",
            "statement": "The flat-FRW Einstein tensor gives an exact nonproxy residual witness: only the static Minkowski baseline has zero vacuum EH residual among the audited candidates; every nonflat candidate needs a stress-energy/cosmological source.  Current axiom-branch and strict artifacts provide local matter/gauge receivers, conditional units, and nonproxy residual scaffolds, but zero rows provide the full metric bundle, stress-energy tensor, Newton/action unit, nonproxy variation, and noncircular strict source package required for EH coupling closure.",
            "finite_counts": counts,
        },
        "decision": {
            "bounded_result": "P3153 constructs the requested GR/EH bridge object and finds a concrete source gap rather than closure: nonflat metric receivers need a source, and current artifacts do not export the full EH coupling package.",
            "why_not_strict": "Minkowski is a zero vacuum baseline but not observed GR dynamics.  Nonflat FRW rows have nonzero Einstein residuals; cancelling them would require a sourced T_mu_nu/Lambda plus units and a metric variation theorem not exported by P3149/P3146/P3094/P2686/P3145.",
            "next_honest_step": "Construct P3154 as exactly one stress-energy source candidate for the axiom branch: derive a symbolic T_mu_nu from the P3149 Higgs/matter local Lagrangian in the same convention and test conservation plus coupling dimensions.  If no noncircular state/VEV/unit source is exported, preserve the P3153 EH source-gap boundary.",
            "negative_export_flags": {
                "einstein_hilbert_coupling_exported": False,
                "stress_energy_source_exported": False,
                "newton_constant_source_exported": False,
                "unit_bearing_L_total_exported": False,
                "strict_SM_generation_exported": False,
                "strict_GR_generation_exported": False,
                "selector_closure_exported": False,
                "bridge_or_role_transfer_exported": False,
                "ToE_closure_exported": False,
            },
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = ["# P3153/S2103 axiom-branch GR/EH nonproxy coupling audit", "", f"Status: `{payload['status']}`", "", "## Constructed object", f"- `{payload['constructed_object']['name']}`", f"- Classification: `{payload['constructed_object']['classification']}`", f"- Scope: `{payload['constructed_object']['scope']}`", "", "## Finite theorem", f"`{th['name']}`: {th['statement']}", "", "## Finite counts"]
    for k, v in th["finite_counts"].items():
        lines.append(f"- `{k}`: `{v}`")
    lines.extend(["", "## FRW Einstein residual rows"])
    for row in payload["frw_einstein_residual_rows"]:
        lines.append(f"- `{row['candidate_metric']}`: `a(t)={row['scale_factor_a_t']}`, `G_00={row['G_00']}`, `G_ii={row['G_ii_covariant_each']}`, vacuum zero `{row['vacuum_EH_residual_zero']}`")
    lines.extend(["", "## Source interface rows"])
    for row in payload["source_interface_rows"]:
        lines.append(f"- `{row['candidate']}`: metric `{row['metric_bundle_exported']}`, Tmunu `{row['stress_energy_tensor_exported']}`, unit `{row['newton_or_action_unit_exported']}`, variation `{row['nonproxy_variation_exported']}`, strict source `{row['noncircular_strict_source']}`; {row['verdict']}")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Why this is not strict", payload["decision"]["why_not_strict"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3153/S2103 axiom-branch GR/EH nonproxy coupling audit", "## P3153/S2103 axiom-branch GR/EH nonproxy coupling audit\n\n`P3153/S2103` constructs `G_EH^ax`, a finite nonproxy GR/EH coupling audit after the P3150-P3152 SM source-selection stack.  Exact flat-FRW Einstein residuals show only the static Minkowski baseline has zero vacuum residual; all audited nonflat receivers need stress-energy or cosmological source data.  Current repo rows provide local gauge receivers, conditional units, and residual scaffolds, but `0` rows export the full metric bundle, `T_mu_nu`, Newton/action unit, nonproxy variation, and noncircular strict source package.  No EH coupling closure, unit-bearing `L_total`, strict GR/SM generation, bridge/role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3153/S2103 EH coupling source gap", "## P3153/S2103 EH coupling source gap\n\n`P3153/S2103` verifies that nonflat FRW metric receivers carry nonzero Einstein tensor residuals unless a stress-energy/cosmological source is supplied.  The axiom branch currently lacks the noncircular state/VEV/unit/metric-variation theorem needed to turn local SM receivers into a unit-bearing Einstein-Hilbert coupling.\n")
    append_once(AGENTS, "Current axiom GR/EH nonproxy coupling guardrail (P3153/S2103, 2026-07-13)", "## Current axiom GR/EH nonproxy coupling guardrail (P3153/S2103, 2026-07-13)\n\n- P3153 constructs `G_EH^ax`, a finite FRW/EH nonproxy coupling audit for the axiom branch after P3152.\n- Exact flat-FRW residual rows show only the static Minkowski baseline is vacuum-zero; audited nonflat metric receivers require a stress-energy or cosmological source.\n- Current artifacts do not export a full EH coupling package: metric bundle, `T_mu_nu`, Newton/action unit, nonproxy variation, and noncircular strict source are not jointly present in any row.\n- Do not promote `G_EH^ax`, local SM receivers, conditional units, or residual scaffolds to GR closure, unit-bearing `L_total`, strict SM/GR generation, selector closure, bridge/role transfer, or ToE.\n- Next honest move: P3154 should construct exactly one symbolic stress-energy candidate from the P3149 matter/Higgs local Lagrangian and test conservation plus coupling dimensions; otherwise preserve the EH source-gap boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
