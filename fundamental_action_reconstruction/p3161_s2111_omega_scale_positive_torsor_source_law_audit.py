"""P3161/S2111: Omega_scale positive torsor source-law audit.

P3160 closed the alpha_geo/pi phase-locking route and recommended exactly one
new source object outside that route.  P3161 constructs that object explicitly:
Omega_scale, a proposed strict source law selecting a representative of the
positive mass/action unit torsor Omega_M/K_dim.

The central theorem is an equivariant-section obstruction.  Current strict
inputs used in the unit chain are dimensionless/invariant under the positive
scale action, while Omega_M is a free R_{>0}-torsor.  A source law from invariant
data to a free torsor cannot be R_{>0}-equivariant: if f(x) is selected and x is
scale-invariant, equivariance requires f(x)=c*f(x) for every c>0, impossible.
Therefore any accepted Omega_scale must introduce a genuinely scale-charged
strict datum, not another dimensionless alpha/phase/entropy/spectral number.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3161_s2111_omega_scale_positive_torsor_source_law_audit.json"
MD = GEN / "p3161_s2111_omega_scale_positive_torsor_source_law_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "P3157": GEN / "p3157_s2107_omega_dim_mass_unit_torsor_audit.json",
    "P3158": GEN / "p3158_s2108_post_p3157_unit_source_dependency_reconciliation.json",
    "P3160": GEN / "p3160_s2110_alpha_geo_phase_locking_no_root_of_unity_theorem.json",
    "P3116": GEN / "p3116_s2066_k_dim_dimension_source_functor_audit.json",
    "P3117": GEN / "p3117_s2067_omega_dim_dimension_character_source_audit.json",
}
ALPHA_GEO = 4.0 * math.log(2.0)
A_PHI = 2.0 * math.pi / ALPHA_GEO
SCALE_FACTORS = [1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8]
GATES = (
    "new_object_not_alpha_pi_replay",
    "strict_nadsoliton_data_only",
    "nonzero_source_value",
    "positive_torsor_target",
    "scale_charged_input_exported",
    "equivariant_section_possible",
    "selects_unique_positive_representative",
    "couples_to_K_dim_or_Omega_M",
    "breaks_scale_covariance_nonconventionally",
    "unit_chain_import_free",
    "selector_bridge_ltotal_toe_free",
)


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def rg_hits() -> dict[str, Any]:
    patterns = {
        "omega_scale_frontier": r"Omega_scale|positive_torsor_source_law|Omega_M|K_dim|Omega_dim|positive scale torsor",
        "scale_orbit_obstruction": r"scale orbit|scale-covariance|scale covariance|rescaling|torsor|equivariant|unit source",
        "alpha_pi_closed": r"P3159|P3160|alpha_geo/pi|alpha_geo/\(2\*pi\)|root of unity|phase locking",
        "forbidden_imports": r"Planck|apparatus|observed light|selector closure|bridge/role-transfer|role transfer|L_total|ToE",
    }
    out: dict[str, Any] = {}
    for name, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-g", "*.py", "-g", "*.md", "-g", "*.json"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        lines = [line for line in proc.stdout.splitlines() if line]
        out[name] = {"count": len(lines), "samples": lines[:20]}
    return out


def candidate_sources() -> list[dict[str, Any]]:
    def row(candidate: str, formula: str, blocker: str, **flags: bool) -> dict[str, Any]:
        base = {gate: False for gate in GATES}
        base.update({
            "new_object_not_alpha_pi_replay": True,
            "strict_nadsoliton_data_only": True,
            "positive_torsor_target": True,
            "unit_chain_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
        })
        base.update(flags)
        base["candidate"] = candidate
        base["formula"] = formula
        base["blocker"] = blocker
        return base

    return [
        row("dimensionless_alpha_phi_section", "Omega_scale := A_phi or alpha_geo", "dimensionless invariant; P3160 forbids phase-locking as source", nonzero_source_value=True, couples_to_K_dim_or_Omega_M=True),
        row("entropy_cardinality_scale", "Omega_scale := exp(alpha_geo)=16", "cardinality is positive but invariant under unit rescaling", nonzero_source_value=True),
        row("z12_laplacian_gap_scale", "Omega_scale := first nonzero Z12 Laplacian gap", "finite spectral number is dimensionless unless a metric/unit is already supplied", nonzero_source_value=True),
        row("damping_tail_ratio_scale", "Omega_scale := beta/tail-ratio normalizer", "P2692/P2957 leave beta scale target-dependent or orbit-gauge", nonzero_source_value=True, couples_to_K_dim_or_Omega_M=True),
        row("higgs_mu2_formal_scale", "Omega_scale := sqrt(mu2/alpha_geo)", "uses the formal unknown it is supposed to source", nonzero_source_value=True, positive_torsor_target=True, couples_to_K_dim_or_Omega_M=True),
        row("cohomology_period_scale", "Omega_scale := nonzero Z12 cocycle period", "period is a count/phase period, not a scale-charged physical unit", nonzero_source_value=True),
        row("spectral_action_logdet_scale", "Omega_scale := exp(logdet(L+m^2))", "logdet branch has dimensionless coupling-scale orbit and no physical unit source", nonzero_source_value=True),
        row("receiver_readout_unit_scale", "Omega_scale := U_readout from receiver density", "readout unit is a placeholder/gauge representative, not a strict source", nonzero_source_value=True, strict_nadsoliton_data_only=True),
        row("planck_imported_scale", "Omega_scale := m_Pl or hbar/c template", "external Planck template import", nonzero_source_value=True, strict_nadsoliton_data_only=False, unit_chain_import_free=False),
        row("apparatus_calibrated_scale", "Omega_scale := detector calibration", "apparatus import", nonzero_source_value=True, strict_nadsoliton_data_only=False, unit_chain_import_free=False),
        row("selector_chosen_scale", "Omega_scale := selector-fixed representative", "selector replay is blocked", nonzero_source_value=True, selector_bridge_ltotal_toe_free=False),
        row("hypothetical_scale_charged_source", "Omega_scale := value of new strict scale-charged datum S_+", "admissible shape only; no current artifact exports S_+", nonzero_source_value=False, scale_charged_input_exported=False, equivariant_section_possible=True, selects_unique_positive_representative=False, couples_to_K_dim_or_Omega_M=True, breaks_scale_covariance_nonconventionally=True),
    ]


def equivariant_obstruction_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for cand in candidates:
        for c in SCALE_FACTORS:
            invariant_input = not cand["scale_charged_input_exported"]
            selected = cand["nonzero_source_value"]
            transformed = c if selected else None
            residual = 0.0 if (cand["scale_charged_input_exported"] and cand["equivariant_section_possible"] and selected) else (abs(transformed - 1.0) if selected else None)
            rows.append({
                "candidate": cand["candidate"],
                "scale_factor_c": c,
                "input_fixed_by_scale_action": invariant_input,
                "selected_representative_normalized": 1.0 if selected else None,
                "torsor_action_c_times_representative": transformed,
                "equivariance_residual_abs": residual,
                "equivariant_section_obstructed": bool(invariant_input and c != 1 and selected),
                "blocker": cand["blocker"],
            })
    return rows


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for cand in candidates:
        for gate in GATES:
            rows.append({"candidate": cand["candidate"], "gate": gate, "passed": bool(cand[gate]), "detail": "passed" if cand[gate] else cand["blocker"]})
    return rows


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for candidate in sorted({row["candidate"] for row in gates}):
        sub = [row for row in gates if row["candidate"] == candidate]
        out.append({"candidate": candidate, "passed_gates": sum(row["passed"] for row in sub), "failed_gates": sum(not row["passed"] for row in sub), "accepted_Omega_scale_source": all(row["passed"] for row in sub)})
    return out


def theorem_steps() -> list[dict[str, str]]:
    return [
        {"step": "T1", "claim": "Let X be the currently exported strict input package for the unit chain; all available alpha/phase/entropy/spectral candidates are invariant under positive unit rescaling.", "status": "premise_from_P3116_P3158_P3160"},
        {"step": "T2", "claim": "Let T=Omega_M be a free R_{>0}-torsor with action t -> c t.", "status": "torsor_definition"},
        {"step": "T3", "claim": "If f:X->T is an R_{>0}-equivariant source law and X is fixed by the scale action, then f(x)=f(c.x)=c.f(x) for every c>0.", "status": "equivariance_equation"},
        {"step": "T4", "claim": "The free positive torsor has no element fixed by every c>0, so such f cannot exist.", "status": "obstruction"},
        {"step": "T5", "claim": "An accepted Omega_scale therefore requires a new scale-charged strict datum S_+ or a nonconventional scale-covariance-breaking source theorem, neither exported by current artifacts.", "status": "frontier"},
    ]


def build_payload() -> dict[str, Any]:
    p3157 = load_json(INPUTS["P3157"])
    p3158 = load_json(INPUTS["P3158"])
    p3160 = load_json(INPUTS["P3160"])
    candidates = candidate_sources()
    obstructions = equivariant_obstruction_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Omega_scale_source"]]
    return {
        "status": "P3161_OMEGA_SCALE_POSITIVE_TORSOR_SOURCE_LAW_EQUIVARIANT_OBSTRUCTION",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "input_statuses": {"P3157": p3157.get("status"), "P3158": p3158.get("status"), "P3160": p3160.get("status")},
        "content_grep": rg_hits(),
        "constructed_theoretical_objects": {
            "audit_object": "OmegaScalePositiveTorsorSourceLawAudit",
            "target": "select a strict positive representative of Omega_M/K_dim without alpha_geo/pi replay",
            "exact_equivariant_obstruction_steps": theorem_steps(),
            "candidate_Omega_scale_sources": candidates,
            "scale_action_factors": SCALE_FACTORS,
            "equivariant_obstruction_rows": obstructions,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "candidate_Omega_scale_sources": len(candidates),
            "scale_action_factors": len(SCALE_FACTORS),
            "equivariant_obstruction_rows": len(obstructions),
            "candidate_gate_rows": len(gates),
            "accepted_Omega_scale_sources": len(accepted),
            "obstructed_invariant_rows": sum(row["equivariant_section_obstructed"] for row in obstructions),
            "alpha_geo": ALPHA_GEO,
            "A_phi": A_PHI,
        },
        "finite_theorem": {
            "name": "P3161_T1_no_Omega_scale_from_invariant_strict_data",
            "statement": "No currently exported dimensionless/invariant strict datum can source-select a representative of the positive Omega_M/K_dim torsor by an equivariant source law.  Invariant input plus free positive torsor action forces f(x)=c f(x) for all c>0, impossible.  All concrete current candidates either remain dimensionless/orbit-gauge, import external calibration, replay selector/closed lanes, or merely name a hypothetical scale-charged source not exported by artifacts.",
        },
        "decision": {
            "bounded_result": "P3161 constructs the missing Omega_scale source-law object and proves a sharper obstruction: alpha/pi, entropy, finite spectra, damping ratios, logdets, readout placeholders, and formal Higgs scales cannot select the positive mass/action unit because they are scale-invariant or imported.  The admissible shape is now explicit: a new strict scale-charged datum S_+ with a nonzero value and coupling to Omega_M/K_dim.",
            "next_honest_step": "The next honest move is not another alpha_geo/pi or dimensionless invariant scan.  It must either construct one genuinely scale-charged strict source datum S_+ and test the Omega_scale gates again, or pivot to the other open leaf, Lambda_origin_source_localizer, with a non-alpha/pi origin and explicit Phi_Info/A_phi coupling.  Without S_+ or Lambda_origin, preserve the no-strict-unit/no-new-live-frontier certificate.",
            "negative_export_flags": {key: False for key in ["Omega_scale_source_exported", "Omega_M_scale_fixed", "K_dim_functor_exported", "Omega_dim_source_exported", "physical_unit_source_exported", "Higgs_VEV_source_exported", "EH_coupling_exported", "QW_2191_discharged", "selector_closure_exported", "bridge_completion_exported", "role_transfer_exported", "unit_bearing_L_total_exported", "ToE_closure_exported"]},
            "positive_scoped_flags": {"Omega_scale_object_constructed": True, "equivariant_torsor_obstruction_proved": True, "candidate_source_matrix_built": True, "next_source_shape_S_plus_identified": True, "repo_grep_performed": True},
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = [
        "# P3161/S2111 Omega_scale positive torsor source-law audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact obstruction",
    ]
    for step in payload["constructed_theoretical_objects"]["exact_equivariant_obstruction_steps"]:
        lines.append(f"- `{step['step']}` `{step['status']}`: {step['claim']}")
    lines.extend([
        "",
        "## Finite certificate",
        f"- candidate Omega_scale sources: `{cert['candidate_Omega_scale_sources']}`",
        f"- scale action factors: `{cert['scale_action_factors']}`",
        f"- equivariant obstruction rows: `{cert['equivariant_obstruction_rows']}`",
        f"- candidate gate rows: `{cert['candidate_gate_rows']}`",
        f"- accepted Omega_scale sources: `{cert['accepted_Omega_scale_sources']}`",
        f"- obstructed invariant rows: `{cert['obstructed_invariant_rows']}`",
        f"- A_phi retained: `{cert['A_phi']}`",
        "",
        "## Finite theorem",
        f"`{payload['finite_theorem']['name']}`: {payload['finite_theorem']['statement']}",
        "",
        "## Decision",
        payload["decision"]["bounded_result"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3161/S2111 Omega_scale positive torsor source-law audit", "## P3161/S2111 Omega_scale positive torsor source-law audit\n\n`P3161/S2111` constructs the missing `Omega_scale` source-law audit for selecting a representative of the positive `Omega_M/K_dim` torsor after P3160 closed alpha/pi phase locking.  The exact equivariant-section obstruction says that invariant strict data cannot map equivariantly into a free positive torsor: `f(x)=f(c.x)=c f(x)` for all `c>0` is impossible.  It audits `12` candidate source classes over `7` scale factors, builds `84` obstruction rows and `132` gate rows, and finds `0` accepted `Omega_scale` sources.  No strict physical unit, Higgs VEV, EH coupling, selector closure, bridge/role-transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3161/S2111 Omega_scale torsor source remains unsourced", "## P3161/S2111 Omega_scale torsor source remains unsourced\n\n`P3161/S2111` proves that current dimensionless strict data cannot equivariantly select the positive mass/action unit torsor.  A Lagrangian/action-unit source would require a genuinely scale-charged strict datum `S_+` or an equivalent nonconventional scale-covariance-breaking theorem, not another alpha/pi, entropy, spectral, damping, logdet, readout, Planck, apparatus, or selector representative.\n")
    append_once(AGENTS, "Current Omega_scale positive torsor source-law guardrail (P3161/S2111, 2026-07-13)", "## Current Omega_scale positive torsor source-law guardrail (P3161/S2111, 2026-07-13)\n\n- P3161 constructs the `Omega_scale` positive torsor source-law audit requested after P3160, targeting a strict representative of `Omega_M/K_dim` without alpha_geo/pi replay.\n- The exact equivariant-section obstruction shows that invariant strict data cannot select a free positive torsor representative: `f(x)=f(c.x)=c f(x)` for all `c>0` is impossible.\n- The audit tests `12` candidate source classes over `7` scale factors, with `84` obstruction rows and `132` gate rows; `0` candidates export an import-free strict `Omega_scale` source.\n- Do not promote alpha/pi sections, entropy/cardinality, Z12 spectra, damping ratios, formal Higgs scales, cohomology periods, logdet actions, readout placeholders, Planck/apparatus imports, selector choices, or hypothetical unexported `S_+` to physical units, Higgs VEV, EH coupling, selector closure, bridge completion, role transfer, `L_total`, or ToE.\n- Next honest move: either construct a real scale-charged strict source datum `S_+` with nonzero value and `Omega_M/K_dim` coupling, or pivot to `Lambda_origin_source_localizer`; otherwise preserve the no-strict-unit/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
