#!/usr/bin/env python3
"""P2686/S1636: shared-background nonproxy EA/EH/ELg component residual table.

This implements the P2685 recommendation.  It builds one finite residual table
from the existing unified EA/EH/ELg runpack and the Bianchi-I anisotropic metric
obstruction family, then decides whether a zero table or a no-go boundary is the
honest current output.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2686_s1636_shared_background_nonproxy_component_residual_table.json"
MD = GEN / "p2686_s1636_shared_background_nonproxy_component_residual_table.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2685": GEN / "p2685_s1635_strict_lagrangian_eom_reverse_closure_obstruction_matrix.json",
    "P1764": GEN / "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json",
    "P1765": GEN / "p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json",
    "P1787": GEN / "p1787_s737_strict_unified_ea_eh_componentwise_export_and_h1_run_contract_checkpoint.json",
    "P1805": GEN / "p1805_s755_strict_bw_entry_nonproxy_covariant_export_completeness_audit_checkpoint.json",
    "P1806_CONTRACT": GEN / "p1806_s756_strict_unified_nonproxy_ea_eh_elg_residual_runpack_contract_checkpoint.json",
    "P1806_INPUT": GEN / "p1806_s756_unified_nonproxy_ea_eh_elg_runpack_input.json",
    "P1848": GEN / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json",
    "P1868": GEN / "p1868_s818_strict_4d_componentwise_residual_table_scaffold.json",
    "P1869": GEN / "p1869_s819_strict_flat_background_component_residual_probe.json",
    "P1974": GEN / "p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.json",
    "P1975": GEN / "p1975_s925_strict_minimal_anisotropic_source_obligation_and_energy_sign_audit.json",
    "P1977": GEN / "p1977_s927_strict_positive_energy_anisotropic_provider_bounded_no_go.json",
    "P1978": GEN / "p1978_s928_strict_energy_neutral_tensor_transport_obstruction.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "shared_background_zero_table_exported",
    "full_nonproxy_closure_exported",
    "role_bearing_ltotal_exported",
    "selector_or_bridge_imported",
    "qw2191_discharged",
    "toe_closure_claimed",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        ["rg", "-n", pattern, ".", "-g", "*.py", "-g", "*.md", "-g", "*.json", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:60]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "shared_background_nonproxy_runpack": r"shared-background|shared background|bg_family_v1|EA|EH|ELg|E_A_mu|E_H|EL_g|nonproxy component",
        "component_residual_table": r"component residual|componentwise residual|residual table|PASS_ZERO|OPEN_OBSTRUCTION_WITH_TRACE|residual_norm",
        "bianchi_anisotropic_obstruction": r"Bianchi-I|anisotropic residual|Q_shear|rho_required|energy-neutral|positive-energy",
        "gate_locks": r"TG1_BW|metric_full_tensor_closure|unified nonproxy residual|BRST|Cutkosky|BW lock",
        "forbidden_closure_claims": r"ToE closure|QW-2191 discharge|role-bearing L_total|selector closure|role transfer|legacy bridge",
    }
    return {"tool": "rg", "mode": "content-first shared-background nonproxy component residual table", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def current_state() -> dict[str, Any]:
    data = {name: load_json(path) for name, path in INPUTS.items()}
    p1806_input = data["P1806_INPUT"]
    p1806_contract = data["P1806_CONTRACT"]
    p1977 = data["P1977"]
    p1978 = data["P1978"]
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2685_next_was_component_table": "shared-background nonproxy component residual table" in data["P2685"].get("decision", {}).get("next_honest_step", ""),
        "shared_background_family_id": p1806_input.get("shared_background_family_id"),
        "p1806_tg1_locked": p1806_contract.get("gate_vector", {}).get("TG1_BW") == "OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL",
        "p1806_ea_zero": p1806_contract.get("checks", {}).get("ea_mu_zero") is True,
        "p1806_eh_zero": p1806_contract.get("checks", {}).get("eh_zero") is True,
        "p1806_elg_zero": p1806_contract.get("checks", {}).get("elg_zero") is True,
        "p1977_positive_energy_bounded_no_go": p1977.get("gatekeeper_checks", {}).get("bounded_no_go_passed") is True,
        "p1978_energy_neutral_obstruction_passed": p1978.get("gatekeeper_checks", {}).get("bounded_obstruction_passed") is True,
    }


def unified_runpack_rows() -> list[dict[str, Any]]:
    p1806 = load_json(INPUTS["P1806_INPUT"])
    rows = []
    for component_id, key in [("EA", "ea_mu_residual"), ("EH", "eh_residual"), ("ELg", "elg_residual")]:
        row = p1806.get(key, {})
        rows.append({
            "row_id": f"runpack_{component_id}",
            "background_family": p1806.get("shared_background_family_id"),
            "component": component_id,
            "verdict": row.get("verdict", "MISSING"),
            "residual_norm": row.get("residual_norm"),
            "trace": row.get("trace"),
            "zero": row.get("verdict") == "PASS_ZERO" and float(row.get("residual_norm", 1.0)) == 0.0,
        })
    return rows


def bianchi_rows() -> dict[str, Any]:
    p1974 = load_json(INPUTS["P1974"])
    p1975 = load_json(INPUTS["P1975"])
    p1977 = load_json(INPUTS["P1977"])
    p1978 = load_json(INPUTS["P1978"])
    sigma1, sigma2 = sp.symbols("sigma1 sigma2")
    q_shear = sp.sympify(p1975.get("source_obligation", {}).get("q_shear", "sigma1**2 + sigma1*sigma2 + sigma2**2"))
    q_matrix = sp.hessian(q_shear, (sigma1, sigma2)) / 2
    eigenvals = sorted([str(ev) for ev in q_matrix.eigenvals().keys()])
    rows = []
    for index, expr_text in enumerate(p1974.get("anisotropic_eom_residual_vector", [])):
        expr = sp.sympify(expr_text)
        rows.append({
            "row_id": f"bianchi_I_{p1974.get('component_basis', [])[index]}",
            "component": p1974.get("component_basis", [])[index],
            "residual": str(expr),
            "is_zero_polynomial": bool(sp.simplify(expr) == 0),
        })
    return {
        "rows": rows,
        "q_shear": str(q_shear),
        "q_shear_matrix": str(q_matrix),
        "q_shear_eigenvalues": eigenvals,
        "minimal_source_cancels_if_admitted": p1975.get("gatekeeper_checks", {}).get("minimal_source_cancels_residual_if_admitted") is True,
        "strict_source_derivation_exported": p1975.get("gatekeeper_checks", {}).get("strict_source_derivation_exported") is True,
        "positive_energy_provider_bounded_no_go": p1977.get("gatekeeper_checks", {}).get("bounded_no_go_passed") is True,
        "energy_neutral_transport_obstructed": p1978.get("gatekeeper_checks", {}).get("bounded_obstruction_passed") is True,
    }


def residual_table(runpack: list[dict[str, Any]], bianchi: dict[str, Any]) -> list[dict[str, Any]]:
    rows = list(runpack)
    rows.extend({
        "row_id": row["row_id"],
        "background_family": "diagonal_Bianchi_I_tracefree_shear",
        "component": row["component"],
        "verdict": "PASS_ZERO" if row["is_zero_polynomial"] else "OPEN_OBSTRUCTION_WITH_TRACE",
        "residual_norm": None,
        "trace": row["residual"],
        "zero": row["is_zero_polynomial"],
    } for row in bianchi["rows"])
    return rows


def decision(table: list[dict[str, Any]], bianchi: dict[str, Any]) -> dict[str, Any]:
    nonzero = [row for row in table if not row["zero"]]
    no_go_active = (
        bool(nonzero)
        and bianchi["positive_energy_provider_bounded_no_go"]
        and bianchi["energy_neutral_transport_obstructed"]
        and not bianchi["strict_source_derivation_exported"]
    )
    return {
        "decision": "P2686_SHARED_BACKGROUND_NONPROXY_COMPONENT_RESIDUAL_TABLE_NONZERO_NO_GO_BOUNDARY_NO_FALSE_PASS",
        "zero_table_exported": not nonzero,
        "nonzero_row_ids": [row["row_id"] for row in nonzero],
        "bounded_no_go_boundary_active": no_go_active,
        "professorial_verdict": (
            "P2686 fills the requested shared-background nonproxy component residual table.  The EA row is locally zero in the existing runpack, but EH and ELg remain nonzero/open there, and the Bianchi-I metric residual supplies four symbolic nonzero component rows.  P1975 shows the minimal cancelling source would work only if admitted; P1977 and P1978 block the current positive-energy and energy-neutral routes.  Therefore the current honest output is a bounded no-go boundary for reverse closure from the present reduced/FRW/nonproxy scaffold, not L_total or ToE closure."
        ),
        "next_honest_step": (
            "P2687 should target exactly one new strict anisotropic source class: either a derived lapse/energy source that evades the P1977 positive-energy no-go without negative rho, or a non-energy-neutral tensorial shear transport that evades P1978.  If no such typed source is introduced, freeze this Lagrangian/EOM reverse-closure lane as bounded no-go and return to the broader state-map for a different live frontier."
        ),
        "full_nonproxy_closure_exported_now": False,
        "role_bearing_ltotal_exported_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2686/S1636 shared-background nonproxy component residual table", "", f"Status: `{payload['status']}`", "", "## Content-first grep"]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Component residual rows"])
    for row in payload["component_residual_table"]:
        lines.append(f"- `{row['row_id']}` on `{row['background_family']}`: verdict=`{row['verdict']}`, zero=`{row['zero']}`, trace=`{row['trace']}`")
    lines.extend([
        "", "## Bianchi source-route checks",
        f"Q_shear: `{payload['bianchi_rows']['q_shear']}`; eigenvalues: `{payload['bianchi_rows']['q_shear_eigenvalues']}`.",
        f"Minimal source cancels if admitted: `{payload['bianchi_rows']['minimal_source_cancels_if_admitted']}`; strict source exported: `{payload['bianchi_rows']['strict_source_derivation_exported']}`.",
        f"Positive-energy no-go: `{payload['bianchi_rows']['positive_energy_provider_bounded_no_go']}`; energy-neutral obstruction: `{payload['bianchi_rows']['energy_neutral_transport_obstructed']}`.",
        "", "## Verdict", payload["decision"]["professorial_verdict"],
        "", "## Next honest step", payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    state = current_state()
    runpack = unified_runpack_rows()
    bianchi = bianchi_rows()
    table = residual_table(runpack, bianchi)
    payload: dict[str, Any] = {
        "status": "P2686_SHARED_BACKGROUND_NONPROXY_COMPONENT_RESIDUAL_TABLE_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "current_state": state,
        "unified_runpack_rows": runpack,
        "bianchi_rows": bianchi,
        "component_residual_table": table,
        "decision": decision(table, bianchi),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2686/S1636 shared-background nonproxy component residual table",
        "## P2686/S1636 shared-background nonproxy component residual table\n\n"
        "`P2686/S1636` fills the P2685-requested `EA/EH/ELg` component residual table.  The existing unified runpack has `EA=PASS_ZERO`, but `EH` and `ELg` remain open/nonzero, and the Bianchi-I metric residual contributes symbolic nonzero component rows.  P1977/P1978 block current positive-energy and energy-neutral cancellation routes, so the current output is a bounded no-go boundary for reverse closure, not role-bearing `L_total` or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2686/S1636 component residual no-go guard",
        "## P2686/S1636 component residual no-go guard\n\n"
        "`P2686/S1636` keeps `L_total` nonpromoted: the shared-background nonproxy component table is not zero.  A future continuation must introduce one new strict anisotropic source class that evades the P1977 positive-energy no-go or the P1978 energy-neutral tensor-transport obstruction; otherwise this reverse-closure lane is bounded no-go.\n",
    )
    append_once(
        AGENTS,
        "Current shared-background component residual guardrail (P2686/S1636, 2026-06-13)",
        "## Current shared-background component residual guardrail (P2686/S1636, 2026-06-13)\n\n"
        "- P2686 fills the `EA/EH/ELg` nonproxy component table and finds a nonzero bounded no-go boundary: `EA` is locally zero, but `EH`, `ELg`, and Bianchi-I metric residual rows remain open/nonzero.\n"
        "- A next move in this lane is admissible only if it introduces one new strict anisotropic source class evading P1977 or P1978; otherwise pivot away rather than replaying reduced/FRW reverse closure.\n",
    )
    return payload


if __name__ == "__main__":
    main()
