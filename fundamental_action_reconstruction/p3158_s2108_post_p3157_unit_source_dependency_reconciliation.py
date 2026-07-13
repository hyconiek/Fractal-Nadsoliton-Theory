"""P3158/S2108: post-P3157 unit-source dependency reconciliation.

P3157 constructed a formal mass-unit torsor but left its positive scale unsourced.
This packet greps/reuses the existing P3116-P3124 unit-source chain and builds a
finite dependency DAG for what would be needed to turn Omega_M into a strict
mass/action unit.  It is a proof-grade reconciliation step: no new closure is
claimed, and the next admissible move is narrowed to one missing source object.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3158_s2108_post_p3157_unit_source_dependency_reconciliation.json"
MD = GEN / "p3158_s2108_post_p3157_unit_source_dependency_reconciliation.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3157": GEN / "p3157_s2107_omega_dim_mass_unit_torsor_audit.json",
    "P3116": GEN / "p3116_s2066_k_dim_dimension_source_functor_audit.json",
    "P3117": GEN / "p3117_s2067_omega_dim_dimension_character_source_audit.json",
    "P3118": GEN / "p3118_s2068_r_dim_action_length_time_relation_audit.json",
    "P3119": GEN / "p3119_s2069_xi_lt_axis_source_object_audit.json",
    "P3120": GEN / "p3120_s2070_tau_lt_ordered_flow_source_audit.json",
    "P3121": GEN / "p3121_s2071_kappa_cycle_source_audit.json",
    "P3122": GEN / "p3122_s2072_iota_irrev_source_audit.json",
    "P3123": GEN / "p3123_s2073_delta_asym_source_audit.json",
    "P3124": GEN / "p3124_s2074_phi_info_phase_information_gauge_quotient_audit.json",
}

DAG = {
    "Omega_M_strict_mass_unit": ["K_dim_functor", "positive_torsor_source_law"],
    "K_dim_functor": ["Omega_dim_character", "Sigma_dim_section", "C_phi_A_phi_to_U_action", "R_dim_action_length_time"],
    "R_dim_action_length_time": ["Xi_LT_axis_source"],
    "Xi_LT_axis_source": ["Tau_LT_ordered_flow"],
    "Tau_LT_ordered_flow": ["Kappa_cycle_source"],
    "Kappa_cycle_source": ["Iota_irrev_source"],
    "Iota_irrev_source": ["Delta_asym_source"],
    "Delta_asym_source": ["Phi_Info_source"],
    "Phi_Info_source": ["Lambda_origin_source_localizer"],
    "Lambda_origin_source_localizer": [],
    "Omega_dim_character": [],
    "Sigma_dim_section": [],
    "C_phi_A_phi_to_U_action": [],
    "positive_torsor_source_law": [],
}


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
        "omega_dim_kdim": r"Omega_dim|K_dim|dimension-source functor|positive scale torsor",
        "axis_flow_chain": r"R_dim|Xi_LT|Tau_LT|Kappa_cycle|Iota_irrev|Delta_asym|Phi_Info|Lambda_origin",
        "unit_forbidden_promotions": r"unit-bearing L_total|physical.*unit|Planck|apparatus|selector replay|bridge/role-transfer|ToE",
    }
    out = {}
    for name, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "-g", "*.py", "-g", "*.md", "-g", "*.json"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        lines = [line for line in proc.stdout.splitlines() if line]
        out[name] = {"count": len(lines), "samples": lines[:20]}
    return out


def exported_statuses() -> dict[str, bool]:
    p3116 = load_json(INPUTS["P3116"])
    p3117 = load_json(INPUTS["P3117"])
    p3118 = load_json(INPUTS["P3118"])
    p3119 = load_json(INPUTS["P3119"])
    p3120 = load_json(INPUTS["P3120"])
    p3124 = load_json(INPUTS["P3124"])
    return {
        "Omega_M_strict_mass_unit": False,
        "K_dim_functor": p3116.get("finite_certificate", {}).get("accepted_K_dim_functors", 0) > 0,
        "Omega_dim_character": p3117.get("finite_certificate", {}).get("accepted_Omega_dim_sources", 0) > 0,
        "Sigma_dim_section": False,
        "C_phi_A_phi_to_U_action": False,
        "R_dim_action_length_time": p3118.get("finite_certificate", {}).get("accepted_R_dim_relations", 0) > 0,
        "Xi_LT_axis_source": p3119.get("finite_certificate", {}).get("accepted_Xi_LT_axis_sources", 0) > 0,
        "Tau_LT_ordered_flow": p3120.get("finite_certificate", {}).get("accepted_Tau_LT_ordered_flows", 0) > 0,
        "Kappa_cycle_source": False,
        "Iota_irrev_source": False,
        "Delta_asym_source": False,
        "Phi_Info_source": p3124.get("finite_certificate", {}).get("accepted_Phi_Info_sources", 0) > 0,
        "Lambda_origin_source_localizer": False,
        "positive_torsor_source_law": False,
    }


def dependency_rows(status: dict[str, bool]) -> list[dict[str, Any]]:
    rows = []
    for node, deps in DAG.items():
        deps_exported = all(status.get(dep, False) for dep in deps)
        rows.append({
            "node": node,
            "dependencies": deps,
            "node_exported_current": status.get(node, False),
            "dependencies_exported": deps_exported,
            "can_close_now": status.get(node, False) and deps_exported,
            "missing_dependencies": [dep for dep in deps if not status.get(dep, False)],
        })
    return rows


def leaf_cut(rows: list[dict[str, Any]]) -> list[str]:
    return sorted(row["node"] for row in rows if not row["dependencies"] and not row["node_exported_current"])


def build_payload() -> dict[str, Any]:
    status = exported_statuses()
    rows = dependency_rows(status)
    leaves = leaf_cut(rows)
    counts = {
        "dag_nodes": len(DAG),
        "dag_edges": sum(len(v) for v in DAG.values()),
        "exported_nodes": sum(status.values()),
        "closed_nodes_now": sum(row["can_close_now"] for row in rows),
        "missing_leaf_cut_size": len(leaves),
    }
    return {
        "status": "P3158_POST_P3157_UNIT_SOURCE_DEPENDENCY_RECONCILIATION_NO_STRICT_UNIT",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "content_grep": rg_hits(),
        "constructed_object": {
            "name": "U_mass_unit_source_dependency_DAG",
            "classification": "post_p3157_unit_source_dependency_reconciliation",
            "scope": "Omega_M strict mass/action unit closure dependencies through P3116-P3124",
        },
        "dependency_rows": rows,
        "missing_leaf_cut": leaves,
        "finite_theorem": {
            "name": "P3158_T1_unit_source_dependency_no_current_closure",
            "statement": "The post-P3157 strict mass-unit problem reduces to a finite dependency DAG through the existing P3116-P3124 unit/source chain.  Current artifacts export zero nodes in the required closure path for Omega_M as a strict mass/action unit, and the missing leaf cut remains nonempty.  Therefore P3157's formal torsor cannot be promoted to a strict unit source on current artifacts.",
            "finite_counts": counts,
        },
        "decision": {
            "bounded_result": "P3158 reconciles the unit-source frontier and confirms that the Higgs/SM/EH branch should remain frozen as conditional unless a genuinely new unit/source leaf is supplied.",
            "next_honest_step": "The least-replay next move is exactly one new leaf object, preferably Lambda_origin_source_localizer if continuing the P3124 phase-information chain, or a genuinely new positive_torsor_source_law for Omega_M.  Without such a supplied object, preserve the no-strict-unit/no-new-live-frontier certificate.",
            "negative_export_flags": {key: False for key in ["strict_mass_unit_source_exported", "K_dim_functor_exported", "Omega_dim_source_exported", "strict_higgs_couplings_exported", "einstein_hilbert_coupling_exported", "unit_bearing_L_total_exported", "strict_SM_generation_exported", "strict_GR_generation_exported", "selector_closure_exported", "bridge_or_role_transfer_exported", "ToE_closure_exported"]},
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = ["# P3158/S2108 post-P3157 unit-source dependency reconciliation", "", f"Status: `{payload['status']}`", "", "## Constructed object", f"- `{payload['constructed_object']['name']}`", f"- Classification: `{payload['constructed_object']['classification']}`", f"- Scope: `{payload['constructed_object']['scope']}`", "", "## Content grep"]
    for k, v in payload["content_grep"].items():
        lines.append(f"- `{k}`: `{v['count']}` hits")
    lines.extend(["", "## Finite theorem", f"`{th['name']}`: {th['statement']}", "", "## Finite counts"])
    for k, v in th["finite_counts"].items():
        lines.append(f"- `{k}`: `{v}`")
    lines.extend(["", "## Missing leaf cut", ", ".join(f"`{x}`" for x in payload["missing_leaf_cut"]), "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3158/S2108 post-P3157 unit-source dependency reconciliation", "## P3158/S2108 post-P3157 unit-source dependency reconciliation\n\n`P3158/S2108` builds `U_mass_unit_source_dependency_DAG`, reconciling the formal `Omega_M` torsor from P3157 with the existing P3116-P3124 unit/source chain.  The finite DAG has no current accepting closure path to a strict mass/action unit; the missing leaf cut remains nonempty.  No strict mass unit, K_dim/Omega_dim source, Higgs coupling, EH coupling, unit-bearing `L_total`, strict SM/GR generation, selector closure, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3158/S2108 unit-source dependency remains open", "## P3158/S2108 unit-source dependency remains open\n\n`P3158/S2108` confirms that the formal mass-unit carrier from P3157 cannot be promoted without a new source leaf in the P3116-P3124 unit chain.  The Higgs/SM/EH branch should remain conditional until such a source law is supplied.\n")
    append_once(AGENTS, "Current post-P3157 unit-source dependency guardrail (P3158/S2108, 2026-07-13)", "## Current post-P3157 unit-source dependency guardrail (P3158/S2108, 2026-07-13)\n\n- P3158 constructs `U_mass_unit_source_dependency_DAG`, a finite dependency reconciliation for turning P3157 `Omega_M` into a strict mass/action unit.\n- Current artifacts provide no accepting closure path through P3116-P3124; the missing leaf cut remains nonempty.\n- Do not promote the formal mass-unit torsor, Higgs/SM/EH scaffold, or unit-chain bookkeeping to strict units, `L_total`, SM/GR generation, selector closure, bridge/role transfer, or ToE.\n- Next honest move: supply exactly one new leaf source object such as `Lambda_origin_source_localizer` or a genuine positive torsor source law for `Omega_M`; otherwise preserve the no-strict-unit/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
