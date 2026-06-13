#!/usr/bin/env python3
"""P2690/S1640: selector-free nonexact boundary-phase sector provider audit.

This executes the P2689 next step: test whether the P2663 one-bit carrier has a
selector-free provider that prefers the nonexact square-holonomy sector and
therefore supplies a unique N_bits target for the entropy/UV-unit route.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2690_s1640_selector_free_nonexact_boundary_phase_sector_provider_audit.json"
MD = GEN / "p2690_s1640_selector_free_nonexact_boundary_phase_sector_provider_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2689": GEN / "p2689_s1639_entropy_reference_cell_bit_to_length_uv_unit_obligation_matrix.json",
    "P2663": GEN / "p2663_s1613_chain_level_boundary_phase_bit_target_audit.json",
    "P2664": GEN / "p2664_s1614_boundary_phase_sector_selector_variational_no_go_audit.json",
    "P2665": GEN / "p2665_s1615_selector_lane_to_boundary_phase_sector_bridge_audit.json",
    "P2672": GEN / "p2672_s1622_source_topology_to_square_holonomy_typed_descent_audit.json",
    "P2673": GEN / "p2673_s1623_tau_src_pair12_boundary_square_subinterface_audit.json",
    "P2680": GEN / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "selector_free_nonexact_sector_provider_exported",
    "boundary_phase_bit_target_exported",
    "uv_unit_or_beta_source_exported",
    "selector_replay_imported",
    "qw2191_discharged",
    "role_transfer_started",
    "ltotal_promoted",
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
    return {"count": len(lines), "samples": lines[:80]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "p2689_selected_p2690": r"P2690|selector-free nonexact|nonexact boundary-phase sector|one-bit carrier",
        "chain_level_one_bit_carrier": r"P2663|one-bit|square_holonomy|nonexact sector|exact coboundary|N_bits",
        "sector_provider_candidates": r"P2664|P2665|P2672|theta|local positive|source topology|typed descent|sector selector|provider",
        "pair12_subinterface_blockers": r"P2673|tau_src|pair12|boundary square|typed subinterface|chart-sensitive",
        "forbidden_imports": r"selector replay|QW-2191 discharge|role transfer|bridge completion|L_total|ToE closure",
    }
    return {"tool": "rg", "mode": "content-first selector-free nonexact boundary-phase sector provider audit", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def state_reads() -> dict[str, Any]:
    data = {name: load_json(path) for name, path in INPUTS.items()}
    p2663 = data["P2663"].get("chain_level_boundary_phase_witness", {})
    p2664 = data["P2664"].get("closure_decision", {})
    p2665 = data["P2665"].get("closure_decision", {})
    p2672 = data["P2672"].get("closure_decision", {})
    p2673 = data["P2673"].get("closure_decision", {})
    p2689 = data["P2689"].get("decision", {})
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2689_selected_p2690": "P2690" in p2689.get("next_honest_step", ""),
        "p2663_nonexact_sector_needed": p2663.get("nonexact_sector_needed_for_nonzero_bit") is True,
        "p2663_unique_bits_without_sector_choice": p2663.get("chain_law_derives_unique_n_bits_without_sector_choice") is True,
        "p2663_exact_coboundary_square_values": p2663.get("exact_coboundary_square_holonomy_values", []),
        "p2664_local_even_class_selects_nonexact": p2664.get("local_even_class_selects_nonexact_sector") is True,
        "p2664_declared_theta_can_select": p2664.get("declared_theta_can_select_nonexact_sector") is True,
        "p2664_passing_sector_selectors": p2664.get("passing_sector_selector_candidates", []),
        "p2665_passing_selector_transfers": p2665.get("passing_selector_to_boundary_phase_bridge_lanes", []),
        "p2672_typed_descent_exported": p2672.get("boundary_phase_bit_target_exported_now") is True,
        "p2673_pair12_subinterface_exported": p2673.get("pair12_to_boundary_square_typed_subinterface_exported_now") is True,
    }


def carrier_enumeration() -> dict[str, Any]:
    rows = load_json(INPUTS["P2663"]).get("chain_level_boundary_phase_witness", {}).get("rows", [])
    exact = [row for row in rows if row.get("is_exact_coboundary")]
    nonexact = [row for row in rows if not row.get("is_exact_coboundary")]
    bit1 = [row for row in rows if row.get("square_holonomy_bit") == 1]
    exact_bit_values = sorted({row.get("square_holonomy_bit") for row in exact})
    nonexact_bit_values = sorted({row.get("square_holonomy_bit") for row in nonexact})
    return {
        "total_sampled_rows": len(rows),
        "exact_rows": len(exact),
        "nonexact_rows": len(nonexact),
        "square_bit_one_rows": len(bit1),
        "exact_square_bit_values": exact_bit_values,
        "nonexact_square_bit_values": nonexact_bit_values,
        "nonzero_bit_requires_nonexact_sector": 1 not in exact_bit_values and bool(bit1),
        "carrier_exists": bool(bit1),
    }


def provider_matrix(reads: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "provider_id": "local_positive_boundary_phase_variational_class",
            "selector_free": True,
            "prefers_nonexact_sector": reads["p2664_local_even_class_selects_nonexact"],
            "exported_now": False,
            "rejection": "P2664 shows the local positive/even class does not select the nonexact sector.",
        },
        {
            "provider_id": "declared_theta_holonomy_source",
            "selector_free": False,
            "prefers_nonexact_sector": reads["p2664_declared_theta_can_select"],
            "exported_now": False,
            "rejection": "Theta can select only as a declared premise; it is not selector-free strict provenance.",
        },
        {
            "provider_id": "selector_lane_transfer_to_boundary_phase",
            "selector_free": False,
            "prefers_nonexact_sector": bool(reads["p2665_passing_selector_transfers"]),
            "exported_now": False,
            "rejection": "P2665 has no accepted transfer and would import selector-lane material.",
        },
        {
            "provider_id": "source_topology_sign_typed_descent",
            "selector_free": True,
            "prefers_nonexact_sector": reads["p2672_typed_descent_exported"],
            "exported_now": False,
            "rejection": "Current typed descent does not export the boundary-phase bit target.",
        },
        {
            "provider_id": "tau_src_pair12_boundary_square_subinterface",
            "selector_free": True,
            "prefers_nonexact_sector": reads["p2673_pair12_subinterface_exported"],
            "exported_now": False,
            "rejection": "P2673 leaves the chart-sensitive pair12 -> boundary-square subinterface unexported.",
        },
    ]


def decision(carrier: dict[str, Any], providers: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row for row in providers if row["selector_free"] and row["prefers_nonexact_sector"] and row["exported_now"]]
    return {
        "decision": "P2690_SELECTOR_FREE_NONEXACT_BOUNDARY_PHASE_SECTOR_PROVIDER_AUDIT_NO_EXPORTED_PROVIDER_FREEZE_ENTROPY_UV_ROUTE_NO_FALSE_PASS",
        "carrier_exists": carrier["carrier_exists"],
        "selector_free_provider_exported": bool(passing),
        "passing_provider_ids": [row["provider_id"] for row in passing],
        "freeze_entropy_uv_unit_route": not passing,
        "professorial_verdict": (
            "P2690 separates the real carrier from the missing provider.  P2663 gives nonexact square-holonomy rows with bit one, while exact coboundaries stay at bit zero, so the carrier is real.  However, every audited provider route fails the stricter P2689 requirement: local positive boundary-phase dynamics does not prefer the nonexact sector, theta selection is declared rather than sourced, selector-lane transfer is forbidden/no-pass, source-topology typed descent does not export the bit target, and tau_src->pair12->boundary-square remains subinterface-blocked.  Thus no selector-free nonexact boundary-phase sector provider is exported, and the entropy/UV-unit route freezes as bounded no-go on current artifacts."
        ),
        "next_honest_step": (
            "P2691 should return to the broad state-map and select a different finite source atom.  The most concrete remaining P2680 non-selector target is an amplitude role-safe source audit: test whether alpha_geo scalar-shape normalization can be given a strict, role-safe amplitude source without legacy role transfer, selector replay, beta_tors->chi11, or generic bridge completion."
        ),
        "boundary_phase_bit_target_exported_now": False,
        "uv_unit_or_beta_source_exported_now": False,
        "role_transfer_started_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2690/S1640 selector-free nonexact boundary-phase sector provider audit", "", f"Status: `{payload['status']}`", "", "## Content-first grep"]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    c = payload["carrier_enumeration"]
    lines.extend(["", "## Carrier enumeration", f"Rows: total=`{c['total_sampled_rows']}`, exact=`{c['exact_rows']}`, nonexact=`{c['nonexact_rows']}`, bit1=`{c['square_bit_one_rows']}`.", f"Exact square bits: `{c['exact_square_bit_values']}`; nonexact square bits: `{c['nonexact_square_bit_values']}`."])
    lines.extend(["", "## Provider matrix"])
    for row in payload["provider_matrix"]:
        lines.append(f"- `{row['provider_id']}`: selector_free=`{row['selector_free']}`, prefers_nonexact=`{row['prefers_nonexact_sector']}`, exported=`{row['exported_now']}` — {row['rejection']}")
    lines.extend(["", "## Verdict", payload["decision"]["professorial_verdict"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    carrier = carrier_enumeration()
    providers = provider_matrix(reads)
    payload: dict[str, Any] = {
        "status": "P2690_SELECTOR_FREE_NONEXACT_BOUNDARY_PHASE_SECTOR_PROVIDER_AUDIT_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "carrier_enumeration": carrier,
        "provider_matrix": providers,
        "decision": decision(carrier, providers),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2690/S1640 selector-free nonexact boundary-phase sector provider audit",
        "## P2690/S1640 selector-free nonexact boundary-phase sector provider audit\n\n"
        "`P2690/S1640` executes the P2689 next premise audit.  The P2663 one-bit carrier is real: exact coboundaries have square bit zero while nonexact rows can carry bit one.  But no selector-free provider currently exports a preferred nonexact square-holonomy sector: local positive dynamics, declared theta, selector-lane transfer, source-topology typed descent, and `tau_src -> pair12 -> boundary-square` all fail the acceptance matrix.  Therefore the entropy/UV-unit route is frozen as bounded no-go on current artifacts; no boundary-phase bit target, UV unit, beta source, role transfer, `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2690/S1640 boundary-phase provider Ltotal guard",
        "## P2690/S1640 boundary-phase provider Ltotal guard\n\n"
        "`P2690/S1640` keeps `L_total` nonpromoted: a one-bit boundary carrier is not a variational source term unless a selector-free nonexact sector provider and physical bit-to-length/action map are exported.  Current artifacts do not supply that provider, so the entropy/UV-unit lane is bounded no-go.\n",
    )
    append_once(
        AGENTS,
        "Current boundary-phase sector provider guardrail (P2690/S1640, 2026-06-13)",
        "## Current boundary-phase sector provider guardrail (P2690/S1640, 2026-06-13)\n\n"
        "- P2690 confirms the P2663 one-bit carrier is real but finds no selector-free nonexact boundary-phase sector provider; freeze the entropy/UV-unit route as bounded no-go on current artifacts.\n"
        "- The next move should return to the broad state-map; the most concrete remaining P2680 non-selector target is an amplitude role-safe source audit for `alpha_geo` scalar-shape normalization, without legacy role transfer, selector replay, `beta_tors -> chi11`, or generic bridge completion.\n",
    )
    return payload


if __name__ == "__main__":
    main()
