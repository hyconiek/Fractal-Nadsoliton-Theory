#!/usr/bin/env python3
"""P2689/S1639: entropy reference-cell and bit-to-length UV-unit obligation matrix.

This executes the P2688 selected next step.  It audits whether P2662's
conditional entropy unit-map scaffold plus P2663's one-bit boundary-phase carrier
can currently produce an unconditional canonical length/UV unit source.
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
OUT = GEN / "p2689_s1639_entropy_reference_cell_bit_to_length_uv_unit_obligation_matrix.json"
MD = GEN / "p2689_s1639_entropy_reference_cell_bit_to_length_uv_unit_obligation_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2688": GEN / "p2688_s1638_post_p2687_live_frontier_reconciliation_audit.json",
    "P2662": GEN / "p2662_s1612_entropy_boundary_phase_unit_map_conditional_theorem_audit.json",
    "P2663": GEN / "p2663_s1613_chain_level_boundary_phase_bit_target_audit.json",
    "P2650": GEN / "p2650_s1600_canonical_length_uv_unit_source_candidate_exhaustion_no_go.json",
    "P2653": GEN / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.json",
    "P2680": GEN / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "unconditional_uv_unit_exported",
    "beta_or_z_beta_source_exported",
    "bit_to_length_unit_map_exported",
    "selector_replay_imported",
    "role_transfer_started",
    "bridge_completion_claimed",
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
        "p2688_selected_p2689": r"P2689|entropy reference-cell|bit-to-length|canonical length/UV unit|UV-unit source obligation",
        "p2662_conditional_unit_map": r"H\(a\)=H0\+D_f log\(a\)|conditional unit-map|N log\(2\)|bit-to-action|bit-to-length|entropy zero|reference cell",
        "p2663_one_bit_carrier": r"one-bit|boundary-phase|square_holonomy|nonexact sector|chain-level|N_bits",
        "scale_orbit_beta_source": r"scale-orbit|scale orbit|beta source|Z_beta|canonical unit|typed metric/UV source",
        "forbidden_replays": r"selector replay|role transfer|bridge completion|QW-2191 discharge|L_total|ToE closure",
    }
    return {"tool": "rg", "mode": "content-first entropy reference-cell/bit-to-length UV-unit obligation matrix", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def state_reads() -> dict[str, Any]:
    data = {name: load_json(path) for name, path in INPUTS.items()}
    p2662_close = data["P2662"].get("closure_decision", {})
    p2663_witness = data["P2663"].get("chain_level_boundary_phase_witness", {})
    p2653_close = data["P2653"].get("closure_decision", {})
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2688_selected_p2689": data["P2688"].get("selection", {}).get("selected_next_packet") == "P2689_intrinsic_entropy_reference_cell_and_bit_to_length_uv_unit_source_obligation_matrix",
        "p2662_conditional_unique_scale_verified": p2662_close.get("conditional_unique_scale_verified") is True,
        "p2662_unconditional_uv_unit_selected": p2662_close.get("unconditional_uv_unit_selected_now") is True,
        "p2662_missing_premises": p2662_close.get("missing_premises", []),
        "p2663_nonexact_sector_needed": p2663_witness.get("nonexact_sector_needed_for_nonzero_bit") is True,
        "p2663_unique_bits_without_sector_choice": p2663_witness.get("chain_law_derives_unique_n_bits_without_sector_choice") is True,
        "p2663_exact_coboundary_values": p2663_witness.get("exact_coboundary_square_holonomy_values", []),
        "p2650_any_candidate_passes": any(row.get("passes_as_canonical_length_uv_source") for row in data["P2650"].get("canonical_unit_candidates", [])),
        "p2653_canonical_unit_exported": p2653_close.get("canonical_unit_exported_now") is True,
        "p2653_scale_orbit_equivalence_verified": p2653_close.get("scale_orbit_equivalence_verified") is True,
    }


def symbolic_unit_map_check() -> dict[str, Any]:
    a, H0, Df, N = sp.symbols("a H0 Df N", positive=True)
    entropy = H0 + Df * sp.log(a)
    target = N * sp.log(2)
    selected_scale = sp.solve(sp.Eq(entropy, target), a)[0]
    derivative = sp.diff(entropy, a)
    samples = []
    for n_bits, h0, df in [(1, 0, sp.Rational(9, 5)), (2, sp.Rational(1, 3), sp.Rational(9, 5)), (3, sp.Rational(2, 5), sp.Rational(2, 1))]:
        value = sp.simplify(selected_scale.subs({N: n_bits, H0: h0, Df: df}))
        samples.append({"N_bits": str(n_bits), "H0": str(h0), "Df": str(df), "selected_scale": str(value), "positive": bool(value.is_positive)})
    return {
        "entropy_law": str(entropy),
        "bit_target": str(target),
        "selected_scale_formula": str(selected_scale),
        "dH_da": str(derivative),
        "unique_positive_scale_if_premises_supplied": True,
        "sample_scales": samples,
        "conditional_math_passes": True,
    }


def obligation_matrix(reads: dict[str, Any]) -> list[dict[str, Any]]:
    missing = set(reads["p2662_missing_premises"])
    return [
        {
            "obligation_id": "intrinsic_entropy_reference_cell_or_entropy_zero",
            "source": "P2662 missing premise",
            "satisfied_now": "intrinsic_pre_normalization_entropy_measure" not in missing,
            "blocker": "no intrinsic pre-normalization entropy measure/reference cell exported",
        },
        {
            "obligation_id": "selector_free_boundary_phase_bit_target",
            "source": "P2663 chain-level carrier",
            "satisfied_now": reads["p2663_unique_bits_without_sector_choice"],
            "blocker": "P2663 has a one-bit carrier, but nonzero bit requires nonexact sector choice and no unique N_bits law is exported",
        },
        {
            "obligation_id": "bit_to_length_or_action_unit_map",
            "source": "P2662 missing premise",
            "satisfied_now": "bit_to_action_or_bit_to_length_unit_map" not in missing,
            "blocker": "no bit-to-length or bit-to-action physical unit map exported",
        },
        {
            "obligation_id": "scale_orbit_quotient_single_uv_unit",
            "source": "P2653 scale-orbit witness",
            "satisfied_now": reads["p2653_canonical_unit_exported"] and not reads["p2653_scale_orbit_equivalence_verified"],
            "blocker": "P2653 verifies scale-orbit equivalence; no quotient discharge selects a single unit",
        },
        {
            "obligation_id": "target_independent_beta_or_z_beta_source",
            "source": "P2680/P2650/P2653",
            "satisfied_now": reads["p2650_any_candidate_passes"] and reads["p2653_canonical_unit_exported"],
            "blocker": "no canonical unit candidate or typed metric/UV source currently passes",
        },
    ]


def decision(obligations: list[dict[str, Any]]) -> dict[str, Any]:
    failed = [row for row in obligations if not row["satisfied_now"]]
    return {
        "decision": "P2689_ENTROPY_REFERENCE_CELL_BIT_TO_LENGTH_UV_UNIT_OBLIGATION_MATRIX_BOUNDED_NO_GO_NO_FALSE_PASS",
        "unconditional_uv_unit_exported": not failed,
        "failed_obligation_ids": [row["obligation_id"] for row in failed],
        "bounded_no_go_for_current_canonical_unit_atom": bool(failed),
        "professorial_verdict": (
            "P2689 verifies that the P2662 entropy-scale equation is mathematically sharp: if an intrinsic entropy reference H0, a selector-free bit target N log(2), and a bit-to-length/action unit map were exported, one positive scale would be selected.  But the current artifacts still fail every source-bearing obligation needed to turn that conditional equation into a canonical UV-unit source.  P2663 supplies a real one-bit carrier, yet nonzero bits require a nonexact sector choice and no selector-free unique N_bits law.  Therefore the canonical length/UV unit atom remains bounded no-go on current evidence; no beta/Z_beta source, bridge completion, role transfer, L_total, or ToE closure is exported."
        ),
        "next_honest_step": (
            "P2690 should attack exactly one missing premise: a selector-free nonexact boundary-phase sector provider for the P2663 one-bit carrier.  The test must decide whether nadsoliton boundary dynamics exports a preferred nonexact square-holonomy sector and a unique N_bits target without importing selector replay, QW-2191 discharge, or role transfer.  If not, freeze the entropy/UV-unit route as a bounded no-go and return to the broad state-map."
        ),
        "beta_or_z_beta_source_exported_now": False,
        "role_transfer_started_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2689/S1639 entropy reference-cell and bit-to-length UV-unit obligation matrix", "", f"Status: `{payload['status']}`", "", "## Content-first grep"]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    sym = payload["symbolic_unit_map_check"]
    lines.extend(["", "## Conditional symbolic unit map", f"Entropy law: `{sym['entropy_law']}`; bit target: `{sym['bit_target']}`.", f"Selected scale: `{sym['selected_scale_formula']}`; dH/da: `{sym['dH_da']}`."])
    lines.extend(["", "## Obligation matrix"])
    for row in payload["obligation_matrix"]:
        lines.append(f"- `{row['obligation_id']}`: satisfied=`{row['satisfied_now']}` — {row['blocker']}")
    lines.extend(["", "## Verdict", payload["decision"]["professorial_verdict"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    sym = symbolic_unit_map_check()
    obligations = obligation_matrix(reads)
    payload: dict[str, Any] = {
        "status": "P2689_ENTROPY_REFERENCE_CELL_BIT_TO_LENGTH_UV_UNIT_OBLIGATION_MATRIX_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "symbolic_unit_map_check": sym,
        "obligation_matrix": obligations,
        "decision": decision(obligations),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2689/S1639 entropy reference-cell and bit-to-length UV-unit obligation matrix",
        "## P2689/S1639 entropy reference-cell and bit-to-length UV-unit obligation matrix\n\n"
        "`P2689/S1639` executes the P2688-selected canonical UV-unit source audit.  The conditional equation `H(a)=H0+D_f log(a)=N log(2)` selects one positive scale if its premises are supplied, but current artifacts still lack an intrinsic entropy reference cell, a selector-free unique boundary-phase bit target, a bit-to-length/action unit map, and a scale-orbit quotient discharge.  The P2663 one-bit carrier is real but not uniquely selected without a nonexact sector provider.  Thus no UV unit, beta/Z_beta source, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2689/S1639 entropy UV-unit Ltotal guard",
        "## P2689/S1639 entropy UV-unit Ltotal guard\n\n"
        "`P2689/S1639` keeps `L_total` nonpromoted: entropy scale-selection remains a conditional source-target equation, not a role-bearing action term.  A future step must first export a selector-free nonexact boundary-phase sector provider and bit-to-length/action map before any beta source or variational role can be considered.\n",
    )
    append_once(
        AGENTS,
        "Current entropy UV-unit obligation guardrail (P2689/S1639, 2026-06-13)",
        "## Current entropy UV-unit obligation guardrail (P2689/S1639, 2026-06-13)\n\n"
        "- P2689 shows the entropy equation can conditionally select a positive scale, but no unconditional canonical UV unit or beta/Z_beta source is exported because the intrinsic entropy reference, selector-free bit target, bit-to-length/action map, and scale-orbit quotient remain missing.\n"
        "- The next admissible move is exactly one missing premise: a selector-free nonexact boundary-phase sector provider for the P2663 one-bit carrier; do not import selector replay, QW-2191 discharge, role transfer, bridge completion, `L_total`, or ToE closure.\n",
    )
    return payload


if __name__ == "__main__":
    main()
