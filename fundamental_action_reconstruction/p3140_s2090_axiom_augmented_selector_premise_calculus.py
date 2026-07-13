"""P3140/S2090: axiom-augmented selector premise calculus.

User question: what if symmetry breaking / selector is accepted axiomatically?
This packet constructs the finite non-strict axiom calculus for the D_HL joint
origin/polarity target.  It explicitly separates conditional axiom closure from
strict-source closure, as required by the repository guardrails.
"""

from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3140_s2090_axiom_augmented_selector_premise_calculus.json"
MD = GEN / "p3140_s2090_axiom_augmented_selector_premise_calculus.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "D1": GEN / "d1_selector_axiom_necessity_current_best_supported_conclusion_summary.json",
    "P1521": ROOT / "P1521_S471_MINIMAL_STRICT_OSPLIT_IO_CONTRACT_PACKET_PL.md",
    "P2164": GEN / "p2164_s1114_strict_qw2191_blocker_lane_baselined_entry_packet.json",
    "P3139": GEN / "p3139_s2089_dhl_lane_no_new_frontier_reconciliation.json",
}

N = 12
LAMBDAS = [-1, 1]
UNITS = [1, 5, 7, 11]

AXIOMS = {
    "A_origin": "postulate one absolute support representative r0 in Z12",
    "A_lambda": "postulate one polarity lambda0 in {+1,-1}",
    "A_coupling": "postulate that the selected pair (r0,lambda0) couples to the P3134 D_HL formula as the selector channel",
}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load_json(path: Path) -> dict[str, Any]:
    if path.exists() and path.suffix == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    return {}


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def pair_space() -> list[tuple[int, int]]:
    return [(r, lam) for r in range(N) for lam in LAMBDAS]


def transform(pair: tuple[int, int], t: int = 0, unit: int = 1, flip_lambda: bool = False) -> tuple[int, int]:
    r, lam = pair
    return ((unit * r + t) % N, -lam if flip_lambda else lam)


def orbit(seed: tuple[int, int], allow_translation: bool = True, allow_aut: bool = True, allow_lambda_flip: bool = True) -> set[tuple[int, int]]:
    translations = range(N) if allow_translation else [0]
    units = UNITS if allow_aut else [1]
    flips = [False, True] if allow_lambda_flip else [False]
    return {transform(seed, t, u, f) for t in translations for u in units for f in flips}


def survivors(axioms: set[str], r0: int = 0, lambda0: int = 1) -> list[tuple[int, int]]:
    out = []
    for r, lam in pair_space():
        if "A_origin" in axioms and r != r0:
            continue
        if "A_lambda" in axioms and lam != lambda0:
            continue
        out.append((r, lam))
    return out


def axiom_row(names: tuple[str, ...]) -> dict[str, Any]:
    ax = set(names)
    kept = survivors(ax)
    unique_pair = len(kept) == 1
    has_coupling = "A_coupling" in ax
    conditional_selector = unique_pair and has_coupling
    return {
        "axioms": list(names),
        "axiom_count": len(names),
        "surviving_pairs": len(kept),
        "selects_unique_origin": len({r for r, _ in kept}) == 1,
        "selects_unique_lambda": len({lam for _, lam in kept}) == 1,
        "selects_unique_pair": unique_pair,
        "has_DHL_coupling_axiom": has_coupling,
        "conditional_non_strict_J_DHL_exported": conditional_selector,
        "strict_source_exported": False,
        "QW_2191_discharged_strictly": False,
        "classification": "non_strict_axiom_augmented_selector" if conditional_selector else "incomplete_axiom_package",
    }


def build_payload() -> dict[str, Any]:
    inputs = {name: load_json(path) for name, path in INPUTS.items()}
    all_rows = []
    names = list(AXIOMS)
    for k in range(0, len(names) + 1):
        for combo in combinations(names, k):
            all_rows.append(axiom_row(combo))

    no_axiom_orbit = orbit((0, 1))
    minimal_closing = [row for row in all_rows if row["conditional_non_strict_J_DHL_exported"] and row["axiom_count"] == 3]
    obligation_rows = [
        {"obligation": "absolute_origin", "strict_current": False, "axiom_needed": "A_origin", "non_strict_if_assumed": True},
        {"obligation": "unpaired_lambda", "strict_current": False, "axiom_needed": "A_lambda", "non_strict_if_assumed": True},
        {"obligation": "DHL_selector_coupling", "strict_current": False, "axiom_needed": "A_coupling", "non_strict_if_assumed": True},
        {"obligation": "strict_source_provenance", "strict_current": False, "axiom_needed": None, "non_strict_if_assumed": False},
        {"obligation": "variational_unit_Ltotal_bridge", "strict_current": False, "axiom_needed": "downstream theorem even after selector axioms", "non_strict_if_assumed": False},
    ]
    theorem = {
        "name": "P3140_T1_axiom_augmented_DHL_selector_calculus",
        "statement": "On the finite pair space Z12 x {±1}, no current strict artifact selects a unique joint pair.  Adding A_origin and A_lambda reduces the pair space to one stipulated pair, and adding A_coupling turns that pair into a conditional D_HL selector channel.  This is a valid non-strict axiom-augmented closure, but it is not a strict-source theorem and does not discharge QW-2191 in the strict core.",
        "finite_counts": {
            "pair_space_size": len(pair_space()),
            "no_axiom_orbit_size_under_translation_aut_lambda_flip": len(no_axiom_orbit),
            "axiom_packages_tested": len(all_rows),
            "minimal_non_strict_closing_packages": len(minimal_closing),
        },
    }
    return {
        "status": "P3140_AXIOM_AUGMENTED_SELECTOR_PREMISE_CALCULUS_NON_STRICT",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "repo_context": {
            "D1_project_conclusion": inputs["D1"].get("project_conclusion"),
            "P2164_selected_lane": inputs["P2164"].get("strict_qw2191_blocker_lane_baselined_entry_packet", {}).get("selected_next_blocker_lane"),
            "P3139_status": inputs["P3139"].get("status"),
        },
        "constructed_object": {
            "name": "S_sel^ax(D_HL) axiom-augmented selector premise calculus",
            "axioms": AXIOMS,
            "selected_pair_under_full_axiom_package": {"r0": 0, "lambda0": 1, "note": "0,+1 are conventionally named representatives inside the axiom package; strict provenance is not claimed"},
        },
        "finite_theorem": theorem,
        "axiom_package_rows": all_rows,
        "obligation_rows": obligation_rows,
        "decision": {
            "bounded_result": "Assuming axioms can honestly break the selector symmetry, but only as a labelled non-strict extension.  The minimal D_HL package is A_origin + A_lambda + A_coupling: it selects one stipulated pair and couples it to D_HL.  It does not prove strict provenance, does not erase QW-2191 from the strict core, and does not by itself supply variational/unit/L_total or ToE closure.",
            "positive_scoped_flags": {
                "non_strict_axiom_calculus_exported": True,
                "minimal_axiom_package_identified": True,
                "conditional_non_strict_J_DHL_exported_under_full_axioms": bool(minimal_closing),
            },
            "negative_export_flags": {
                "strict_selector_source_exported": False,
                "QW_2191_discharged_strictly": False,
                "D_HL_strict_source_exported": False,
                "Zeta_OS_exported_strictly": False,
                "Gamma_SO_exported_strictly": False,
                "bridge_completion_exported": False,
                "legacy_role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "next_honest_step": "If the project accepts an axiom route, formalize a separate non-strict axiom-augmented branch with A_origin, A_lambda, and A_coupling as explicit assumptions, then audit downstream consequences without calling them strict.  If the project stays strict-core, do not use the axiom package as evidence; pivot to one genuinely new strict source object that can derive one of these axioms rather than assume it.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = [
        "# P3140/S2090 Axiom-augmented selector premise calculus",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        "- Axioms:",
    ]
    for name, desc in payload["constructed_object"]["axioms"].items():
        lines.append(f"  - `{name}`: {desc}")
    lines.extend([
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Finite counts",
    ])
    for key, value in th["finite_counts"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Obligations"])
    for row in payload["obligation_rows"]:
        lines.append(f"- `{row['obligation']}`: strict_current=`{row['strict_current']}`, axiom_needed=`{row['axiom_needed']}`, non_strict_if_assumed=`{row['non_strict_if_assumed']}`")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3140/S2090 axiom-augmented selector premise calculus", "## P3140/S2090 axiom-augmented selector premise calculus\n\n`P3140/S2090` answers the axiom route explicitly.  On `Z12 x {±1}`, the non-strict axiom package `A_origin + A_lambda + A_coupling` selects one stipulated `(r0,lambda0)` pair and couples it to the P3134 `D_HL` channel.  This is a valid axiom-augmented selector branch, but it is explicitly non-strict: no strict selector source, strict `QW-2191` discharge, `D_HL` source, `Zeta_OS`, `Gamma_SO`, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3140/S2090 axiom branch is non-strict and non-variational by itself", "## P3140/S2090 axiom branch is non-strict and non-variational by itself\n\n`P3140/S2090` permits a labelled non-strict branch in which selector axioms choose and couple a `D_HL` pair.  The branch cannot be used as strict variational provenance: unit-bearing action, Hamiltonian normalization, spacetime EOM, `L_total`, bridge-completion theorem, role-transfer theorem, and ToE closure remain downstream and unexported.\n")
    append_once(AGENTS, "Current axiom-augmented selector premise calculus guardrail (P3140/S2090, 2026-07-13)", "## Current axiom-augmented selector premise calculus guardrail (P3140/S2090, 2026-07-13)\n\n- P3140 constructs the explicit non-strict axiom route for the `D_HL` selector lane: `A_origin` chooses an absolute support representative, `A_lambda` chooses polarity, and `A_coupling` couples the chosen pair to `D_HL`.\n- The full package conditionally selects one `(r0,lambda0)` pair on `Z12 x {±1}`, but this is an axiom-augmented branch, not strict-source closure.\n- Do not promote the axiom package to strict `QW-2191` discharge, strict selector closure, `D_HL`/`Zeta_OS`/`Gamma_SO` sourcehood, bridge completion, role transfer, `L_total`, or ToE closure.\n- Next honest move: either maintain a separate non-strict axiom branch and audit downstream consequences under explicit assumptions, or, for strict-core work, construct a genuinely new source theorem deriving one of `A_origin`, `A_lambda`, or `A_coupling` rather than assuming it.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
