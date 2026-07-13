"""P3146/S2096: axiomatic unit/action-measure bridge candidate.

P3145 left physical units/action measure as the most bottlenecked reverse-layout
row.  This packet takes the next honest axiomatic step: construct the smallest
conditional bridge object that can make a strict finite receiver look like a
unit-bearing action density, then compute which assumptions are mathematically
necessary.  It is an axiom-branch object only, not a strict source theorem.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3146_s2096_axiom_unit_action_measure_bridge.json"
MD = GEN / "p3146_s2096_axiom_unit_action_measure_bridge.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3145": GEN / "p3145_s2095_strict_kernel_reverse_sm_gr_layout.json",
    "P3144": GEN / "p3144_s2094_upsilon_sel_unit_measure_obstruction.json",
    "P3143": GEN / "p3143_s2093_vlift_weight_scale_source_audit.json",
    "P3104": GEN / "p3104_s2054_spectral_triple_geometry_interface_obstruction_audit.json",
    "P3094": GEN / "p3094_s2044_stress_energy_metric_response_obstruction_audit.json",
}

# Dimension basis: action H, length L, time T.  The receiver R is dimensionless.
AXIOMS = {
    "A_cell": {
        "meaning": "postulate a physical cell length ell_* converting Z12 receiver steps to length",
        "dimension_vector": [0, 1, 0],
        "strict_sourced": False,
        "imports_selector_axioms": False,
    },
    "A_clock": {
        "meaning": "postulate a clock scale tau_* converting gate/evolution labels to time",
        "dimension_vector": [0, 0, 1],
        "strict_sourced": False,
        "imports_selector_axioms": False,
    },
    "A_action": {
        "meaning": "postulate an action quantum hbar_* converting dimensionless finite sums to action",
        "dimension_vector": [1, 0, 0],
        "strict_sourced": False,
        "imports_selector_axioms": False,
    },
}

TARGETS = {
    "length_measure": [0, 1, 0],
    "time_measure": [0, 0, 1],
    "action_weight": [1, 0, 0],
    "lagrangian_density_1d": [1, -1, -1],  # action per length per time in the minimal 1D chart
}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def rank_int(rows: list[list[int]]) -> int:
    # Small exact rational Gaussian elimination using fractions.
    from fractions import Fraction

    mat = [[Fraction(x) for x in row] for row in rows if any(row)]
    if not mat:
        return 0
    m, n = len(mat), len(mat[0])
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, m) if mat[i][c]), None)
        if pivot is None:
            continue
        mat[r], mat[pivot] = mat[pivot], mat[r]
        div = mat[r][c]
        mat[r] = [v / div for v in mat[r]]
        for i in range(m):
            if i != r and mat[i][c]:
                factor = mat[i][c]
                mat[i] = [a - factor * b for a, b in zip(mat[i], mat[r])]
        r += 1
        if r == m:
            break
    return r


def in_span(vectors: list[list[int]], target: list[int]) -> bool:
    return rank_int(vectors) == rank_int(vectors + [target])


def audit_axiom_sets() -> list[dict[str, Any]]:
    names = list(AXIOMS)
    rows: list[dict[str, Any]] = []
    for size in range(len(names) + 1):
        for combo in itertools.combinations(names, size):
            vectors = [AXIOMS[name]["dimension_vector"] for name in combo]
            span = {target: in_span(vectors, vec) for target, vec in TARGETS.items()}
            rows.append({
                "axiom_set": list(combo),
                "rank_HTL": rank_int(vectors),
                "covered_targets": span,
                "all_unit_targets_covered": all(span.values()),
                "strict_source_exported": False,
                "unit_bearing_action_measure_available_conditionally": all(span.values()),
                "blocked_as_strict_by": [
                    "cell length, clock scale, and action quantum are postulates, not strict-source exports",
                    "no nonproxy EH/ELg/gauge EOM coupling is derived from the finite receiver",
                    "does not resolve selector/orientation source obligations",
                ] if all(span.values()) else ["dimension rank is insufficient for a unit-bearing action density"],
            })
    return rows


def build_payload() -> dict[str, Any]:
    rows = audit_axiom_sets()
    full = [row for row in rows if row["all_unit_targets_covered"]]
    minimal_full = [row for row in full if len(row["axiom_set"]) == min(len(r["axiom_set"]) for r in full)] if full else []
    theorem = {
        "name": "P3146_T1_minimal_axiom_unit_action_measure_bridge",
        "statement": "In the dimension basis (H,L,T), a unit-bearing 1D selector/action-measure bridge from dimensionless strict finite receivers requires independent length, time, and action scale inputs.  The finite audit over all 8 subsets of {A_cell,A_clock,A_action} shows that only the full triple spans length, time, action, and the minimal Lagrangian-density dimension H L^-1 T^-1.  Thus the axiomatic bridge is internally coherent but non-strict: it installs units only by three explicit postulates and exports no non-premise source law.",
        "finite_counts": {
            "axiom_subsets_audited": len(rows),
            "dimension_basis_size": 3,
            "fully_covering_subsets": len(full),
            "minimal_fully_covering_subsets": len(minimal_full),
            "strict_source_subsets": 0,
            "conditional_unit_action_measure_subsets": len(full),
        },
    }
    return {
        "status": "P3146_AXIOM_UNIT_ACTION_MEASURE_BRIDGE_CONDITIONAL_NON_STRICT",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "constructed_object": {
            "name": "Lambda_unit^ax = (A_cell, A_clock, A_action) conditional unit/action-measure bridge",
            "classification": "axiom_augmented_dimensionful_measure_bridge_for_strict_finite_receivers",
            "formula": "S_ax[R]=hbar_* sum_{x in Z12 x {±1}} mu_unit(x) R(x), with d\u2113=ell_*, dt=tau_*, and density dimension hbar_*/(ell_* tau_*)",
            "selector_safety": "does not import A_origin/A_lambda; consequently it remains selector-blind and cannot localize a branch",
        },
        "axioms": AXIOMS,
        "targets": TARGETS,
        "audit_rows": rows,
        "finite_theorem": theorem,
        "decision": {
            "bounded_result": "The missing physical unit/action-measure object can be built as a clean conditional axiom bridge.  Computation shows the bridge needs exactly the full triple A_cell + A_clock + A_action to span the minimal unit targets; no proper subset suffices.",
            "why_not_strict": "All three scales are imported postulates.  The construction does not source ell_*, tau_*, or hbar_* from nadsoliton data, does not derive nonproxy metric/gauge EOM, and does not solve selector/orientation.",
            "next_honest_step": "Pick exactly one of the three postulates, preferably A_action because it is the least duplicated by prior length/entropy lanes, and audit whether any current or newly proposed strict object can source a nonzero action quantum hbar_* without legacy role transfer, selector replay, or Upsilon invariant-measure replay.",
            "negative_export_flags": {
                "strict_unit_source_exported": False,
                "unit_bearing_L_total_exported": False,
                "nonproxy_EH_ELg_EOM_exported": False,
                "strict_selector_closure_exported": False,
                "SM_GR_reduction_exported": False,
                "ToE_closure_exported": False,
            },
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = [
        "# P3146/S2096 Axiom unit/action-measure bridge",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Classification: `{payload['constructed_object']['classification']}`",
        f"- Formula: `{payload['constructed_object']['formula']}`",
        f"- Selector safety: {payload['constructed_object']['selector_safety']}",
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Finite counts",
    ]
    for key, value in th["finite_counts"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Axiom-subset audit"])
    for row in payload["audit_rows"]:
        lines.append(f"- `{row['axiom_set']}`: rank `{row['rank_HTL']}`, all targets covered `{row['all_unit_targets_covered']}`, strict source `{row['strict_source_exported']}`")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Why this is not strict", payload["decision"]["why_not_strict"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3146/S2096 axiom unit/action-measure bridge", "## P3146/S2096 axiom unit/action-measure bridge\n\n`P3146/S2096` constructs `Lambda_unit^ax=(A_cell,A_clock,A_action)`, a conditional dimensionful bridge from dimensionless strict finite receivers to an action-measure expression `S_ax[R]=hbar_* sum mu_unit R`.  The exact finite dimension audit over all `8` axiom subsets shows that only the full triple spans length, time, action, and the minimal density dimension `H L^-1 T^-1`; no strict unit source, `L_total`, SM/GR reduction, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3146/S2096 Lambda_unit^ax remains axiomatic", "## P3146/S2096 Lambda_unit^ax remains axiomatic\n\n`P3146/S2096` supplies a clean conditional unit/action-measure scaffold for axiom-branch model building, but it imports `ell_*`, `tau_*`, and `hbar_*` as postulates.  It therefore does not derive a unit-bearing physical Lagrangian, nonproxy Einstein/gauge EOM, or strict source of action units from the nadsoliton.\n")
    append_once(AGENTS, "Current axiom unit/action-measure bridge guardrail (P3146/S2096, 2026-07-13)", "## Current axiom unit/action-measure bridge guardrail (P3146/S2096, 2026-07-13)\n\n- P3146 constructs `Lambda_unit^ax=(A_cell,A_clock,A_action)`, a conditional bridge assigning length, time, and action scales to dimensionless strict finite receivers.\n- The exact audit over all `8` subsets proves that only the full triple spans length, time, action, and minimal `H L^-1 T^-1` density; no proper subset gives a unit-bearing action measure.\n- This is axiomatic and non-strict: `ell_*`, `tau_*`, and `hbar_*` are imported postulates, no nonproxy EOM/gauge coupling is derived, and selector/orientation remains unresolved.\n- Do not promote `Lambda_unit^ax` to strict unit source, `L_total`, SM/GR reduction, bridge completion, role transfer, or ToE closure.  Next honest move: audit exactly one postulate, preferably `A_action`, for a non-premise strict source of `hbar_*`.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
