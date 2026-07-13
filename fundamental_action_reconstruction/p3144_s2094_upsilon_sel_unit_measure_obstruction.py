"""P3144/S2094: Upsilon_sel unit-measure obstruction.

P3143 requested one genuinely new typed unit-measure object coupling the selector
local chart to an action measure without importing A_origin/A_lambda.  This
packet constructs the natural candidate Upsilon_sel^unit on X = Z12 x {±1} and
computes the finite symmetry obstruction: under the source-free diagonal
translation plus polarity-pairing action, the invariant unit measure is uniform
on the single 24-point orbit and cannot localize the selector branch.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3144_s2094_upsilon_sel_unit_measure_obstruction.json"
MD = GEN / "p3144_s2094_upsilon_sel_unit_measure_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3143": GEN / "p3143_s2093_vlift_weight_scale_source_audit.json",
    "P3142": GEN / "p3142_s2092_axiom_selector_field_variational_lift.json",
    "P2887": GEN / "p2887_s1837_c12_unit_measure_localized_density_9_over_5_no_go_audit.json",
    "P3139": GEN / "p3139_s2089_dhl_lane_no_new_frontier_reconciliation.json",
}

N = 12
LAMBDAS = (-1, 1)
SELECTED = (0, 1)
Point = tuple[int, int]
Action = Callable[[Point], Point]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def space() -> list[Point]:
    return [(r, lam) for r in range(N) for lam in LAMBDAS]


def translate(t: int) -> Action:
    return lambda p: ((p[0] + t) % N, p[1])


def invert(p: Point) -> Point:
    return ((-p[0]) % N, -p[1])


def flip_lambda(p: Point) -> Point:
    return (p[0], -p[1])


def generators() -> dict[str, Action]:
    return {"T_plus_1": translate(1), "Aut_inversion_with_polarity": invert, "lambda_pairing": flip_lambda}


def orbit(seed: Point, gens: dict[str, Action]) -> list[Point]:
    seen = {seed}
    frontier = [seed]
    while frontier:
        p = frontier.pop()
        for fn in gens.values():
            q = fn(p)
            if q not in seen:
                seen.add(q)
                frontier.append(q)
    return sorted(seen)


def all_orbits(gens: dict[str, Action]) -> list[list[Point]]:
    remaining = set(space())
    orbits: list[list[Point]] = []
    while remaining:
        seed = next(iter(remaining))
        orb = orbit(seed, gens)
        orbits.append(orb)
        remaining -= set(orb)
    return sorted(orbits, key=len, reverse=True)


def measure_uniform(points: list[Point]) -> dict[str, str]:
    mass = Fraction(1, len(points))
    return {f"{r},{lam}": str(mass) for r, lam in points}


def invariance_failures(measure: dict[Point, Fraction], gens: dict[str, Action]) -> list[dict[str, Any]]:
    failures = []
    for name, fn in gens.items():
        for p, mass in measure.items():
            q = fn(p)
            if measure[q] != mass:
                failures.append({"generator": name, "point": list(p), "image": list(q), "mass": str(mass), "image_mass": str(measure[q])})
    return failures


def candidate_measures(gens: dict[str, Action]) -> list[dict[str, Any]]:
    pts = space()
    uniform_mass = {p: Fraction(1, 24) for p in pts}
    delta_mass = {p: Fraction(0, 1) for p in pts}
    delta_mass[SELECTED] = Fraction(1, 1)
    polarity_mass = {p: Fraction(1, 12) if p[1] == 1 else Fraction(0, 1) for p in pts}
    candidates = [
        {"id": "U1_uniform_full_orbit", "measure": uniform_mass, "imports_selector_axiom": False, "unit_normalized": True, "localized_on_selected": False},
        {"id": "U2_delta_selected_branch", "measure": delta_mass, "imports_selector_axiom": True, "unit_normalized": True, "localized_on_selected": True},
        {"id": "U3_positive_polarity_orbit", "measure": polarity_mass, "imports_selector_axiom": True, "unit_normalized": True, "localized_on_selected": False},
    ]
    rows = []
    for cand in candidates:
        failures = invariance_failures(cand["measure"], gens)
        full_invariant = len(failures) == 0
        accepted = full_invariant and cand["unit_normalized"] and cand["localized_on_selected"] and not cand["imports_selector_axiom"]
        rows.append({
            "candidate_id": cand["id"],
            "unit_normalized": cand["unit_normalized"],
            "full_symmetry_invariant": full_invariant,
            "localized_on_selected_branch": cand["localized_on_selected"],
            "imports_selector_axiom": cand["imports_selector_axiom"],
            "invariance_failure_count": len(failures),
            "sample_failures": failures[:6],
            "accepted_as_upsilon_sel_unit_measure": accepted,
        })
    return rows


def build_payload() -> dict[str, Any]:
    inputs = {name: load_json(path) for name, path in INPUTS.items()}
    gens = generators()
    orbits = all_orbits(gens)
    uniform = measure_uniform(orbits[0])
    candidates = candidate_measures(gens)
    accepted = [c for c in candidates if c["accepted_as_upsilon_sel_unit_measure"]]
    theorem = {
        "name": "P3144_T1_upsilon_sel_invariant_unit_measure_no_local_selector",
        "statement": "On X=Z12 x {±1}, the source-free action generated by diagonal translation, inversion-with-polarity, and lambda pairing has one 24-point orbit.  Therefore any invariant unit measure is constant on all selector branches; the unique full-orbit probability measure assigns mass 1/24 to the axiom-selected branch and cannot localize it.  Delta/localized measures localize only by importing selector axioms and fail symmetry invariance.",
        "finite_counts": {
            "point_count": len(space()),
            "generator_count": len(gens),
            "orbit_count": len(orbits),
            "largest_orbit_size": len(orbits[0]),
            "candidate_measures_tested": len(candidates),
            "accepted_unit_selector_measures": len(accepted),
        },
    }
    return {
        "status": "P3144_UPSILON_SEL_UNIT_MEASURE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "repo_context": {name: inputs[name].get("status") for name in inputs},
        "constructed_object": {
            "name": "Upsilon_sel^unit source-free invariant unit measure candidate",
            "space": "Z12 x {±1}",
            "selected_branch": list(SELECTED),
            "classification": "unit_measure_candidate_and_symmetry_obstruction",
        },
        "symmetry_generators": list(gens.keys()),
        "orbit_rows": [{"orbit_index": idx, "size": len(orb), "points_sample": [list(p) for p in orb[:8]]} for idx, orb in enumerate(orbits)],
        "uniform_measure_sample": dict(list(uniform.items())[:8]),
        "candidate_measure_rows": candidates,
        "finite_theorem": theorem,
        "decision": {
            "bounded_result": "P3144 constructs the requested Upsilon_sel unit-measure object.  A source-free unit measure exists, but it is the uniform full-orbit measure and is selector-blind.  Measures that localize the selected branch import A_origin/A_lambda and fail the required symmetry invariance.  Thus no strict unit-bearing selector action is installed.",
            "axiomatic_route_recommendation": "The axiom route is worth keeping only as a separate conditional model-building branch for downstream consequences.  It is not worth presenting as the main strict-core route until a non-premise symmetry-breaking source or a non-uniform unit measure theorem is exported.",
            "negative_export_flags": {
                "nonpremise_upsilon_sel_unit_measure_exported": False,
                "unit_bearing_selector_action_exported": False,
                "selector_localization_without_axioms_exported": False,
                "strict_QW_2191_discharged": False,
                "strict_selector_closure_exported": False,
                "bridge_completion_exported": False,
                "legacy_role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "next_honest_step": "Pivot away from the axiom branch for strict-core work unless a new nonuniform invariant measure/source is supplied.  The next proof-grade strict move should be a new non-premise symmetry-breaking source candidate outside invariant unit measures; the axiom branch may continue only as explicitly labelled conditional phenomenology, not as ToE evidence.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = [
        "# P3144/S2094 Upsilon_sel unit-measure obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Space: `{payload['constructed_object']['space']}`",
        f"- Selected branch: `{payload['constructed_object']['selected_branch']}`",
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Finite counts",
    ]
    for key, value in th["finite_counts"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Candidate measures"])
    for row in payload["candidate_measure_rows"]:
        lines.append(f"- `{row['candidate_id']}`: unit `{row['unit_normalized']}`, invariant `{row['full_symmetry_invariant']}`, localized `{row['localized_on_selected_branch']}`, imports selector `{row['imports_selector_axiom']}`, accepted `{row['accepted_as_upsilon_sel_unit_measure']}`")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Axiomatic-route recommendation", payload["decision"]["axiomatic_route_recommendation"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3144/S2094 Upsilon_sel unit-measure obstruction", "## P3144/S2094 Upsilon_sel unit-measure obstruction\n\n`P3144/S2094` constructs the requested `Upsilon_sel^unit` candidate on `Z12 x {±1}`.  The source-free action generated by translation, inversion-with-polarity, and lambda pairing has one `24`-point orbit, so the invariant unit measure is uniform and selector-blind (`1/24` on the axiom-selected branch).  Localized delta/polarity measures import selector axioms and fail full symmetry invariance.  No strict unit-bearing selector action, `QW-2191` discharge, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3144/S2094 Upsilon_sel measure remains selector-blind", "## P3144/S2094 Upsilon_sel measure remains selector-blind\n\n`P3144/S2094` supplies the finite unit-measure obstruction requested after P3143: the only source-free invariant unit measure on the selector pair space is uniform over the full `24`-point orbit, while any localized selector measure imports the axiom branch.  Thus no physical action measure or selector-localized EOM term is installed.\n")
    append_once(AGENTS, "Current Upsilon_sel unit-measure obstruction guardrail (P3144/S2094, 2026-07-13)", "## Current Upsilon_sel unit-measure obstruction guardrail (P3144/S2094, 2026-07-13)\n\n- P3144 constructs the requested `Upsilon_sel^unit` candidate on `Z12 x {±1}` and proves the source-free invariant measure is uniform on one `24`-point orbit.\n- The uniform measure is unit-normalized but selector-blind; localized delta/polarity measures import `A_origin/A_lambda` and fail full symmetry invariance.\n- Do not promote invariant unit measure, localized axiom measures, or the P3140-P3144 axiom branch to strict `QW-2191` discharge, selector closure, unit-bearing action/EOM, bridge completion, role transfer, `L_total`, or ToE closure.\n- Recommendation: keep the axiom route only as explicitly labelled conditional phenomenology; for strict-core work, pivot to a genuinely new non-premise symmetry-breaking source outside invariant unit measures.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
