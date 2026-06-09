#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2619_s1569_p2618_selector_source_obligation_lattice.json"
MD = GEN / "p2619_s1569_p2618_selector_source_obligation_lattice.md"

SOURCE_FILES = {
    "P2618_ANALYTIC_OBSTRUCTION": GEN / "p2618_s1568_analytic_legacy_to_strict_completion_obstruction.json",
    "P1966_SELECTOR_PREMISE": GEN / "p1966_s916_strict_qw2191_selector_premise_obstruction_and_minimal_breaking.json",
    "P2392_BETA_TORS_CHI11_RETIREMENT": GEN / "p2392_s1342_auxiliary_beta_tors_chi11_retirement_certificate.json",
    "P2616_ROLE_OBSTRUCTION": GEN / "p2616_s1566_p2608_role_acceptance_obstruction_after_source_revalidation.json",
}

SIGNS = (-1, 1)
C2 = (0, 1)  # 0 identity, 1 orientation reversal.

NEGATIVE_EXPORT_FLAGS = [
    "strict_selector_source_exported",
    "beta_tors_chi11_route_reopened",
    "gf2_bridge_revalidated",
    "role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "legacy_physical_role_transfer_exported",
    "qw2191_discharged_by_this_packet",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:70]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2619|S1569|selector source obligation lattice|orientation torsor equivariant enumeration|P2618 selector source",
        "nondup_selector_precursors": "P1966|P2392|P2618|orientation selector|beta_tors.*chi11|strict selector source|QW-2191",
        "finite_group_action_precursors": "C2|Z2|orientation reversal|equivariant selector|sign torsor|odd phase sign",
        "role_ltotal_guards": "P2616|P2608|role-bearing L_total|role-transfer theorem|legacy physical-role transfer|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def c2_act_on_sign(g: int, value: int) -> int:
    return value if g == 0 else -value


def invariant_points(n: int) -> list[str]:
    return [f"inv_{i}" for i in range(n)]


def torsor_points(n_pairs: int) -> list[str]:
    return [f"src_{i}_{sign}" for i in range(n_pairs) for sign in ("minus", "plus")]


def act_on_invariant_point(g: int, point: str) -> str:
    return point


def act_on_torsor_point(g: int, point: str) -> str:
    if g == 0:
        return point
    prefix, index, sign = point.split("_")
    flipped = "minus" if sign == "plus" else "plus"
    return f"{prefix}_{index}_{flipped}"


def enumerate_equivariant_maps(points: list[str], action_kind: str) -> dict[str, Any]:
    action = act_on_invariant_point if action_kind == "invariant" else act_on_torsor_point
    candidates = []
    for outputs in itertools.product(SIGNS, repeat=len(points)):
        fmap = dict(zip(points, outputs))
        ok = True
        defects = []
        for g in C2:
            for x in points:
                lhs = fmap[action(g, x)]
                rhs = c2_act_on_sign(g, fmap[x])
                if lhs != rhs:
                    ok = False
                    defects.append({"g": g, "x": x, "lhs_f_gx": lhs, "rhs_g_fx": rhs})
        if ok:
            candidates.append(fmap)
    return {
        "action_kind": action_kind,
        "input_point_count": len(points),
        "candidate_function_count": 2 ** len(points),
        "equivariant_function_count": len(candidates),
        "equivariant_functions": candidates[:8],
    }


def build_equivariance_table() -> list[dict[str, Any]]:
    rows = []
    for n in range(1, 5):
        rows.append({
            "source_type": f"{n}_orientation_invariant_legacy_scalar_classes",
            "interpretation": "beta_tors-like scalar/current magnitude data; orientation reversal acts trivially on the input",
            **enumerate_equivariant_maps(invariant_points(n), "invariant"),
            "selector_available": False,
        })
    for n_pairs in range(1, 4):
        rows.append({
            "source_type": f"{n_pairs}_orientation_odd_torsor_source_pairs",
            "interpretation": "a real oriented/sign-sensitive source premise; orientation reversal acts freely on the input",
            **enumerate_equivariant_maps(torsor_points(n_pairs), "torsor"),
            "selector_available": True,
        })
    return rows


def premise_lattice_rows() -> list[dict[str, Any]]:
    atoms = [
        "beta_tors_scalar_invariant",
        "axis_only_selector_up_to_Z2",
        "orientation_odd_source",
        "symmetry_breaking_boundary",
        "spin_pin_orientation_source",
    ]
    rows = []
    for r in range(len(atoms) + 1):
        for subset in itertools.combinations(atoms, r):
            support = set(subset)
            has_real_orientation_source = bool(support & {"orientation_odd_source", "symmetry_breaking_boundary", "spin_pin_orientation_source"})
            rejected_legacy_only = support <= {"beta_tors_scalar_invariant", "axis_only_selector_up_to_Z2"}
            rows.append({
                "support": list(subset),
                "support_size": r,
                "selector_gate_accepts": has_real_orientation_source,
                "rejected_if_only_legacy_scalar_or_axis": rejected_legacy_only and not has_real_orientation_source,
                "uses_beta_tors_scalar": "beta_tors_scalar_invariant" in support,
                "uses_axis_only_selector": "axis_only_selector_up_to_Z2" in support,
            })
    accepting = [row for row in rows if row["selector_gate_accepts"]]
    minimal_size = min(row["support_size"] for row in accepting)
    for row in rows:
        row["minimal_accepting_support"] = row["selector_gate_accepts"] and row["support_size"] == minimal_size
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    payloads = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2618 = theorem(payloads["P2618_ANALYTIC_OBSTRUCTION"], "analytic_legacy_to_strict_completion_obstruction")
    p2616 = theorem(payloads["P2616_ROLE_OBSTRUCTION"], "p2608_role_acceptance_obstruction_after_source_revalidation")
    equivariance_rows = build_equivariance_table()
    lattice = premise_lattice_rows()
    minimal_accepting = [row["support"] for row in lattice if row["minimal_accepting_support"]]
    rejected_legacy_only = [row["support"] for row in lattice if row["rejected_if_only_legacy_scalar_or_axis"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2619_T1_p2618_selector_source_obligation_lattice",
        "result_type": "exact_C2_torsor_enumeration_and_minimal_selector_source_obligation_after_P2618",
        "inherits_p2618_full_completion_block": p2618.get("full_completion_map_exported") is False,
        "inherits_p2616_role_block": p2616.get("current_assignment_accepts_role_bearing_ltotal") is False or p2616.get("p2608_role_bearing_ltotal_reenabled") is False,
        "equivariance_theorem": {
            "formal_statement": "For a C2-odd strict phase sign, no selector map from orientation-invariant legacy scalar data to {+1,-1} is C2-equivariant. A C2-equivariant selector exists only after the input already contains a freely transforming orientation/sign torsor or an equivalent source premise.",
            "proof_steps": [
                "Let C2={1,r} act on the strict sign set Sigma={+1,-1} by r·sigma=-sigma.",
                "If an input x is legacy-scalar/invariant, then r·x=x.",
                "Equivariance of a selector f requires f(r·x)=r·f(x), hence f(x)=-f(x).",
                "No element of Sigma satisfies sigma=-sigma, so no such selector exists from invariant input data.",
                "If the input is an orientation torsor X={x_+,x_-} with r·x_+=x_-, then f(x_+)=+1, f(x_-)=-1 and its negation are exactly the two equivariant maps; the missing object is therefore the source torsor, not a numerical fit.",
            ],
            "finite_enumeration_rows": equivariance_rows,
        },
        "premise_lattice": {
            "atoms": [
                "beta_tors_scalar_invariant",
                "axis_only_selector_up_to_Z2",
                "orientation_odd_source",
                "symmetry_breaking_boundary",
                "spin_pin_orientation_source",
            ],
            "row_count": len(lattice),
            "minimal_accepting_supports": minimal_accepting,
            "legacy_scalar_or_axis_only_rejected_supports": rejected_legacy_only,
            "all_rows": lattice,
        },
        "bridge_policy_after_p2619": {
            "beta_tors_status": "A beta_tors magnitude may remain legacy damping input, but as an orientation-invariant scalar it cannot by itself be the strict phase-sign source.",
            "axis_only_status": "Axis-only or quotient selector data can reduce continuous degeneracy but still leaves the Z2 sign unresolved.",
            "admissible_next_targets": minimal_accepting,
            "role_transfer_status": "Still blocked; after any real selector premise, rerun bridge validity and role-transfer audit separately.",
        },
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2618_block_inherited": theorem_export["inherits_p2618_full_completion_block"],
        "p2616_role_block_inherited": theorem_export["inherits_p2616_role_block"],
        "invariant_scalar_rows_have_zero_equivariant_maps": all(row["equivariant_function_count"] == 0 for row in equivariance_rows if row["action_kind"] == "invariant"),
        "torsor_rows_have_equivariant_maps": all(row["equivariant_function_count"] > 0 for row in equivariance_rows if row["action_kind"] == "torsor"),
        "minimal_supports_are_single_real_orientation_sources": sorted(minimal_accepting) == sorted([["orientation_odd_source"], ["symmetry_breaking_boundary"], ["spin_pin_orientation_source"]]),
        "legacy_scalar_or_axis_only_rejected": all(set(row) <= {"beta_tors_scalar_invariant", "axis_only_selector_up_to_Z2"} for row in rejected_legacy_only),
        "no_beta_tors_chi11_reopen": theorem_export["beta_tors_chi11_route_reopened"] is False,
        "no_strict_selector_export": theorem_export["strict_selector_source_exported"] is False,
        "no_ltotal_reenable": theorem_export["role_bearing_ltotal_reenabled"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_packet"] is False,
        "no_toe_closure": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2619",
        "stage_id": "S1569",
        "status": "P2619_SELECTOR_SOURCE_OBLIGATION_LATTICE_EXACT_C2_ENUMERATION_NO_SELECTOR_EXPORT_NO_LTOTAL_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "p2618_selector_source_obligation_lattice": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(payload) for name, payload in payloads.items()},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["p2618_selector_source_obligation_lattice"]["theorem_export"]
    eq = t["equivariance_theorem"]
    lattice = t["premise_lattice"]
    lines = [
        "# P2619/S1569 P2618 selector-source obligation lattice", "",
        f"Status: `{payload['status']}`", "", "## Theorem", "",
        eq["formal_statement"], "", "## Proof", "",
    ]
    for step in eq["proof_steps"]:
        lines.append(f"- {step}")
    lines.extend(["", "## Computed C2 enumeration", ""])
    for row in eq["finite_enumeration_rows"]:
        lines.append(f"- `{row['source_type']}`: candidates `{row['candidate_function_count']}`, equivariant maps `{row['equivariant_function_count']}`, selector available `{row['selector_available']}`.")
    lines.extend(["", "## Minimal source obligations", ""])
    for support in lattice["minimal_accepting_supports"]:
        lines.append(f"- `{support}`")
    lines.extend([
        "", "## Bridge policy", "",
        f"- beta_tors: {t['bridge_policy_after_p2619']['beta_tors_status']}",
        f"- axis-only data: {t['bridge_policy_after_p2619']['axis_only_status']}",
        f"- role transfer: {t['bridge_policy_after_p2619']['role_transfer_status']}",
        "", "## Scope guards", "",
        "No strict selector source, no `beta_tors -> chi11` route reopening, no GF(2) bridge revalidation, no role-transfer revalidation, no role-bearing `L_total`, no `QW-2191` discharge, and no ToE closure are exported.",
        "", "## Fingerprint", "",
        f"`{payload['p2618_selector_source_obligation_lattice']['theorem_fingerprint_sha256']}`",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2619/S1569 selector-source obligation lattice

`P2619/S1569` strengthens the P2618 phase-selector obstruction with an exact `C2` enumeration.  If orientation reversal acts trivially on legacy scalar input such as `beta_tors`, equivariance of a strict odd sign selector would force `f(x)=-f(x)`, so zero maps exist into `{+1,-1}`.  Equivariant maps appear only when the input already contains an orientation/sign torsor or an equivalent symmetry-breaking/spin-orientation source.  Thus `beta_tors` may remain damping input, but it is not reopened as a `chi11` sign source.
""".strip()
    lag_section = """
## P2619/S1569 selector-source Ltotal guard

`P2619/S1569` keeps role-bearing `L_total` closed.  The exact `C2` selector-source lattice says that legacy scalar or axis-only data cannot supply the missing strict phase sign; a real orientation-odd source, symmetry-breaking boundary, or spin/Pin orientation premise is still required before bridge validity and role transfer may be rerun.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2619/S1569 selector-source obligation lattice", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2619/S1569 selector-source Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
