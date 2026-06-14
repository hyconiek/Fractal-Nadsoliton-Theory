#!/usr/bin/env python3
"""P2715/S1665: Aut-equivariant scalar-source to orientation-torsor no-go.

P2714 showed that the P2708 orientation torsor has no Aut(Z12)-compatible
section.  The next professorial/theoretical-physics question is whether an
ordinary strict scalar source (entropy scale, alpha normalization, beta damping,
or scalar Lagrangian datum) can induce such a section.  P2715 checks the finite
equivariance condition: a source with trivial Aut action can map to the
orientation torsor only if the torsor has a fixed point.  It does not.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2715_s1665_aut_equivariant_scalar_source_to_orientation_torsor_no_go.json"
MD = GEN / "p2715_s1665_aut_equivariant_scalar_source_to_orientation_torsor_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

UNITS = [1, 5, 7, 11]
TORSOR = [-1, 1]

INPUTS = {
    "P2713_INTAKE_GATE": GEN / "p2713_s1663_post_p2712_new_typed_object_intake_gate_certificate.json",
    "P2714_ORIENTATION_TORSOR": GEN / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.json",
    "P2689_UV_UNIT": GEN / "p2689_s1639_entropy_uv_unit_obligation_audit.json",
    "P2691_ALPHA_GEO": GEN / "p2691_s1641_alpha_geo_amplitude_role_safe_source_audit.json",
    "P2692_BETA_SOURCE": GEN / "p2692_s1642_target_independent_positive_beta_z_beta_source_audit.json",
    "P2687_ANISOTROPIC_SOURCE": GEN / "p2687_s1637_anisotropic_source_class_audit.json",
}

SCALAR_SOURCE_CLASSES = [
    {
        "source_class": "entropy_or_uv_scale_scalar",
        "representative_artifact": "P2689_UV_UNIT",
        "aut_action": "trivial_scalar",
        "physical_reason": "A positive scale/reference-cell datum has no orientation parity by itself.",
    },
    {
        "source_class": "alpha_geo_amplitude_normalization_scalar",
        "representative_artifact": "P2691_ALPHA_GEO",
        "aut_action": "trivial_scalar",
        "physical_reason": "A scalar amplitude normalization is invariant under inversion unless a separate pseudoscalar law is supplied.",
    },
    {
        "source_class": "positive_beta_or_z_beta_damping_scalar",
        "representative_artifact": "P2692_BETA_SOURCE",
        "aut_action": "trivial_scalar",
        "physical_reason": "Positive damping/compression parameters are orientation-blind scalars.",
    },
    {
        "source_class": "scalar_lagrangian_density_or_metric_residual",
        "representative_artifact": "P2687_ANISOTROPIC_SOURCE",
        "aut_action": "trivial_scalar",
        "physical_reason": "A scalar Lagrangian density does not pick a Z12 orientation unless augmented by a strict pseudoscalar/chiral source.",
    },
]

NEGATIVE_EXPORT_FLAGS = [
    "scalar_source_breaks_orientation_torsor",
    "aut_equivariant_scalar_to_torsor_map_exported",
    "strict_pseudoscalar_or_chiral_source_exported",
    "strict_mechanism_fixing_lambda_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def torsor_action(unit: int, orientation: int) -> int:
    if unit in (1, 5):
        return orientation
    if unit in (7, 11):
        return -orientation
    raise ValueError(f"not a Z12 unit: {unit}")


def trivial_domain_action(unit: int, point: str) -> str:
    if unit not in UNITS:
        raise ValueError(f"not a Z12 unit: {unit}")
    return point


def equivariant_maps_from_domain(domain_points: list[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for images in itertools.product(TORSOR, repeat=len(domain_points)):
        mapping = dict(zip(domain_points, images))
        failures = []
        for unit in UNITS:
            for point in domain_points:
                lhs = mapping[trivial_domain_action(unit, point)]
                rhs = torsor_action(unit, mapping[point])
                if lhs != rhs:
                    failures.append({
                        "unit": unit,
                        "domain_point": point,
                        "map_value": "+omega" if mapping[point] == 1 else "-omega",
                        "torsor_action_value": "+omega" if rhs == 1 else "-omega",
                    })
        rows.append({
            "mapping": {point: "+omega" if image == 1 else "-omega" for point, image in mapping.items()},
            "aut_equivariant": not failures,
            "failure_count": len(failures),
            "failure_sample": failures[:4],
        })
    return rows


def scalar_source_matrix() -> list[dict[str, Any]]:
    rows = []
    one_point_maps = equivariant_maps_from_domain(["scalar_value"])
    two_point_maps = equivariant_maps_from_domain(["scalar_branch_a", "scalar_branch_b"])
    one_point_admitted = any(row["aut_equivariant"] for row in one_point_maps)
    two_point_admitted = any(row["aut_equivariant"] for row in two_point_maps)
    for source in SCALAR_SOURCE_CLASSES:
        rows.append({
            **source,
            "domain_model": "Aut-trivial scalar domain",
            "one_point_equivariant_map_exists": one_point_admitted,
            "two_point_trivial_branch_equivariant_map_exists": two_point_admitted,
            "breaks_orientation_torsor": False,
            "blocker": "Aut-trivial scalar data can map equivariantly to the orientation torsor only through an Aut-fixed torsor point; P2714/P2715 find none.",
        })
    return rows


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2715/S1665 Aut-equivariant scalar-source to orientation-torsor no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite equivariance result",
        f"- one_point_equivariant_maps={payload['finite_equivariance_summary']['one_point_equivariant_maps']}",
        f"- two_point_trivial_branch_equivariant_maps={payload['finite_equivariance_summary']['two_point_trivial_branch_equivariant_maps']}",
        "",
        "## Scalar source matrix",
    ]
    for row in payload["scalar_source_matrix"]:
        lines.append(f"- `{row['source_class']}`: breaks_orientation_torsor={row['breaks_orientation_torsor']}. {row['blocker']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    one_point_maps = equivariant_maps_from_domain(["scalar_value"])
    two_point_maps = equivariant_maps_from_domain(["scalar_branch_a", "scalar_branch_b"])
    matrix = scalar_source_matrix()
    one_count = sum(1 for row in one_point_maps if row["aut_equivariant"])
    two_count = sum(1 for row in two_point_maps if row["aut_equivariant"])
    no_unlock = one_count == 0 and two_count == 0 and all(not row["breaks_orientation_torsor"] for row in matrix)
    payload = {
        "status": "P2715_AUT_TRIVIAL_SCALAR_SOURCE_TO_ORIENTATION_TORSOR_NO_GO" if no_unlock else "P2715_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: read_json(path).get("status") for key, path in INPUTS.items()},
        "finite_equivariance_summary": {
            "aut_group": "Aut(Z12)=U(12)={1,5,7,11}",
            "torsor": ["+omega", "-omega"],
            "orientation_reversing_units": [7, 11],
            "one_point_candidate_maps": len(one_point_maps),
            "one_point_equivariant_maps": one_count,
            "two_point_trivial_branch_candidate_maps": len(two_point_maps),
            "two_point_trivial_branch_equivariant_maps": two_count,
        },
        "one_point_map_rows": one_point_maps,
        "two_point_trivial_branch_map_rows": two_point_maps,
        "scalar_source_matrix": matrix,
        "decision": {
            "new_source_class_tested": "Aut-trivial strict scalar sources as torsor breakers",
            "scalar_source_breaks_orientation_torsor": False,
            "aut_equivariant_scalar_to_torsor_map_exported": False,
            "strict_pseudoscalar_or_chiral_source_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2715 tests whether ordinary strict scalar data can break the P2714 orientation torsor.  The finite equivariance calculation finds zero Aut-equivariant maps from a trivial scalar source domain to the orientation torsor, because orientation-reversing units 7 and 11 require an Aut-fixed torsor point and none exists.  Entropy/UV scale, alpha_geo amplitude, positive beta damping, and scalar Lagrangian data therefore remain orientation-blind unless a new strict pseudoscalar/chiral source is supplied.",
            "next_honest_step": "Do not seek the missing sign in Aut-trivial scalar quantities.  The next admissible move must supply a genuinely strict inversion-odd pseudoscalar/chiral source with a nonzero signed value coupled to the orientation torsor, or pivot to a different typed object outside the closed lanes.  Without that, preserve the P2697-P2715 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2715/S1665 Aut-equivariant scalar-source to orientation-torsor no-go", "## P2715/S1665 Aut-equivariant scalar-source to orientation-torsor no-go\n\n`P2715/S1665` tests whether ordinary strict scalar data can break the P2714 orientation torsor.  The finite equivariance calculation finds zero `Aut(Z12)`-equivariant maps from Aut-trivial scalar domains to the `+omega/-omega` torsor because orientation-reversing units `7` and `11` require an Aut-fixed torsor point.  Entropy/UV scale, `alpha_geo` amplitude, positive `beta` damping, and scalar Lagrangian data therefore remain orientation-blind unless a new strict inversion-odd pseudoscalar/chiral source is supplied; no `QW-2191`, pair12 strict-core, bridge, role-transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2715/S1665 scalar-source torsor no-go Ltotal guard", "## P2715/S1665 scalar-source torsor no-go Ltotal guard\n\n`P2715/S1665` is an equivariance obstruction for scalar sources, not a variational source construction.  Scalar Lagrangian data are Aut-trivial for the orientation torsor unless a strict pseudoscalar/chiral term is newly exported, so this does not promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current Aut-equivariant scalar-source torsor no-go guardrail (P2715/S1665, 2026-06-14)", "## Current Aut-equivariant scalar-source torsor no-go guardrail (P2715/S1665, 2026-06-14)\n\n- P2715 tests whether Aut-trivial strict scalar data can break the P2714 orientation torsor and finds zero Aut-equivariant maps to the `+omega/-omega` torsor.\n- Entropy/UV scale, `alpha_geo` amplitude, positive `beta` damping, and scalar Lagrangian data remain orientation-blind unless a genuinely strict inversion-odd pseudoscalar/chiral source with a nonzero signed value is exported.\n- Do not replay scalar-source attempts as `QW-2191` discharge, selector closure, pair12 strict-core, role transfer, bridge closure, `L_total`, or ToE; a next admissible move requires such a pseudoscalar/chiral source or a different genuinely new typed object outside closed lanes.\n")
    return payload


if __name__ == "__main__":
    main()
