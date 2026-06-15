#!/usr/bin/env python3
"""P2766/S1716: post-P2761-to-P2765 provenance-state reconciliation.

P2765 recommended a reconciliation rather than another hidden promotion.  This
script treats P2761-P2765 as a finite provenance ledger for the strict moment-map
coupling candidates lambda_sm_eff, kappa_gr_eff, and epsilon_mix_eff.  It checks
whether any named provenance atom has become theorem-level closed.  None has:
reference-cell/action-density, sign convention, field/curvature normalization,
and 4D nonproxy variational residual insertion all remain blocked on current
artifacts.  The result is a no-closure certificate, not an L_total promotion.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2761 = GEN / "p2761_s1711_kernel_moment_coupling_provenance_obstruction.json"
P2762 = GEN / "p2762_s1712_reference_cell_action_density_normalization_obstruction.json"
P2763 = GEN / "p2763_s1713_moment_coupling_sign_convention_conditional_obstruction.json"
P2764 = GEN / "p2764_s1714_field_curvature_normalization_compatibility_obstruction.json"
P2765 = GEN / "p2765_s1715_nonproxy_variational_insertion_residual_obstruction.json"
OUT = GEN / "p2766_s1716_post_moment_provenance_state_reconciliation.json"
MD = GEN / "p2766_s1716_post_moment_provenance_state_reconciliation.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "p2765_recommendation": r"post-P2761-to-P2765|P2697-P2765|genuinely new theorem fixing one named provenance atom",
    "reference_cell_atom": r"P2762|reference-cell action-density|canonical_reference_cell|action_density_normalization",
    "sign_atom": r"P2763|sign-convention|four-branch|unique_sign_branch",
    "field_curvature_atom": r"P2764|field-curvature|field/curvature normalization|positive-normalization",
    "nonproxy_atom": r"P2765|nonproxy variational|4D nonproxy|four_d_covariant",
    "closure_boundary": r"physical-coupling provenance|role-bearing L_total|ToE closure|bridge closure|role transfer|selector closure",
}

ATOM_SPECS = {
    "reference_cell_action_density": {
        "source": "P2762",
        "accepted_path": ["acceptance_matrix", "accepted_as_reference_cell_action_density_theorem"],
        "theorem_flag": "canonical_reference_cell/action_density_theorem_exported",
    },
    "sign_convention": {
        "source": "P2763",
        "accepted_path": ["acceptance_matrix", "accepted_as_sign_convention_theorem"],
        "theorem_flag": "sign_convention_theorem_exported",
    },
    "field_curvature_normalization": {
        "source": "P2764",
        "accepted_path": ["acceptance_matrix", "accepted_as_field_curvature_normalization_theorem"],
        "theorem_flag": "field_curvature_normalization_theorem_exported",
    },
    "nonproxy_variational_insertion": {
        "source": "P2765",
        "accepted_path": ["acceptance_matrix", "accepted_as_nonproxy_variational_insertion_theorem"],
        "theorem_flag": "nonproxy_variational_insertion_theorem_exported",
    },
}

NEGATIVE_EXPORT_FLAGS = [
    "new_provenance_atom_closed",
    "physical_coupling_provenance_theorem_exported",
    "role_bearing_ltotal_promoted",
    "selector_closure_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def get_path(obj: dict[str, Any], path: list[str]) -> Any:
    cur: Any = obj
    for key in path:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    return cur


def run_rg(pattern: str) -> list[str]:
    cmd = ["rg", "-n", "--glob", "!generated/*.json", pattern, "AGENTS.md", "fundamental_action_reconstruction"]
    proc = subprocess.run(cmd, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed for {pattern!r}: {proc.stderr}")
    return [line for line in proc.stdout.splitlines() if line.strip()]


def evidence_scan() -> dict[str, Any]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        hits = run_rg(pattern)
        rows.append({"lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return {"row_count": len(rows), "rows": rows, "hit_counts": {r["lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def provenance_atom_matrix(inputs: dict[str, dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for atom, spec in ATOM_SPECS.items():
        payload = inputs[spec["source"]]
        accepted_value = get_path(payload, spec["accepted_path"])
        rows.append({
            "atom": atom,
            "source": spec["source"],
            "input_status": payload.get("status"),
            "accepted_path": ".".join(spec["accepted_path"]),
            "accepted_value": accepted_value,
            "closed_by_current_artifacts": accepted_value is True,
            "blocker_summary": payload.get("decision", {}).get("reason") or payload.get("acceptance_matrix", {}).get("blocker"),
        })
    closed = [row["atom"] for row in rows if row["closed_by_current_artifacts"]]
    open_atoms = [row["atom"] for row in rows if not row["closed_by_current_artifacts"]]
    return {
        "rows": rows,
        "row_count": len(rows),
        "closed_atom_count": len(closed),
        "open_atom_count": len(open_atoms),
        "closed_atoms": closed,
        "open_atoms": open_atoms,
        "all_named_atoms_closed": len(open_atoms) == 0,
    }


def finite_closure_lattice(matrix: dict[str, Any]) -> dict[str, Any]:
    atoms = [row["atom"] for row in matrix["rows"]]
    current_vector = {row["atom"]: row["closed_by_current_artifacts"] for row in matrix["rows"]}
    profiles = []
    for values in itertools.product([False, True], repeat=len(atoms)):
        vector = dict(zip(atoms, values))
        profiles.append({
            "closed_count": sum(1 for v in values if v),
            "vector": vector,
            "would_license_physical_coupling_provenance": all(values),
        })
    current_profile = {
        "vector": current_vector,
        "closed_count": sum(1 for value in current_vector.values() if value),
        "would_license_physical_coupling_provenance": all(current_vector.values()),
    }
    return {
        "atom_order": atoms,
        "profile_count": len(profiles),
        "license_profile_count": sum(1 for p in profiles if p["would_license_physical_coupling_provenance"]),
        "current_profile": current_profile,
        "current_profile_is_license_profile": current_profile["would_license_physical_coupling_provenance"],
    }


def acceptance_matrix(scan: dict[str, Any], atom_matrix: dict[str, Any], lattice: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_evidence_present": scan["all_patterns_have_hits"],
        "all_four_named_atoms_audited": atom_matrix["row_count"] == 4,
        "no_named_atom_closed": atom_matrix["closed_atom_count"] == 0,
        "finite_lattice_exhausted": lattice["profile_count"] == 16,
        "only_full_four_atom_profile_licenses_physical_provenance": lattice["license_profile_count"] == 1,
        "current_profile_does_not_license_physical_provenance": not lattice["current_profile_is_license_profile"],
        "new_theorem_supplied_in_this_reconciliation": False,
        "physical_coupling_provenance_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_physical_coupling_provenance_theorem": False,
        "accepted_as_ltotal_promotion": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The reconciliation supplies no new theorem and all four named provenance atoms remain open; the finite closure lattice licenses physical-coupling provenance only in the all-atoms-closed profile, which is not the current profile.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    matrix = payload["provenance_atom_matrix"]
    lines = [
        "# P2766/S1716 post-P2761-to-P2765 provenance-state reconciliation",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Provenance atom ledger",
    ]
    for row in matrix["rows"]:
        lines.append(f"- {row['atom']} ({row['source']}): closed={row['closed_by_current_artifacts']}; accepted_value={row['accepted_value']}")
    lines.extend([
        "",
        "## Finite closure lattice",
        f"- profile_count={payload['finite_closure_lattice']['profile_count']}",
        f"- license_profile_count={payload['finite_closure_lattice']['license_profile_count']}",
        f"- current_profile_is_license_profile={payload['finite_closure_lattice']['current_profile_is_license_profile']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {"P2761": read_json(P2761), "P2762": read_json(P2762), "P2763": read_json(P2763), "P2764": read_json(P2764), "P2765": read_json(P2765)}
    scan = evidence_scan()
    atoms = provenance_atom_matrix(inputs)
    lattice = finite_closure_lattice(atoms)
    acceptance = acceptance_matrix(scan, atoms, lattice)
    payload = {
        "status": "P2766_POST_MOMENT_PROVENANCE_STATE_RECONCILIATION_NO_CLOSURE",
        "input_hashes": {key: sha(path) for key, path in {"P2761": P2761, "P2762": P2762, "P2763": P2763, "P2764": P2764, "P2765": P2765}.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_scope": "post-P2761-to-P2765 moment-map physical-coupling provenance state",
        "content_evidence_scan": scan,
        "provenance_atom_matrix": atoms,
        "finite_closure_lattice": lattice,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote lambda_sm_eff, kappa_gr_eff, or epsilon_mix_eff to physical Lagrangian couplings.  The next honest move must introduce one genuinely new theorem or artifact closing exactly one named atom (reference-cell/action-density, sign convention, field/curvature normalization, or 4D nonproxy residual rows).  If no such new object is supplied, preserve the P2697-P2766 no-closure certificate and pivot to a fresh broad state-map intake rather than replaying L_total promotion.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2766/S1716 post-moment provenance-state reconciliation", "## P2766/S1716 post-moment provenance-state reconciliation\n\n`P2766/S1716` reconciles the P2761-P2765 moment-map provenance sequence instead of adding an implicit `L_total` promotion.  The finite atom ledger has four named open atoms: reference-cell/action-density normalization, sign convention, field/curvature normalization, and 4D nonproxy variational residual insertion.  A 16-profile Boolean closure lattice confirms that physical-coupling provenance would be licensed only by the all-atoms-closed profile; the current profile closes zero atoms.  No physical-coupling provenance theorem, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2766/S1716 post-moment provenance-state Ltotal guard", "## P2766/S1716 post-moment provenance-state Ltotal guard\n\n`P2766/S1716` adds no variational source term.  It records that P2761-P2765 collectively leave all named physical-coupling provenance atoms open, and a finite closure lattice does not license `lambda_sm_eff`, `kappa_gr_eff`, or `epsilon_mix_eff` as physical Lagrangian couplings.  Therefore it cannot promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current post-moment provenance-state reconciliation guardrail (P2766/S1716, 2026-06-15)", "## Current post-moment provenance-state reconciliation guardrail (P2766/S1716, 2026-06-15)\n\n- P2766 reconciles the P2761-P2765 physical-coupling provenance sequence for `lambda_sm_eff`, `kappa_gr_eff`, and `epsilon_mix_eff`.\n- The finite four-atom ledger keeps reference-cell/action-density, sign convention, field/curvature normalization, and 4D nonproxy variational residual insertion all open; the 16-profile closure lattice licenses physical-coupling provenance only in the all-atoms-closed profile, not in the current zero-closed-atom profile.\n- Do not promote moment-map coefficients to physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.  A next admissible move must introduce one genuinely new theorem/artifact closing exactly one named provenance atom, or pivot to a fresh broad state-map intake while preserving the P2697-P2766 no-closure certificate.\n")
    return payload


if __name__ == "__main__":
    main()
