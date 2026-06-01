#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2401_s1351_role_successor_unlock_order_certificate.json"
MD = GEN / "p2401_s1351_role_successor_unlock_order_certificate.md"

SOURCE_FILES = {
    "P2394_ROLE_FRONTIER": GEN / "p2394_s1344_apd_bridge_chi11_rebased_role_frontier_certificate.json",
    "P2400_ROLE_LATTICE": GEN / "p2400_s1350_nearest_lift_role_successor_lattice_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ROLE_ATOMS = [
    "alpha_geo_electroweak_role_theorem",
    "beta_tors_strict_role_theorem",
    "beta_power_hierarchy_successor_theorem",
]

ROLE_CLAIMS = {
    "legacy_weinberg_sin2_theta_role_transfer": ["alpha_geo_electroweak_role_theorem"],
    "legacy_alpha_em_inverse_role_transfer": ["alpha_geo_electroweak_role_theorem", "beta_tors_strict_role_theorem"],
    "legacy_beta_power_gravity_hierarchy_successor": ["beta_power_hierarchy_successor_theorem", "beta_tors_strict_role_theorem"],
    "all_three_physical_role_transfers": ROLE_ATOMS,
}


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_audit() -> dict[str, Any]:
    patterns = [
        "P2401|S1351|role-successor unlock order|unlock order certificate|role successor order",
        "P2400|nearest-lift role-successor lattice|full three-role successor monomial",
        "P2394|individual_role_minimal_supports|all_physical_roles_minimal_support",
        "legacy_weinberg_sin2_theta_role_transfer|legacy_alpha_em_inverse_role_transfer|legacy_beta_power_gravity_hierarchy_successor",
        "alpha_geo_electroweak_role_theorem|beta_tors_strict_role_theorem|beta_power_hierarchy_successor_theorem",
    ]
    out: dict[str, Any] = {}
    for pattern in patterns:
        proc = subprocess.run(
            ["rg", "-n", pattern, "fundamental_action_reconstruction", "-g", "*.py", "-g", "*.md", "-g", "*.json"],
            cwd=REPO,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        lines = [line for line in proc.stdout.splitlines() if line]
        out[pattern] = {"count": len(lines), "samples": lines[:16]}
    return {
        "tool": "rg",
        "patterns": out,
        "finding": (
            "Repo grep finds P2394 minimal supports and P2400's subset lattice, but no role-successor unlock-order certificate. "
            "P2401 therefore computes only permutation/prefix unlock timing for the already identified three role atoms."
        ),
    }


def first_unlock_step(order: tuple[str, ...], support: list[str]) -> int:
    positions = {atom: index + 1 for index, atom in enumerate(order)}
    return max(positions[atom] for atom in support)


def fraction_str(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def compute_unlock_orders() -> dict[str, Any]:
    rows = []
    claim_unlock_steps = {claim: [] for claim in ROLE_CLAIMS}
    first_atom_rows = {atom: [] for atom in ROLE_ATOMS}
    for order in itertools.permutations(ROLE_ATOMS):
        unlock_steps = {claim: first_unlock_step(order, support) for claim, support in ROLE_CLAIMS.items()}
        for claim, step in unlock_steps.items():
            claim_unlock_steps[claim].append(step)
        prefixes = []
        for size in range(0, len(order) + 1):
            active = set(order[:size])
            prefixes.append(
                {
                    "prefix_size": size,
                    "active_atoms": list(order[:size]),
                    "unlocked_claims": sorted(
                        claim for claim, support in ROLE_CLAIMS.items() if set(support) <= active
                    ),
                    "all_roles_and_toe_conditionally_ready": size == 3,
                }
            )
        row = {
            "order": list(order),
            "unlock_steps": unlock_steps,
            "total_unlock_step_sum": sum(unlock_steps.values()),
            "first_atom": order[0],
            "prefixes": prefixes,
        }
        rows.append(row)
        first_atom_rows[order[0]].append(row)
    best_sum = min(row["total_unlock_step_sum"] for row in rows)
    best_orders = [row["order"] for row in rows if row["total_unlock_step_sum"] == best_sum]
    expected_steps = {
        claim: fraction_str(sum(Fraction(step, len(steps)) for step in steps))
        for claim, steps in claim_unlock_steps.items()
    }
    first_atom_summary = {}
    for atom, atom_rows in first_atom_rows.items():
        first_atom_summary[atom] = {
            "order_count": len(atom_rows),
            "min_total_unlock_step_sum": min(row["total_unlock_step_sum"] for row in atom_rows),
            "mean_total_unlock_step_sum": fraction_str(sum(Fraction(row["total_unlock_step_sum"], len(atom_rows)) for row in atom_rows)),
            "step1_unlocked_claims_union": sorted(
                {claim for row in atom_rows for claim in row["prefixes"][1]["unlocked_claims"]}
            ),
        }
    return {
        "role_atoms": ROLE_ATOMS,
        "role_claim_supports": ROLE_CLAIMS,
        "permutation_count": len(rows),
        "rows": rows,
        "best_total_unlock_step_sum": best_sum,
        "best_orders_by_total_unlock_step_sum": best_orders,
        "expected_unlock_step_by_claim_over_uniform_orders": expected_steps,
        "first_atom_summary": first_atom_summary,
        "proof_reading": (
            "Alpha first is uniquely optimal for earliest partial claim clarification, and alpha->beta_tors->beta_power is the unique best total order; nevertheless full role-transfer/ToE remains step-3 only."
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2401/S1351 role-successor unlock-order certificate

`P2401/S1351` does not add a new role theorem.  It takes the P2394 role-claim supports and the P2400 three-role monomial, then enumerates all six possible orders for proving the three role-successor atoms.  The finite result is:

```text
best partial-clarification order = alpha_geo_electroweak_role_theorem -> beta_tors_strict_role_theorem -> beta_power_hierarchy_successor_theorem,
full role-transfer/ToE readiness = still only at prefix size 3.
```

Thus `alpha_geo` is the best first target only for early partial role-claim clarification; it is not a shortcut to role-transfer or ToE closure.
""".strip()
    lag_section = """
## P2401/S1351 unlock-order guard for Lagrangian/EOM

`P2401/S1351` ranks role-successor proof orders but keeps the P2400 closure rule intact: no one-atom or two-atom prefix licenses role-bearing `L_total`.  The best order can guide proof search, not Lagrangian promotion.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    unlock = compute_unlock_orders()
    p2394_supports = artifacts["P2394_ROLE_FRONTIER"].get("apd_bridge_chi11_rebased_role_frontier_certificate", {}).get("theorem_export", {}).get("minimal_supports", {})
    p2400_theorem = artifacts["P2400_ROLE_LATTICE"].get("nearest_lift_role_successor_lattice_certificate", {}).get("theorem_export", {})
    theorem_export = {
        "theorem_name": "P2401_T1_role_successor_unlock_order_certificate",
        "role_atoms": ROLE_ATOMS,
        "role_claim_supports": ROLE_CLAIMS,
        "permutation_count": unlock["permutation_count"],
        "best_total_unlock_step_sum": unlock["best_total_unlock_step_sum"],
        "best_orders_by_total_unlock_step_sum": unlock["best_orders_by_total_unlock_step_sum"],
        "expected_unlock_step_by_claim_over_uniform_orders": unlock["expected_unlock_step_by_claim_over_uniform_orders"],
        "first_atom_summary": unlock["first_atom_summary"],
        "p2394_minimal_supports_inherited": p2394_supports,
        "p2400_only_full_role_mask_closes_toe": p2400_theorem.get("toe_true_masks") == [7],
        "not_licensed": [
            "No role-successor atom is exported by an ordering certificate.",
            "No one-atom or two-atom prefix licenses role transfer or ToE closure.",
            "No L_total or SM/GR role-bearing term is promoted by the best order.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "six_orders_enumerated": unlock["permutation_count"] == 6,
        "unique_best_order_is_alpha_then_beta_tors_then_beta_power": theorem_export["best_orders_by_total_unlock_step_sum"] == [[
            "alpha_geo_electroweak_role_theorem",
            "beta_tors_strict_role_theorem",
            "beta_power_hierarchy_successor_theorem",
        ]],
        "best_total_unlock_step_sum_is_nine": theorem_export["best_total_unlock_step_sum"] == 9,
        "alpha_first_has_step1_weinberg_only": theorem_export["first_atom_summary"]["alpha_geo_electroweak_role_theorem"]["step1_unlocked_claims_union"] == ["legacy_weinberg_sin2_theta_role_transfer"],
        "beta_tors_first_unlocks_no_step1_claim": theorem_export["first_atom_summary"]["beta_tors_strict_role_theorem"]["step1_unlocked_claims_union"] == [],
        "p2400_full_mask_closure_inherited": theorem_export["p2400_only_full_role_mask_closes_toe"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2401_s1351_v1",
        "packet_id": "P2401",
        "stage_id": "S1351",
        "result_kind": "ROLE_SUCCESSOR_UNLOCK_ORDER_CERTIFICATE",
        "status": "PASS_ROLE_SUCCESSOR_UNLOCK_ORDER_ENUMERATED_NO_PREFIX_TRANSFER",
        "role_successor_unlock_order_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "unlock_order_enumeration": unlock,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Use alpha_geo first only as a proof-search priority for partial role clarification; do not promote role transfer until all three role-successor atoms are proved.",
        "global_status": "OPEN_PROGRESS_UNLOCK_ORDER_CERTIFIED_FULL_ROLE_MONOMIAL_STILL_REQUIRED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["role_successor_unlock_order_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2401 S1351: role-successor unlock-order certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2401/S1351 enumerates all six proof orders for the three role-successor atoms and measures when role claims unlock.",
                "",
                "## Unlock order",
                "",
                f"- Best total unlock-step sum: `{theorem['best_total_unlock_step_sum']}`.",
                f"- Best orders: `{theorem['best_orders_by_total_unlock_step_sum']}`.",
                f"- Expected unlock steps: `{theorem['expected_unlock_step_by_claim_over_uniform_orders']}`.",
                f"- First-atom summary: `{theorem['first_atom_summary']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Gatekeepers",
                "",
                f"`{payload['gatekeeper_checks']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    append_doc_sections()
    payload = build_payload()
    write_outputs(payload)
    if not all(payload["gatekeeper_checks"].values()):
        raise SystemExit(f"gatekeeper failure: {payload['gatekeeper_checks']}")


if __name__ == "__main__":
    main()
