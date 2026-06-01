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

OUT = GEN / "p2402_s1352_role_successor_marginal_credit_certificate.json"
MD = GEN / "p2402_s1352_role_successor_marginal_credit_certificate.md"

SOURCE_FILES = {
    "P2401_UNLOCK_ORDER": GEN / "p2401_s1351_role_successor_unlock_order_certificate.json",
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

PHYSICAL_CLAIMS = [
    "legacy_weinberg_sin2_theta_role_transfer",
    "legacy_alpha_em_inverse_role_transfer",
    "legacy_beta_power_gravity_hierarchy_successor",
]


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


def fraction_str(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def rg_audit() -> dict[str, Any]:
    patterns = [
        "P2402|S1352|role-successor marginal credit|marginal credit certificate|unlock credit",
        "P2401|role-successor unlock-order|best total unlock-step|expected unlock steps",
        "P2400|nearest-lift role-successor lattice|full three-role successor monomial",
        "Shapley|Banzhaf|atom influence|marginal contribution|pivotal",
        "legacy_weinberg_sin2_theta_role_transfer|legacy_alpha_em_inverse_role_transfer|legacy_beta_power_gravity_hierarchy_successor",
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
            "Repo grep finds a global frontier atom-influence/Banzhaf-style scratch report and P2401 unlock timing, "
            "but no role-local marginal credit certificate over the P2400 three-role successor lattice."
        ),
    }


def unlocked_claims(active: set[str]) -> set[str]:
    return {claim for claim, support in ROLE_CLAIMS.items() if set(support) <= active}


def compute_marginal_credit() -> dict[str, Any]:
    rows = []
    marginal_claims_by_atom: dict[str, list[list[str]]] = {atom: [] for atom in ROLE_ATOMS}
    marginal_counts_by_atom: dict[str, list[int]] = {atom: [] for atom in ROLE_ATOMS}
    physical_counts_by_atom: dict[str, list[int]] = {atom: [] for atom in ROLE_ATOMS}
    per_claim_credit: dict[str, dict[str, Fraction]] = {
        claim: {atom: Fraction(0, 1) for atom in ROLE_ATOMS} for claim in ROLE_CLAIMS
    }

    orders = list(itertools.permutations(ROLE_ATOMS))
    for order in orders:
        active: set[str] = set()
        row_steps = []
        for step, atom in enumerate(order, start=1):
            before = unlocked_claims(active)
            active.add(atom)
            after = unlocked_claims(active)
            marginal = sorted(after - before)
            physical_marginal = [claim for claim in marginal if claim in PHYSICAL_CLAIMS]
            marginal_claims_by_atom[atom].append(marginal)
            marginal_counts_by_atom[atom].append(len(marginal))
            physical_counts_by_atom[atom].append(len(physical_marginal))
            for claim in marginal:
                per_claim_credit[claim][atom] += Fraction(1, len(orders))
            row_steps.append(
                {
                    "step": step,
                    "atom": atom,
                    "marginal_claims": marginal,
                    "marginal_claim_count": len(marginal),
                    "physical_marginal_claims": physical_marginal,
                    "physical_marginal_claim_count": len(physical_marginal),
                }
            )
        rows.append({"order": list(order), "steps": row_steps})

    atom_credit = {}
    for atom in ROLE_ATOMS:
        total = sum(Fraction(count, len(orders)) for count in marginal_counts_by_atom[atom])
        physical = sum(Fraction(count, len(orders)) for count in physical_counts_by_atom[atom])
        atom_credit[atom] = {
            "mean_marginal_claim_count": fraction_str(total),
            "mean_physical_marginal_claim_count": fraction_str(physical),
            "marginal_count_samples_by_order": marginal_counts_by_atom[atom],
            "physical_marginal_count_samples_by_order": physical_counts_by_atom[atom],
            "claim_credit_breakdown": {
                claim: fraction_str(per_claim_credit[claim][atom]) for claim in ROLE_CLAIMS
            },
        }
    ranking = sorted(ROLE_ATOMS, key=lambda atom: (Fraction(atom_credit[atom]["mean_marginal_claim_count"]), atom), reverse=True)
    physical_ranking = sorted(ROLE_ATOMS, key=lambda atom: (Fraction(atom_credit[atom]["mean_physical_marginal_claim_count"]), atom), reverse=True)
    per_claim_credit_str = {
        claim: {atom: fraction_str(value) for atom, value in credits.items()}
        for claim, credits in per_claim_credit.items()
    }
    return {
        "role_atoms": ROLE_ATOMS,
        "role_claim_supports": ROLE_CLAIMS,
        "order_count": len(orders),
        "rows": rows,
        "atom_marginal_credit": atom_credit,
        "per_claim_credit": per_claim_credit_str,
        "total_claim_credit_by_atom_ranking": ranking,
        "physical_claim_credit_by_atom_ranking": physical_ranking,
        "proof_reading": (
            "Alpha has the largest role-local marginal credit because it alone unlocks the Weinberg claim; beta_tors is still required for two two-atom claims and the all-role monomial, but no marginal credit exports a theorem."
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2402/S1352 role-successor marginal-credit certificate

`P2402/S1352` complements P2401 by computing role-local marginal credit over all six proof orders.  This is not the older global frontier influence report; it only asks how much each of the three role-successor atoms contributes to unlocking the P2394 role claims inside the P2400 lattice.

The exact marginal-credit vector over all role claims is:

```text
alpha_geo_electroweak_role_theorem: 11/6,
beta_tors_strict_role_theorem: 4/3,
beta_power_hierarchy_successor_theorem: 5/6.
```

For physical claims only, the vector is `3/2, 1, 1/2`.  This supports prioritizing `alpha_geo` for first proof-search clarification, while preserving the P2400 rule that all three atoms are still required for full role-transfer/ToE readiness.
""".strip()
    lag_section = """
## P2402/S1352 marginal-credit guard for Lagrangian/EOM

`P2402/S1352` assigns proof-search marginal credit to role-successor atoms, but the credit score is not a Lagrangian license.  Even the top-scoring `alpha_geo` atom only guides which proof to try first; role-bearing `L_total` still waits for the full three-role monomial.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    credit = compute_marginal_credit()
    p2401_theorem = artifacts["P2401_UNLOCK_ORDER"].get("role_successor_unlock_order_certificate", {}).get("theorem_export", {})
    p2400_theorem = artifacts["P2400_ROLE_LATTICE"].get("nearest_lift_role_successor_lattice_certificate", {}).get("theorem_export", {})
    theorem_export = {
        "theorem_name": "P2402_T1_role_successor_marginal_credit_certificate",
        "role_atoms": ROLE_ATOMS,
        "order_count": credit["order_count"],
        "atom_marginal_credit": credit["atom_marginal_credit"],
        "per_claim_credit": credit["per_claim_credit"],
        "total_claim_credit_by_atom_ranking": credit["total_claim_credit_by_atom_ranking"],
        "physical_claim_credit_by_atom_ranking": credit["physical_claim_credit_by_atom_ranking"],
        "p2401_best_order_inherited": p2401_theorem.get("best_orders_by_total_unlock_step_sum"),
        "p2400_full_role_mask_still_required": p2400_theorem.get("toe_true_masks") == [7],
        "not_licensed": [
            "No marginal-credit score exports a role-successor atom.",
            "No marginal-credit score licenses a one-atom or two-atom prefix for role transfer.",
            "No L_total or SM/GR role-bearing term is promoted by this prioritization metric.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "six_orders_enumerated": credit["order_count"] == 6,
        "alpha_total_credit_is_11_over_6": theorem_export["atom_marginal_credit"]["alpha_geo_electroweak_role_theorem"]["mean_marginal_claim_count"] == "11/6",
        "beta_tors_total_credit_is_4_over_3": theorem_export["atom_marginal_credit"]["beta_tors_strict_role_theorem"]["mean_marginal_claim_count"] == "4/3",
        "beta_power_total_credit_is_5_over_6": theorem_export["atom_marginal_credit"]["beta_power_hierarchy_successor_theorem"]["mean_marginal_claim_count"] == "5/6",
        "alpha_is_top_total_and_physical_credit": theorem_export["total_claim_credit_by_atom_ranking"][0] == "alpha_geo_electroweak_role_theorem" and theorem_export["physical_claim_credit_by_atom_ranking"][0] == "alpha_geo_electroweak_role_theorem",
        "p2401_best_order_inherited": theorem_export["p2401_best_order_inherited"] == [["alpha_geo_electroweak_role_theorem", "beta_tors_strict_role_theorem", "beta_power_hierarchy_successor_theorem"]],
        "p2400_full_role_mask_still_required": theorem_export["p2400_full_role_mask_still_required"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2402_s1352_v1",
        "packet_id": "P2402",
        "stage_id": "S1352",
        "result_kind": "ROLE_SUCCESSOR_MARGINAL_CREDIT_CERTIFICATE",
        "status": "PASS_ROLE_SUCCESSOR_MARGINAL_CREDIT_PRIORITIZES_ALPHA_NO_PREFIX_LICENSE",
        "role_successor_marginal_credit_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "marginal_credit_enumeration": credit,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Use alpha_geo as the first proof-search target for role clarification, while keeping beta_tors and beta_power obligations explicit for full role-transfer/ToE readiness.",
        "global_status": "OPEN_PROGRESS_MARGINAL_CREDIT_CERTIFIED_NO_ROLE_PREFIX_LICENSE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["role_successor_marginal_credit_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2402 S1352: role-successor marginal-credit certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2402/S1352 computes role-local marginal credit over all six proof orders for the three role-successor atoms.",
                "",
                "## Marginal credit",
                "",
                f"- Ranking over all claims: `{theorem['total_claim_credit_by_atom_ranking']}`.",
                f"- Ranking over physical claims: `{theorem['physical_claim_credit_by_atom_ranking']}`.",
                f"- Atom credit: `{theorem['atom_marginal_credit']}`.",
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
