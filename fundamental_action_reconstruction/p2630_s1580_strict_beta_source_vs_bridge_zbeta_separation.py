#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import product
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.json"
MD = GEN / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.md"

REPO_ROOT = REPO
P2603_JSON = GEN / "p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem.json"
P2608_JSON = GEN / "p2608_s1558_strict_damping_role_transfer_to_ltotal_theorem.json"
P2610_JSON = GEN / "p2610_s1560_p2601_p2608_critical_revalidation_audit.json"
P2625_JSON = GEN / "p2625_s1575_nonlinear_damping_completion_source_classification.json"
P2629_JSON = GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.json"

NEGATIVE_EXPORT_FLAGS = [
    "legacy_uv_normalization_source_exported",
    "bridge_positive_zbeta_source_exported",
    "nonlinear_damping_completion_source_exported",
    "p2608_role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "qw2191_discharged",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.lean", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO_ROOT,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: research content and source-type terms, not ticket IDs.
    patterns = {
        "strict_internal_beta_source_content": (
            "strict damping beta|beta/eta source|strict beta|fractal codimension slope|"
            "prime-log|slope source|D_f-1|beta=1"
        ),
        "legacy_uv_bridge_coefficient_content": (
            "legacy beta_tors|beta_tors|UV normalization|beta_uv|Z_beta|beta/beta_tors|"
            "positive coefficient source|coefficient bridge"
        ),
        "source_type_separation_content": (
            "source-type separation|internal source|bridge source|normalization-gauge|"
            "target-independent|not a bridge theorem|not.*source theorem"
        ),
        "role_transfer_revalidation_content": (
            "role-transfer theorem|role-bearing L_total|P2608|revalidation|quarantine|"
            "legacy-to-strict bridge|completion map"
        ),
        "closure_guard_content": (
            "QW-2191|ToE closure|full bridge|nonlinear_damping_completion_source|"
            "positive_beta_renormalization_source"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for strict beta source vs bridge Z_beta separation", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def upstream_summary(p2603: dict[str, Any], p2608: dict[str, Any], p2610: dict[str, Any], p2625: dict[str, Any], p2629: dict[str, Any]) -> dict[str, Any]:
    p2603_text = json.dumps(p2603, sort_keys=True)
    p2608_text = json.dumps(p2608, sort_keys=True)
    p2610_text = json.dumps(p2610, sort_keys=True)
    p2625_text = json.dumps(p2625, sort_keys=True)
    return {
        "p2603_claims_strict_internal_beta_eta_source": "strict_damping_beta_eta_source_exported" in p2603_text and "true" in p2603_text.lower(),
        "p2603_scope_excludes_bridge": "legacy-to-strict bridge" in p2603_text or "bridge" in p2603_text,
        "p2608_role_transfer_was_later_quarantined_or_revalidated_needed": "quarantine" in p2610_text.lower() or "revalidation" in p2610_text.lower(),
        "p2625_requires_positive_zbeta_source": "positive_beta_renormalization_source" in p2625_text,
        "p2629_accepts_positive_zbeta_source": bool(p2629.get("exact_source_gate", {}).get("accepts_positive_beta_renormalization_source", False)),
        "p2629_invariant_ratio": p2629.get("normalization_orbit_certificate", {}).get("invariant_ratio"),
        "p2629_failed_gates": p2629.get("exact_source_gate", {}).get("failed_gates", []),
        "p2608_text_mentions_role_bearing_ltotal": "role-bearing" in p2608_text and "L_total" in p2608_text,
    }


def truth_table() -> dict[str, Any]:
    rows = []
    for strict_internal_beta_source, legacy_uv_normalization_source, invariant_match in product([False, True], repeat=3):
        accepts = strict_internal_beta_source and legacy_uv_normalization_source and invariant_match
        rows.append({
            "strict_internal_beta_source": strict_internal_beta_source,
            "legacy_uv_normalization_source": legacy_uv_normalization_source,
            "normalization_invariant_match_beta_micro_over_beta_strict_equals_1": invariant_match,
            "accepts_bridge_positive_zbeta_source": accepts,
            "rejection_reason": [] if accepts else [name for name, value in {
                "missing_strict_internal_beta_source": strict_internal_beta_source,
                "missing_legacy_uv_normalization_source": legacy_uv_normalization_source,
                "normalization_invariant_mismatch": invariant_match,
            }.items() if not value],
        })
    accepting = [row for row in rows if row["accepts_bridge_positive_zbeta_source"]]
    return {
        "row_count": len(rows),
        "accepting_row_count": len(accepting),
        "unique_accepting_row": accepting[0],
        "current_assignment_after_p2629_even_if_p2603_granted": {
            "strict_internal_beta_source": True,
            "legacy_uv_normalization_source": False,
            "normalization_invariant_match_beta_micro_over_beta_strict_equals_1": False,
        },
        "current_accepts_bridge_positive_zbeta_source": False,
        "rows": rows,
        "theorem": (
            "A strict-internal beta source and a legacy-to-strict bridge coefficient source are different typed obligations.  "
            "Even granting an internal strict beta normalization, a bridge Z_beta source additionally requires an independent legacy/UV normalization "
            "source and a normalization-invariant match beta_micro/beta_strict=1.  Therefore P2603-style strict beta sourcehood cannot be silently reused as the P2625 positive_beta_renormalization_source."
        ),
    }


def source_type_matrix(summary: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "claim": "strict_internal_beta_eta_source",
            "best_prior_art": "P2603 strict damping beta/eta source normal form, later scoped by P2610 and bridge guardrails",
            "typed_layer": "strict-internal damping law",
            "can_supply_p2625_positive_zbeta_source": False,
            "reason": "It may identify beta/eta inside the strict kernel normalization, but it does not tie beta_tors=0.01 or beta_uv to beta=1 by an independent bridge normalization theorem.",
        },
        {
            "claim": "legacy_uv_normalization_source",
            "best_prior_art": "P2629 normalization-gauge audit",
            "typed_layer": "legacy/UV-to-strict coefficient bridge",
            "can_supply_p2625_positive_zbeta_source": False,
            "reason": "P2629 shows the absolute number Z_beta=100 is convention-sensitive unless beta_uv=0.01 is independently fixed and the invariant beta_micro/beta_strict mismatch is removed.",
        },
        {
            "claim": "role_transfer_or_ltotal_reuse",
            "best_prior_art": "P2608 role-transfer path, later guarded by P2610 and P2625-P2629",
            "typed_layer": "role-bearing L_total admission",
            "can_supply_p2625_positive_zbeta_source": False,
            "reason": "Role-transfer cannot precede the repaired bridge coefficient source; reusing old P2608 would reverse the dependency order.",
        },
    ]


def recommendation() -> str:
    return (
        "Next honest step: stop trying to recycle P2603/P2608 as the missing bridge coefficient source.  Either prove a new typed legacy/UV normalization theorem "
        "fixing beta_uv=beta_tors=0.01 and beta_micro/beta_strict=1 from nadsoliton dynamics, or formally downgrade the legacy-to-strict damping bridge to an approximate "
        "finite-domain statement before any role-transfer audit is reopened."
    )


def write_markdown(payload: dict[str, Any]) -> None:
    table = payload["bridge_source_truth_table"]
    lines = [
        "# P2630/S1580 strict beta source vs bridge Z_beta separation",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps source-type and bridge-coefficient research content rather than only names or ticket numbers.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Separation theorem",
        "",
        table["theorem"],
        "",
        f"Truth-table rows: `{table['row_count']}`; accepting rows: `{table['accepting_row_count']}`.",
        "",
        "Current assignment even if P2603 is granted:",
        "",
    ])
    for key, value in table["current_assignment_after_p2629_even_if_p2603_granted"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend([
        "",
        "## Verdict",
        "",
        "P2630 does not export `positive_beta_renormalization_source`.  It blocks silent reuse of strict-internal beta sourcehood or old role-transfer text as the missing legacy-to-strict bridge coefficient theorem.",
        "",
        "## Recommended next honest step",
        "",
        payload["recommended_next_honest_step"],
        "",
        f"Fingerprint: `{payload['fingerprint_sha256']}`",
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(parents=True, exist_ok=True)
    p2603 = load_json(P2603_JSON)
    p2608 = load_json(P2608_JSON)
    p2610 = load_json(P2610_JSON)
    p2625 = load_json(P2625_JSON)
    p2629 = load_json(P2629_JSON)
    summary = upstream_summary(p2603, p2608, p2610, p2625, p2629)
    table = truth_table()
    payload: dict[str, Any] = {
        "status": "P2630_STRICT_BETA_SOURCE_VS_BRIDGE_ZBETA_SEPARATION_NO_SOURCE_REUSE_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "sources": {
            "P2603_JSON": rel(P2603_JSON),
            "P2608_JSON": rel(P2608_JSON),
            "P2610_JSON": rel(P2610_JSON),
            "P2625_JSON": rel(P2625_JSON),
            "P2629_JSON": rel(P2629_JSON),
        },
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "upstream_summary": summary,
        "source_type_matrix": source_type_matrix(summary),
        "bridge_source_truth_table": table,
        "recommended_next_honest_step": recommendation(),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2630/S1580 strict beta source vs bridge Z_beta separation guard",
        "\n## P2630/S1580 strict beta source vs bridge Z_beta separation guard\n\n"
        "`P2630/S1580` separates strict-internal beta/eta sourcehood from the legacy-to-strict bridge coefficient obligation.  Even if a P2603-style internal strict beta source is granted, the P2625 `positive_beta_renormalization_source` additionally requires an independent legacy/UV normalization source fixing `beta_uv=beta_tors=0.01` and a normalization-invariant match `beta_micro/beta_strict=1`.  The finite truth table has one accepting row and the current assignment remains rejecting, so P2603/P2608 cannot be silently recycled to reopen bridge completion, role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2630/S1580 strict beta source vs bridge Z_beta Ltotal guard",
        "\n## P2630/S1580 strict beta source vs bridge Z_beta Ltotal guard\n\n"
        "`P2630/S1580` keeps role-bearing `L_total` closed.  A strict-internal beta normalization, an independent legacy/UV normalization, and the normalization-invariant micro/strict match are typed separately; old strict-source or role-transfer packets cannot stand in for the missing bridge coefficient theorem.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "accepting_rows": result["bridge_source_truth_table"]["accepting_row_count"],
        "current_accepts": result["bridge_source_truth_table"]["current_accepts_bridge_positive_zbeta_source"],
        "recommended_next": result["recommended_next_honest_step"],
        "out": rel(OUT),
        "md": rel(MD),
    }, indent=2, sort_keys=True, ensure_ascii=False))
