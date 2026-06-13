#!/usr/bin/env python3
"""P2676/S1626: finite pre-collapse naturality-square audit for S6.

Tests whether the missing Sigma_sel_src_target_v1 -> F301 typed arrow can be
represented by a nonquotient, nonprojector, chart-label-retaining finite
naturality square without importing a convention choice.
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
OUT = GEN / "p2676_s1626_sigma_f301_precollapse_naturality_square_audit.json"
MD = GEN / "p2676_s1626_sigma_f301_precollapse_naturality_square_audit.md"
P2675 = GEN / "p2675_s1625_sigma_to_f301_typed_arrow_s6_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "precollapse_naturality_square_exported",
    "sigma_to_f301_nonconvention_arrow_exported",
    "xor_xnor_orientation_selected_internally",
    "s6_sigma_to_f301_typed_arrow_exported",
    "o3_chart_sensitive_pair12_typed_seed_exported",
    "boundary_square_cycle_typed_arrow_exported",
    "sector_swap_sourced_invariant_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    data = json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(data.encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:60]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "naturality_square_content": "naturality square|commutative square|commutative diagram|natural.*square",
        "sigma_to_f301_content": "Sigma_sel_src_target_v1.*F301|F301.*Sigma_sel_src_target_v1|Sigma_sel_src_target_v1 -> F301|surviving F301",
        "precollapse_content": "pre-collapse|prior to Q_basis|before Q_basis|preLM|basis-free class-reduction|Q_basis_sel_v1 terminal collapse",
        "chart_label_retaining_content": "chart-label-retaining|chart label retaining|chart-sensitive|chart labels",
        "nonprojector_content": "nonprojector|nonprojective|projector-only|without_projector_only_atlas_collapse|local pair1/pair2 atlas",
        "nonconvention_content": "nonconvention|convention|by fiat|must_not_identify|without.*fiat|orientation",
        "selector_source_guard_content": "QW-2191|selector source|source obstruction|sector-swap|sourced invariant|L_total|ToE closure",
    }
    return {
        "tool": "rg",
        "mode": "content-first pre-collapse naturality-square audit for Sigma_sel_src_target_v1 -> F301; not name/number-only",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def upstream_consistency() -> dict[str, Any]:
    p2675 = load_json(P2675)
    decision = p2675.get("closure_decision", {})
    lattice = p2675.get("finite_s6_lattice", {})
    return {
        "p2675_s6_not_passed": decision.get("s6_exported_now") is False,
        "p2675_o3_not_passed": decision.get("o3_exported_now") is False,
        "p2675_missing_pre_collapse": "A5_pre_collapse_nonquotient_descent_exported" in lattice.get("missing_s6_obligations_now", []),
        "p2675_missing_nonprojector": "A6_nonprojector_local_atlas_descent_exported" in lattice.get("missing_s6_obligations_now", []),
        "p2675_missing_no_fiat": "A7_no_fiat_identification_proof_exported" in lattice.get("missing_s6_obligations_now", []),
    }


def truth_value(table: int, sigma: int, chart: int) -> int:
    return (table >> (2 * sigma + chart)) & 1


def finite_boolean_naturality_witness() -> dict[str, Any]:
    """Enumerate all h: Sigma_bit x Chart_bit -> F301_pair_bit maps."""
    rows: list[dict[str, Any]] = []
    for table in range(16):
        values = {(sigma, chart): truth_value(table, sigma, chart) for sigma in (0, 1) for chart in (0, 1)}
        chart_equivariant = all(values[(sigma, 1 - chart)] == 1 - values[(sigma, chart)] for sigma in (0, 1) for chart in (0, 1))
        source_sensitive = any(values[(0, chart)] != values[(1, chart)] for chart in (0, 1))
        chart_sensitive = any(values[(sigma, 0)] != values[(sigma, 1)] for sigma in (0, 1))
        nonconstant = len(set(values.values())) > 1
        passes_formal_square = chart_equivariant and source_sensitive and chart_sensitive and nonconstant
        if passes_formal_square:
            if values[(0, 0)] == 0 and values[(1, 0)] == 1:
                convention_class = "xor_orientation"
            elif values[(0, 0)] == 1 and values[(1, 0)] == 0:
                convention_class = "xnor_reversal_orientation"
            else:
                convention_class = "other"
        else:
            convention_class = "not_formal_square"
        rows.append(
            {
                "table": table,
                "values": {f"sigma{sigma}_chart{chart}": values[(sigma, chart)] for sigma in (0, 1) for chart in (0, 1)},
                "chart_equivariant_under_chart_swap_pair_swap": chart_equivariant,
                "source_sensitive": source_sensitive,
                "chart_sensitive_nonprojector_form": chart_sensitive,
                "nonconstant": nonconstant,
                "passes_formal_square": passes_formal_square,
                "convention_class": convention_class,
                "has_internal_orientation_source": False,
                "passes_export_gate": False,
            }
        )
    formal = [row for row in rows if row["passes_formal_square"]]
    return {
        "domain": "all Boolean maps h(sigma_selector_bit, chart_label_bit) -> f301_pair_bit",
        "total_maps_checked": len(rows),
        "formal_square_candidate_count": len(formal),
        "formal_square_tables": [row["table"] for row in formal],
        "formal_square_convention_classes": sorted({row["convention_class"] for row in formal}),
        "passing_export_gate_count": sum(row["passes_export_gate"] for row in rows),
        "rows": rows,
    }


def obstruction_certificates(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "certificate": "finite_formal_candidates_exist",
            "value": witness["formal_square_candidate_count"] > 0,
            "meaning": "Boolean naturality-like maps exist as formal XOR/XNOR orientation choices.",
        },
        {
            "certificate": "orientation_reversal_pair_survives",
            "value": set(witness["formal_square_convention_classes"]) == {"xnor_reversal_orientation", "xor_orientation"},
            "meaning": "The formal candidates appear in a reversal pair, so finite naturality alone does not pick an internal orientation.",
        },
        {
            "certificate": "internal_orientation_source_absent",
            "value": True,
            "meaning": "No audited row has an internal selector/source bit that chooses XOR over XNOR without convention.",
        },
        {
            "certificate": "export_gate_zero",
            "value": witness["passing_export_gate_count"] == 0,
            "meaning": "No formal Boolean map is promoted to a current Sigma->F301 typed-arrow export.",
        },
    ]


def closure_decision(witness: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2676_PRECOLLAPSE_NATURALITY_SQUARE_AUDIT__FORMAL_XOR_XNOR_PAIR_EXISTS_BUT_NO_INTERNAL_ORIENTATION_SOURCE",
        "professorial_verdict": (
            "P2676 builds the requested finite naturality-square witness for S6. Exhausting all 16 Boolean maps finds two formal chart-equivariant, source-sensitive, nonprojector candidates, "
            "but they are exactly an XOR/XNOR reversal pair. Since the repo does not export an internal orientation/source selecting one member of that pair, the passing export gate remains zero. "
            "Thus the naturality-square route does not currently supply the missing Sigma->F301 typed arrow."
        ),
        "next_honest_step": (
            "The next honest move is not another same-shape Boolean/naturality lift. Either provide a new internal orientation source that breaks the XOR/XNOR reversal inside the Sigma->F301 square, "
            "or record a no-go for the current S6/O3 typed-seed route and return to a different bridge-completion/source class. O4/O5 and L_total reopening remain blocked."
        ),
        "formal_square_candidate_count": witness["formal_square_candidate_count"],
        "passing_export_gate_count": witness["passing_export_gate_count"],
        "precollapse_naturality_square_exported_now": False,
        "sigma_to_f301_typed_arrow_exported_now": False,
        "s6_exported_now": False,
        "o3_exported_now": False,
        "boundary_square_arrow_allowed_next": False,
        "sector_swap_invariant_allowed_next": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["finite_boolean_naturality_witness"]
    decision = payload["closure_decision"]
    lines = [
        "# P2676/S1626 Sigma->F301 pre-collapse naturality-square audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(
        [
            "",
            "## Finite Boolean naturality witness",
            f"Total maps checked: `{witness['total_maps_checked']}`.",
            f"Formal square candidates: `{witness['formal_square_candidate_count']}`.",
            f"Formal candidate tables: `{witness['formal_square_tables']}`.",
            f"Formal convention classes: `{witness['formal_square_convention_classes']}`.",
            f"Passing export gate count: `{witness['passing_export_gate_count']}`.",
            "",
            "## Obstruction certificates",
        ]
    )
    for row in payload["obstruction_certificates"]:
        lines.append(f"- `{row['certificate']}`: `{row['value']}` — {row['meaning']}")
    lines.extend(
        [
            "",
            "## Verdict",
            decision["professorial_verdict"],
            f"Decision: `{decision['decision']}`.",
            "",
            "## Next honest step",
            decision["next_honest_step"],
            "",
            "## Negative exports",
        ]
    )
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    witness = finite_boolean_naturality_witness()
    payload: dict[str, Any] = {
        "status": "P2676_PRECOLLAPSE_NATURALITY_SQUARE_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2675": sha256_file(P2675),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "finite_boolean_naturality_witness": witness,
        "obstruction_certificates": obstruction_certificates(witness),
        "closure_decision": closure_decision(witness),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2676/S1626 Sigma-F301 naturality-square guard",
        "## P2676/S1626 Sigma-F301 naturality-square guard\n\n"
        "`P2676/S1626` exhausts all `16` Boolean maps `h(sigma_selector_bit, chart_label_bit) -> f301_pair_bit` as a finite pre-collapse naturality-square witness for S6.  Two formal chart-equivariant, source-sensitive, nonprojector candidates exist, but they are exactly the XOR/XNOR orientation-reversal pair, and no internal selector/source chooses between them.  Therefore the export gate remains zero: no pre-collapse naturality square, no nonconvention `Sigma_sel_src_target_v1 -> F301` typed arrow, no S6/O3 pass, no O4/O5 admission, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2676/S1626 Sigma-F301 naturality Ltotal guard",
        "## P2676/S1626 Sigma-F301 naturality Ltotal guard\n\n"
        "`P2676/S1626` keeps `L_total` closed after the finite naturality-square audit.  Formal XOR/XNOR Sigma->F301 square candidates do not become variational source terms unless a strict internal orientation source selects one member of the reversal pair and exports an actual pre-collapse typed arrow.\n",
    )
    return payload


if __name__ == "__main__":
    main()
