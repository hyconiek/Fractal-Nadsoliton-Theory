"""P3139/S2089: D_HL lane no-new-live-frontier reconciliation.

This is the proof-level reconciliation requested after P3138.  It does not add
another receiver.  Instead it constructs the missing obstruction object for the
D_HL lane: a compact source-obligation matrix plus a group-action theorem
showing why the tested receiver families cannot become an import-free joint
(r, lambda) source on current artifacts.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3139_s2089_dhl_lane_no_new_frontier_reconciliation.json"
MD = GEN / "p3139_s2089_dhl_lane_no_new_frontier_reconciliation.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3133": GEN / "p3133_s2083_legacy_beta_tors_helical_lock_defect_audit.json",
    "P3134": GEN / "p3134_s2084_legacy_phase_torsion_dhl_candidate_audit.json",
    "P3135": GEN / "p3135_s2085_dhl_origin_polarity_selector_source_matrix.json",
    "P3136": GEN / "p3136_s2086_fourier_phase_dhl_joint_source_candidate.json",
    "P3137": GEN / "p3137_s2087_fourier_frame_source_law_audit.json",
    "P3138": GEN / "p3138_s2088_nonfourier_dhl_extremum_joint_source_audit.json",
}

SEARCH_PATTERNS = [
    "D_HL",
    "J_DHL",
    "F_DHL",
    "Zeta_OS",
    "Gamma_SO",
    "QW-2191",
    "no-new-live-frontier",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def count_pattern(pattern: str) -> int:
    # Local grep equivalent using Python keeps this audit self-contained and reproducible.
    count = 0
    for path in ROOT.rglob("*"):
        if path.is_file() and path.suffix in {".py", ".md", ".json", ".tex"} and ".git" not in path.parts:
            try:
                count += path.read_text(encoding="utf-8", errors="ignore").count(pattern)
            except OSError:
                pass
    return count


def stage_row(stage: str, source: dict[str, Any], constructed: bool, receiver: bool, origin: bool, polarity: bool, import_free: bool, variational: bool, blocker: str) -> dict[str, Any]:
    return {
        "stage": stage,
        "input_status": source.get("status"),
        "constructed_or_audited_object": constructed,
        "has_receiver_structure": receiver,
        "selects_absolute_origin": origin,
        "selects_unpaired_polarity": polarity,
        "import_free_source_exported": import_free,
        "variational_or_Ltotal_ready": variational,
        "blocker": blocker,
    }


def build_payload() -> dict[str, Any]:
    inputs = {name: load_json(path) for name, path in INPUTS.items()}
    p3138_cert = inputs["P3138"].get("finite_certificate", {})
    p3138_neg = inputs["P3138"].get("decision", {}).get("negative_export_flags", {})

    rows = [
        stage_row("P3133 legacy beta_tors/cos phase intake", inputs["P3133"], True, True, False, False, False, False, "legacy ingredients are visible but do not supply support representative plus unpaired polarity after the diagonal quotient"),
        stage_row("P3134 explicit D_HL formula", inputs["P3134"], True, True, False, False, False, False, "D_HL shape is nonzero, but origin r and lambda are carried as parameters"),
        stage_row("P3135 joint (r,lambda) matrix", inputs["P3135"], True, False, False, False, False, False, "Z12 x {±1} remains one orbit under translations/Aut/lambda pairing"),
        stage_row("P3136 Fourier phase receiver", inputs["P3136"], True, True, False, False, False, False, "phase extraction imports Fourier frame, selected character/mode, and polarity convention"),
        stage_row("P3137 Fourier frame source law", inputs["P3137"], True, True, False, False, False, False, "primitive character/pair and phase-zero choices remain orbit/gauge choices; zero source candidates pass all gates"),
        stage_row("P3138 non-Fourier extremum receiver", inputs["P3138"], True, True, False, False, False, False, "local extrema/zero-crossing data translate covariantly and require imported cell-order/orientation conventions"),
    ]

    source_obligations = [
        {"obligation": "O1_explicit_object", "satisfied_by_current_lane": True, "evidence": "P3134 constructs D_HL; P3136-P3138 construct receiver/source candidates"},
        {"obligation": "O2_nonzero_inversion_odd_value", "satisfied_by_current_lane": True, "evidence": "P3134/P3136 confirm nonzero local shape/coefficients in scoped profiles"},
        {"obligation": "O3_absolute_support_origin_after_Z12_quotient", "satisfied_by_current_lane": False, "evidence": "P3134/P3138 require carrying or importing r / a cell order"},
        {"obligation": "O4_unpaired_polarity_lambda", "satisfied_by_current_lane": False, "evidence": "P3135 and P3138 retain lambda or positive/negative polarity pairing"},
        {"obligation": "O5_import_free_strict_source_law", "satisfied_by_current_lane": False, "evidence": "P3135-P3138 source gates accept zero candidates"},
        {"obligation": "O6_variational_unit_coupling", "satisfied_by_current_lane": False, "evidence": "equation sheet / Lagrangian draft mark all rows as non-variational receivers"},
    ]

    theorem = {
        "name": "P3139_T1_DHL_receiver_to_source_obstruction",
        "statement": "On current artifacts, every tested D_HL continuation supplies at most receiver data.  A map that is translation-covariant on Z12 receiver cells but has no additional strict pointed datum cannot descend to an import-free absolute origin selector; an inversion-paired polarity class cannot fix lambda without an independent signed source. Therefore the P3133-P3138 D_HL lane exports no J_DHL/D_HL/Zeta_OS/Gamma_SO source and unlocks no selector, bridge, L_total, or ToE closure.",
        "finite_support": {
            "P3138_profiles_tested": p3138_cert.get("profiles_tested"),
            "P3138_translation_covariant_receiver_rows": p3138_cert.get("translation_t1_covariant_receiver_rows"),
            "P3138_accepted_import_free_joint_sources": p3138_cert.get("accepted_import_free_joint_sources"),
            "P3138_negative_flags_all_false": bool(p3138_neg) and all(value is False for value in p3138_neg.values()),
        },
        "proof_steps": [
            "P3134 constructs a real local D_HL shape but parameterizes it by r and lambda.",
            "P3135 proves the joint pair space is not uniquely selected by current invariant data.",
            "P3136-P3137 prove Fourier receivers/frames do not source the missing frame and phase-zero data.",
            "P3138 proves a non-Fourier local extremum receiver is still translation-covariant and polarity-paired.",
            "No current row satisfies the combined obligations O3-O6, so the lane has no proof-grade live frontier without a genuinely new strict source law.",
        ],
    }

    grep_summary = {pattern: count_pattern(pattern) for pattern in SEARCH_PATTERNS}
    accepted_rows = sum(1 for row in rows if row["import_free_source_exported"])
    blocked_obligations = [row["obligation"] for row in source_obligations if not row["satisfied_by_current_lane"]]

    return {
        "status": "P3139_DHL_LANE_NO_NEW_LIVE_FRONTIER_RECONCILIATION",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "repo_grep_summary": grep_summary,
        "constructed_object": {
            "name": "D_HL source-obligation closure matrix and receiver-to-source obstruction theorem",
            "purpose": "reconcile P3133-P3138 without inventing another receiver variant",
        },
        "stage_rows": rows,
        "source_obligations": source_obligations,
        "theorem": theorem,
        "decision": {
            "accepted_import_free_source_rows": accepted_rows,
            "blocked_obligations": blocked_obligations,
            "bounded_result": "P3139 reconciles P3133-P3138 and closes the current D_HL selector/symmetry-breaking lane as no-new-live-frontier on present artifacts.  The lane has real constructed mathematics and real receivers, but the missing theoretical object is still a strict source law for absolute origin and unpaired polarity, not another extraction rule.",
            "negative_export_flags": {
                "new_DHL_live_frontier_unlocked": False,
                "J_DHL_source_exported": False,
                "D_HL_source_exported": False,
                "Zeta_OS_exported": False,
                "Gamma_SO_exported": False,
                "QW_2191_discharged": False,
                "strict_selector_closure_exported": False,
                "bridge_completion_exported": False,
                "legacy_role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "next_honest_step": "Pivot out of the D_HL receiver family.  The next proof-grade move must introduce one genuinely new strict typed object outside Fourier/extremum/orbit-representative/support-order receivers, or else preserve this no-new-live-frontier certificate.  If a new object is supplied, test exactly one theorem: whether it exports a nonimported absolute origin and a nonzero unpaired polarity coupled to D_HL.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P3139/S2089 D_HL lane no-new-live-frontier reconciliation",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed obstruction object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Purpose: {payload['constructed_object']['purpose']}",
        "",
        "## Source-obligation matrix",
    ]
    for row in payload["source_obligations"]:
        lines.append(f"- `{row['obligation']}`: `{row['satisfied_by_current_lane']}` — {row['evidence']}")
    theorem = payload["theorem"]
    lines.extend([
        "",
        "## Theorem",
        f"`{theorem['name']}`: {theorem['statement']}",
        "",
        "## Finite support",
    ])
    for key, value in theorem["finite_support"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3139/S2089 D_HL lane no-new-live-frontier reconciliation", "## P3139/S2089 D_HL lane no-new-live-frontier reconciliation\n\n`P3139/S2089` reconciles P3133-P3138 by constructing a compact source-obligation matrix and receiver-to-source obstruction theorem for the `D_HL` lane.  The lane has real local mathematics (`D_HL`) and real receivers (Fourier phase/frame and non-Fourier extrema), but no current artifact satisfies the combined obligations of absolute support origin after the `Z12` quotient, unpaired `lambda` polarity, import-free strict source law, and variational/unit coupling.  Therefore the current `D_HL` selector/symmetry-breaking lane is no-new-live-frontier unless a genuinely new strict typed object outside these receiver families is supplied.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3139/S2089 D_HL lane remains non-variational", "## P3139/S2089 D_HL lane remains non-variational\n\n`P3139/S2089` closes the current `D_HL` receiver family as no-new-live-frontier for variational purposes.  Without a strict source law for absolute support origin and unpaired polarity, no `D_HL` receiver can define a Lagrangian density, Hamiltonian normalization, spacetime EOM, physical unit, `L_total`, bridge-completion theorem, role-transfer theorem, or ToE closure.\n")
    append_once(AGENTS, "Current D_HL lane no-new-live-frontier reconciliation guardrail (P3139/S2089, 2026-07-13)", "## Current D_HL lane no-new-live-frontier reconciliation guardrail (P3139/S2089, 2026-07-13)\n\n- P3139 reconciles P3133-P3138 by constructing a source-obligation matrix and receiver-to-source obstruction theorem for the `D_HL` selector/symmetry-breaking lane.\n- The lane has real constructed objects and receivers, but no current artifact supplies all missing strict obligations: absolute support origin after the `Z12` quotient, unpaired `lambda` polarity, import-free strict source law, and variational/unit coupling.\n- Do not continue Fourier-frame, local-extremum, orbit-representative, support-order, or receiver-extraction variants as `J_DHL`, `D_HL`, `Zeta_OS`, `Gamma_SO`, `QW-2191` discharge, strict selector closure, bridge completion, role transfer, `L_total`, or ToE closure.\n- The next proof-grade move must pivot to one genuinely new strict typed object outside the `D_HL` receiver family, or preserve this no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
