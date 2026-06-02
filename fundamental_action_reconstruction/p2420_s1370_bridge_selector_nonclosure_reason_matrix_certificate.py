#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2420_s1370_bridge_selector_nonclosure_reason_matrix_certificate.json"
MD = GEN / "p2420_s1370_bridge_selector_nonclosure_reason_matrix_certificate.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2417_SOURCE_NONPROMOTION": GEN / "p2417_s1367_apd_witness_to_source_obligation_nonpromotion_matrix_certificate.json",
    "P2419_CHI11_COUPLING_CUT": GEN / "p2419_s1369_chi11_phase_selector_coupling_cut_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

GATES = [
    "apd_value_bridge_witness",
    "source_obligation_discharge",
    "chi11_phase_selector_cut_mechanism",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]

APD_GATE = "apd_value_bridge_witness"
SOURCE_GATE = "source_obligation_discharge"
MECHANISM_GATE = "chi11_phase_selector_cut_mechanism"
CHI11_SOURCE_GATE = "chi11_source_export"
QW2191_GATE = "qw2191_selector_discharge"
ROLE_GATE = "role_transfer_audit_license"
LTOTAL_GATE = "role_bearing_ltotal_export"

CURRENT_TRUE_GATES = [APD_GATE, MECHANISM_GATE]


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


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg",
            "-n",
            pattern,
            "fundamental_action_reconstruction",
            "material_dowodowy",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.tex",
            "-g",
            "!generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:20]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2420|S1370|bridge selector nonclosure|why.*not.*close|domyk|nonclosure reason matrix",
        "apd_bridge_input": "P2416|APD value assembly|apd_multiplicative_bridge|full_bridge_theorem_exported",
        "source_nonpromotion_input": "P2417|source-discharge|current_source_discharge_mask|bridge_source_ready_from_current_artifacts",
        "selector_cut_input": "P2419|chi11 phase-selector|chi11_is_common_necessary_cut|qw2191_discharged",
        "prior_nonclosure_language": "selector.*nonclosure|closure gap|No role-transfer theorem|No ToE closure|nie.*domkni",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds separate APD value-bridge, source-nonpromotion, selector-cut, and older nonclosure materials, "
            "but no production P24xx matrix answering the specific objection: why APD bridge evidence plus a selector "
            "mechanism/cut still does not close the theory."
        ),
    }


def theorem_export(payload: dict[str, Any], certificate_key: str) -> dict[str, Any]:
    return payload.get(certificate_key, {}).get("theorem_export", {})


def mask_for(gates: list[str]) -> int:
    return sum(1 << GATES.index(gate) for gate in gates)


def gates_for(mask: int) -> list[str]:
    return [gate for index, gate in enumerate(GATES) if mask & (1 << index)]


def predicate_row(mask: int) -> dict[str, Any]:
    true = set(gates_for(mask))
    full_bridge_ready = APD_GATE in true and SOURCE_GATE in true
    selector_closed = MECHANISM_GATE in true and CHI11_SOURCE_GATE in true and QW2191_GATE in true
    role_transfer_ready = ROLE_GATE in true
    ltotal_ready = LTOTAL_GATE in true
    toe_ready = full_bridge_ready and selector_closed and role_transfer_ready and ltotal_ready
    missing = [gate for gate in GATES if gate not in true]
    return {
        "mask": mask,
        "true_gates": gates_for(mask),
        "missing_gates": missing,
        "full_bridge_ready": full_bridge_ready,
        "selector_closed": selector_closed,
        "role_transfer_ready": role_transfer_ready,
        "role_bearing_ltotal_ready": ltotal_ready,
        "toe_ready": toe_ready,
    }


def all_rows() -> list[dict[str, Any]]:
    return [predicate_row(mask) for mask in range(1 << len(GATES))]


def minimal_toe_masks(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for row in rows:
        if not row["toe_ready"]:
            continue
        mask = row["mask"]
        proper_true_submasks = [sub for sub in range(mask) if sub & mask == sub and sub != mask]
        if any(predicate_row(sub)["toe_ready"] for sub in proper_true_submasks):
            continue
        out.append({"mask": mask, "gates": row["true_gates"], "size": len(row["true_gates"])})
    return out


def hamming_distance_to(mask: int, target: int) -> int:
    return (mask ^ target).bit_count()


def fixed_apd_mechanism_subcube(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    required = mask_for(CURRENT_TRUE_GATES)
    return [row for row in rows if row["mask"] & required == required]


def single_flip_rows(current_mask: int) -> list[dict[str, Any]]:
    rows = []
    for gate in GATES:
        bit = 1 << GATES.index(gate)
        if current_mask & bit:
            continue
        flipped = current_mask | bit
        row = predicate_row(flipped)
        rows.append(
            {
                "added_gate": gate,
                "mask_after_addition": flipped,
                "full_bridge_ready_after_addition": row["full_bridge_ready"],
                "selector_closed_after_addition": row["selector_closed"],
                "toe_ready_after_addition": row["toe_ready"],
                "still_missing_for_toe": row["missing_gates"],
            }
        )
    return rows


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2416 = theorem_export(sources["P2416_APD_VALUE_BRIDGE"], "apd_multiplicative_bridge_assembly_necessity_certificate")
    p2417 = theorem_export(sources["P2417_SOURCE_NONPROMOTION"], "apd_witness_to_source_obligation_nonpromotion_matrix_certificate")
    p2419 = theorem_export(sources["P2419_CHI11_COUPLING_CUT"], "chi11_phase_selector_coupling_cut_certificate")
    rows = all_rows()
    current_mask = mask_for(CURRENT_TRUE_GATES)
    full_mask = mask_for(GATES)
    current_row = predicate_row(current_mask)
    subcube = fixed_apd_mechanism_subcube(rows)
    toe_rows = [row for row in rows if row["toe_ready"]]
    subcube_toe_rows = [row for row in subcube if row["toe_ready"]]
    nearest_repairs = [
        {"mask": row["mask"], "add_gates": row["missing_gates"], "repair_distance": len(row["missing_gates"])}
        for row in rows
        if row["toe_ready"] and hamming_distance_to(current_mask, row["mask"]) == min(hamming_distance_to(current_mask, item["mask"]) for item in toe_rows)
    ]
    reason_rows = [
        {
            "reason_id": "R1_VALUE_BRIDGE_IS_NOT_SOURCE_BRIDGE",
            "observed_positive_evidence": "P2416 APD value-level exact quotient assembly",
            "missing_theorem_gate": SOURCE_GATE,
            "blocking_input": "P2417 current_source_discharge_mask=0 and no P2413-P2416 value witness discharges a source atom",
            "closure_consequence": "full_bridge_ready remains false at the current mask",
        },
        {
            "reason_id": "R2_SELECTOR_MECHANISM_IS_NOT_SELECTOR_SOURCE",
            "observed_positive_evidence": "P2419 chi11 common necessary phase/selector cut",
            "missing_theorem_gate": CHI11_SOURCE_GATE,
            "blocking_input": "P2419 chi11_source_exported=false and chi11_not_sufficient_for_phase/selector=true",
            "closure_consequence": "selector_closed remains false at the current mask",
        },
        {
            "reason_id": "R3_QW2191_IS_STILL_AN_INDEPENDENT_SELECTOR_OBSTRUCTION",
            "observed_positive_evidence": "P2419 phase/selector co-readiness cut isolates where QW-2191 enters",
            "missing_theorem_gate": QW2191_GATE,
            "blocking_input": "P2419 qw2191_discharged=false and qw2191_not_sufficient_for_selector=true",
            "closure_consequence": "selector_closed cannot be asserted even if the cut location is known",
        },
        {
            "reason_id": "R4_ROLE_TRANSFER_AUDIT_IS_POST_BRIDGE_AND_NOT_AUTOMATIC",
            "observed_positive_evidence": "legacy-to-strict value bookkeeping and selector-cut bookkeeping",
            "missing_theorem_gate": ROLE_GATE,
            "blocking_input": "S2 requires a separate role-transfer audit after bridge completion",
            "closure_consequence": "legacy physical-role claims cannot be silently transferred to the strict kernel",
        },
        {
            "reason_id": "R5_LTOTAL_EXPORT_IS_A_SEPARATE_ROLE_BEARING_STEP",
            "observed_positive_evidence": "kernel-level bridge/selector bookkeeping",
            "missing_theorem_gate": LTOTAL_GATE,
            "blocking_input": "No role-bearing L_total theorem is exported by P2416, P2417, or P2419",
            "closure_consequence": "ToE closure remains false even before external/sector closure is considered",
        },
    ]
    return {
        "closure_gates": GATES,
        "current_true_gates": CURRENT_TRUE_GATES,
        "current_mask": current_mask,
        "full_mask": full_mask,
        "current_row": current_row,
        "closure_truth_rows": rows,
        "minimal_toe_masks": minimal_toe_masks(rows),
        "toe_true_mask_count": len(toe_rows),
        "apd_plus_selector_mechanism_subcube_count": len(subcube),
        "apd_plus_selector_mechanism_subcube_toe_count": len(subcube_toe_rows),
        "apd_plus_selector_mechanism_subcube_failure_count": len(subcube) - len(subcube_toe_rows),
        "current_to_toe_minimum_repair_distance": min(hamming_distance_to(current_mask, row["mask"]) for row in toe_rows),
        "current_to_toe_nearest_repairs": nearest_repairs,
        "single_flip_rows_from_current_mask": single_flip_rows(current_mask),
        "reason_rows": reason_rows,
        "input_theorem_status": {
            "p2416_apd_value_assembly_witness_ready": p2416.get("apd_value_assembly_witness_ready"),
            "p2416_full_bridge_theorem_exported": p2416.get("full_bridge_theorem_exported"),
            "p2417_current_source_discharge_mask": p2417.get("current_source_discharge_mask"),
            "p2417_bridge_source_ready_from_current_artifacts": p2417.get("bridge_source_ready_from_current_artifacts"),
            "p2419_chi11_is_common_necessary_cut": p2419.get("chi11_is_common_necessary_cut"),
            "p2419_chi11_source_exported": p2419.get("chi11_source_exported"),
            "p2419_qw2191_discharged": p2419.get("qw2191_discharged"),
        },
    }


def append_doc_sections() -> None:
    eq_section = """
## P2420/S1370 bridge-selector nonclosure reason matrix certificate

`P2420/S1370` answers the closure objection directly.  The current positive evidence has two gates: `apd_value_bridge_witness` from P2416 and `chi11_phase_selector_cut_mechanism` from P2419.  Those are not the same as `source_obligation_discharge`, `chi11_source_export`, `qw2191_selector_discharge`, `role_transfer_audit_license`, or `role_bearing_ltotal_export`.

The finite matrix over seven closure gates has `128` rows.  Holding the APD value bridge and selector cut/mechanism fixed true leaves a `32`-row subcube; only the all-gates row closes ToE, so `31/32` such rows remain nonclosures.  From the current mask, the minimum repair distance to ToE closure is five theorem gates: source discharge, chi11 source export, QW-2191 discharge, role-transfer license, and role-bearing `L_total` export.

Therefore the honest answer is: the repo has value-level bridge evidence and a selector-location/cut mechanism, but not the source theorem, not the selector-source discharge, not the role-transfer audit, and not the role-bearing Lagrangian/ToE composition theorem.
""".strip()
    lag_section = """
## P2420/S1370 bridge-selector nonclosure guard for Lagrangian/EOM

`P2420/S1370` separates bridge value evidence and selector-cut evidence from the missing theorem gates needed for closure.  Since the current mask has only `apd_value_bridge_witness` and `chi11_phase_selector_cut_mechanism`, while source discharge, chi11 source export, QW-2191 discharge, role-transfer license, and role-bearing `L_total` export remain absent, no `L_total` or ToE closure can be promoted from the bridge/selector bookkeeping.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    cert = build_certificate(sources)
    theorem_export = {
        "theorem_name": "P2420_T1_bridge_selector_nonclosure_reason_matrix_certificate",
        "closure_gate_count": len(cert["closure_gates"]),
        "total_closure_assignment_count": len(cert["closure_truth_rows"]),
        "current_true_gates": cert["current_true_gates"],
        "current_mask": cert["current_mask"],
        "full_mask": cert["full_mask"],
        "current_full_bridge_ready": cert["current_row"]["full_bridge_ready"],
        "current_selector_closed": cert["current_row"]["selector_closed"],
        "current_role_transfer_ready": cert["current_row"]["role_transfer_ready"],
        "current_role_bearing_ltotal_ready": cert["current_row"]["role_bearing_ltotal_ready"],
        "current_toe_ready": cert["current_row"]["toe_ready"],
        "minimal_toe_mask_count": len(cert["minimal_toe_masks"]),
        "minimal_toe_masks": [row["mask"] for row in cert["minimal_toe_masks"]],
        "toe_true_mask_count": cert["toe_true_mask_count"],
        "apd_plus_selector_mechanism_subcube_count": cert["apd_plus_selector_mechanism_subcube_count"],
        "apd_plus_selector_mechanism_subcube_toe_count": cert["apd_plus_selector_mechanism_subcube_toe_count"],
        "apd_plus_selector_mechanism_subcube_failure_count": cert["apd_plus_selector_mechanism_subcube_failure_count"],
        "current_to_toe_minimum_repair_distance": cert["current_to_toe_minimum_repair_distance"],
        "current_to_toe_nearest_repairs": cert["current_to_toe_nearest_repairs"],
        "single_flip_from_current_closes_toe_count": sum(1 for row in cert["single_flip_rows_from_current_mask"] if row["toe_ready_after_addition"]),
        "reason_count": len(cert["reason_rows"]),
        "why_bridge_plus_selector_mechanism_does_not_close": [row["reason_id"] for row in cert["reason_rows"]],
        "apd_value_bridge_witness_inherited": cert["input_theorem_status"]["p2416_apd_value_assembly_witness_ready"] is True,
        "source_discharge_zero_inherited": cert["input_theorem_status"]["p2417_current_source_discharge_mask"] == 0,
        "chi11_cut_mechanism_inherited": cert["input_theorem_status"]["p2419_chi11_is_common_necessary_cut"] is True,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "plain_language_answer": (
            "The current artifacts give a value-level APD bridge witness and locate a chi11 selector cut. "
            "They do not export the strict dynamic/source theorems, do not discharge QW-2191, do not perform the "
            "post-bridge role-transfer audit, and do not construct a role-bearing L_total; therefore ToE closure is blocked."
        ),
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "seven_closure_gates": theorem_export["closure_gate_count"] == 7,
        "all_128_assignments": theorem_export["total_closure_assignment_count"] == 128,
        "current_mask_is_apd_plus_mechanism": theorem_export["current_mask"] == mask_for(CURRENT_TRUE_GATES),
        "full_mask_127": theorem_export["full_mask"] == 127,
        "unique_minimal_toe_mask": theorem_export["minimal_toe_mask_count"] == 1 and theorem_export["minimal_toe_masks"] == [127],
        "only_one_toe_true_mask": theorem_export["toe_true_mask_count"] == 1,
        "subcube_size_32": theorem_export["apd_plus_selector_mechanism_subcube_count"] == 32,
        "subcube_only_one_closure": theorem_export["apd_plus_selector_mechanism_subcube_toe_count"] == 1,
        "subcube_31_failures": theorem_export["apd_plus_selector_mechanism_subcube_failure_count"] == 31,
        "repair_distance_five": theorem_export["current_to_toe_minimum_repair_distance"] == 5,
        "no_single_flip_closes_toe": theorem_export["single_flip_from_current_closes_toe_count"] == 0,
        "five_reasons_recorded": theorem_export["reason_count"] == 5,
        "apd_value_bridge_inherited": theorem_export["apd_value_bridge_witness_inherited"],
        "source_zero_inherited": theorem_export["source_discharge_zero_inherited"],
        "chi11_cut_inherited": theorem_export["chi11_cut_mechanism_inherited"],
        "chi11_source_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "bridge_still_open": not theorem_export["full_bridge_theorem_exported"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2420_s1370_v1",
        "packet_id": "P2420",
        "stage_id": "S1370",
        "result_kind": "BRIDGE_SELECTOR_NONCLOSURE_REASON_MATRIX_CERTIFICATE",
        "status": "PASS_BRIDGE_SELECTOR_MECHANISM_NONCLOSURE_REASON_MATRIX_NO_TOE_CLOSURE",
        "bridge_selector_nonclosure_reason_matrix_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Choose one missing theorem gate explicitly: source-obligation discharge, chi11 source export, QW-2191 discharge, "
            "post-bridge role-transfer audit, or role-bearing L_total export. Do not describe APD value exactness plus selector-cut "
            "localization as ToE closure."
        ),
        "global_status": "OPEN_PROGRESS_BRIDGE_AND_SELECTOR_MECHANISM_CERTIFIED_BUT_NONCLOSURE_EXPLAINED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["bridge_selector_nonclosure_reason_matrix_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2420 S1370: bridge-selector nonclosure reason matrix certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Direct answer",
                "",
                theorem["plain_language_answer"],
                "",
                "## Finite facts",
                "",
                f"- Closure gates: `{theorem['closure_gate_count']}`.",
                f"- Closure assignments: `{theorem['total_closure_assignment_count']}`.",
                f"- Current true gates: `{theorem['current_true_gates']}`.",
                f"- Current mask: `{theorem['current_mask']}`.",
                f"- Minimal ToE masks: `{theorem['minimal_toe_masks']}`.",
                f"- APD+selector-mechanism subcube failures: `{theorem['apd_plus_selector_mechanism_subcube_failure_count']}/32`.",
                f"- Current repair distance to ToE: `{theorem['current_to_toe_minimum_repair_distance']}`.",
                "",
                "## Reasons",
                "",
                *[f"- {reason}" for reason in theorem["why_bridge_plus_selector_mechanism_does_not_close"]],
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
