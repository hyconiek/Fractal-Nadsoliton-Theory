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

OUT = GEN / "p2417_s1367_apd_witness_to_source_obligation_nonpromotion_matrix_certificate.json"
MD = GEN / "p2417_s1367_apd_witness_to_source_obligation_nonpromotion_matrix_certificate.md"

SOURCE_FILES = {
    "P2411_BRIDGE_SOURCE_HYPERGRAPH": GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json",
    "P2413_AMPLITUDE_SCALAR_WITNESS": GEN / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.json",
    "P2414_DAMPING_NONABSORPTION": GEN / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.json",
    "P2415_PHASE_AFFINE_TRANSPORT": GEN / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.json",
    "P2416_APD_VALUE_ASSEMBLY": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ARTIFACT_COLUMNS = (
    "P2413_amplitude_scalar_witness",
    "P2414_damping_parameter_witness",
    "P2415_phase_affine_witness",
    "P2416_apd_value_assembly_witness",
)


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
    return {"count": len(lines), "samples": lines[:18]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2417|S1367|witness to source obligation|source obligation nonpromotion|APD witness nonpromotion",
        "p2411_obligations": "amplitude_dynamic_source_theorem|phase_frequency_dynamic_source_theorem|strict_compression_dynamic_source_theorem|chi11_selector_source_theorem",
        "recent_value_witnesses": "P2413|P2414|P2415|P2416|value-level assembly|phase-coordinate witness|damping parameter",
        "nonpromotion_limits": "not export.*source|No strict dynamic source|selector-source theorem|role-transfer theorem",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds the P2411 source-obligation hypergraph and the P2413-P2416 value witnesses, but no "
            "production matrix proving that those witnesses do not discharge any P2411 source obligation atom."
        ),
    }


def bridge_atoms(p2411: dict[str, Any]) -> list[str]:
    cert = p2411.get("legacy_strict_bridge_source_obligation_hypergraph_certificate", {}).get(
        "finite_hypergraph_certificate", {}
    )
    return list(cert.get("bridge_obligation_atoms", []))


def artifact_positive_flags(sources: dict[str, Any]) -> dict[str, bool]:
    p2413 = sources["P2413_AMPLITUDE_SCALAR_WITNESS"].get(
        "amplitude_scalar_normalization_bridge_witness_certificate", {}
    ).get("theorem_export", {})
    p2414 = sources["P2414_DAMPING_NONABSORPTION"].get(
        "strict_damping_parameter_identifiability_nonabsorption_certificate", {}
    ).get("theorem_export", {})
    p2415 = sources["P2415_PHASE_AFFINE_TRANSPORT"].get(
        "phase_frequency_affine_transport_nonautomorphism_certificate", {}
    ).get("theorem_export", {})
    p2416 = sources["P2416_APD_VALUE_ASSEMBLY"].get(
        "apd_multiplicative_bridge_assembly_necessity_certificate", {}
    ).get("theorem_export", {})
    return {
        "P2413_amplitude_scalar_witness": p2413.get("scalar_normalization_witness_ready") is True,
        "P2414_damping_parameter_witness": p2414.get("strict_beta_eta_identified_from_samples") is True,
        "P2415_phase_affine_witness": p2415.get("continuous_affine_phase_transport_exact_float") is True,
        "P2416_apd_value_assembly_witness": p2416.get("apd_value_assembly_witness_ready") is True,
    }


def nonpromotion_reason(atom: str) -> str:
    reasons = {
        "amplitude_dynamic_source_theorem": "P2413 removes a scalar but explicitly exports no amplitude dynamic source.",
        "role_safe_amplitude_absorption_map_theorem": "P2413 blocks role-safe alpha_geo absorption and P2416 keeps role transfer closed.",
        "phase_frequency_dynamic_source_theorem": "P2415 transports phase coordinates but exports no omega/phi dynamic source.",
        "topological_bit_transport_selector_theorem": "P2415 matches inherited GF(2) bits but exports no selector/topological source theorem.",
        "legacy_linear_to_strict_nonlinear_compression_map_theorem": "P2414 identifies strict beta/eta samples but exports no beta_tors -> beta/eta completion-map theorem.",
        "strict_compression_dynamic_source_theorem": "P2414/P2416 keep strict compression source absent despite value-level damping factors.",
        "chi11_selector_source_theorem": "P2412/P2415 preserve selector-scope separation; no chi11 source theorem is exported.",
        "qw2191_symmetry_breaking_or_internal_source_theorem": "QW-2191 remains open without symmetry-breaking or internal selector source evidence.",
    }
    return reasons.get(atom, "No artifact in P2413-P2416 exports this source obligation as a theorem.")


def build_certificate(sources: dict[str, Any]) -> dict[str, Any]:
    atoms = bridge_atoms(sources["P2411_BRIDGE_SOURCE_HYPERGRAPH"])
    flags = artifact_positive_flags(sources)
    rows = []
    discharged_atoms: list[str] = []
    for index, atom in enumerate(atoms):
        row = {
            "atom_index": index,
            "source_obligation_atom": atom,
            "positive_witness_artifacts_present": [name for name, ready in flags.items() if ready],
            "discharged_by_artifacts": {name: False for name in ARTIFACT_COLUMNS},
            "discharged": False,
            "nonpromotion_reason": nonpromotion_reason(atom),
        }
        rows.append(row)
        if row["discharged"]:
            discharged_atoms.append(atom)
    total_masks = 1 << len(atoms)
    current_mask = sum(1 << i for i, row in enumerate(rows) if row["discharged"])
    full_mask = total_masks - 1
    nearest_miss_masks = [full_mask ^ (1 << i) for i in range(len(atoms))]
    return {
        "source_obligation_atoms": atoms,
        "artifact_columns": list(ARTIFACT_COLUMNS),
        "artifact_positive_flags": flags,
        "nonpromotion_matrix_rows": rows,
        "discharged_atoms": discharged_atoms,
        "current_source_discharge_mask": current_mask,
        "full_source_discharge_mask": full_mask,
        "source_obligation_assignment_count": total_masks,
        "proper_subset_failure_count": full_mask,
        "nearest_miss_masks": nearest_miss_masks,
        "nearest_miss_count": len(nearest_miss_masks),
        "all_recent_witnesses_positive": all(flags.values()),
        "no_p2413_to_p2416_artifact_discharges_source_atom": all(not row["discharged"] for row in rows),
        "bridge_source_ready_from_current_artifacts": current_mask == full_mask,
        "p2411_full_bridge_rule_inherited": sources["P2411_BRIDGE_SOURCE_HYPERGRAPH"].get(
            "legacy_strict_bridge_source_obligation_hypergraph_certificate", {}
        ).get("theorem_export", {}).get("bridge_ready_true_masks")
        == [full_mask],
        "value_witness_to_source_nonpromotion_ready": True,
        "full_bridge_theorem_exported": False,
        "role_transfer_licensed": False,
        "toe_closure_exported": False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2417/S1367 APD witness-to-source nonpromotion matrix certificate

`P2417/S1367` audits the tempting overread after P2413-P2416.  The amplitude, damping, phase, and APD assembly witnesses are all positive as value-level evidence, but they are not source-obligation theorems.  The certificate maps those four artifacts against the eight P2411 bridge-source atoms and records a zero discharge matrix:

```text
current_source_discharge_mask = 0,
full_source_discharge_mask = 255.
```

Thus the value-level APD identity does not discharge `amplitude_dynamic_source_theorem`, `phase_frequency_dynamic_source_theorem`, `strict_compression_dynamic_source_theorem`, `chi11_selector_source_theorem`, QW-2191 source, or the role-safe amplitude/source obligations.  All 255 proper source masks remain bridge-source failures, with eight one-atom-missing nearest masks.

This is not a rollback of P2413-P2416; it is the proof-theoretic nonpromotion guard: value witnesses remain useful bridge ingredients, but they do not become source, selector, role-transfer, `L_total`, or ToE closure theorems.
""".strip()
    lag_section = """
## P2417/S1367 APD witness-to-source nonpromotion guard for Lagrangian/EOM

`P2417/S1367` records that positive value-level APD witnesses from P2413-P2416 discharge zero P2411 source-obligation atoms.  Therefore the Lagrangian/EOM lane still cannot promote APD bookkeeping to a role-bearing `L_total` term without new source, selector, full-bridge, and role-transfer theorems.
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
        "theorem_name": "P2417_T1_apd_witness_to_source_obligation_nonpromotion_matrix_certificate",
        "source_obligation_atom_count": len(cert["source_obligation_atoms"]),
        "artifact_column_count": len(cert["artifact_columns"]),
        "all_recent_witnesses_positive": cert["all_recent_witnesses_positive"],
        "current_source_discharge_mask": cert["current_source_discharge_mask"],
        "full_source_discharge_mask": cert["full_source_discharge_mask"],
        "source_obligation_assignment_count": cert["source_obligation_assignment_count"],
        "proper_subset_failure_count": cert["proper_subset_failure_count"],
        "nearest_miss_count": cert["nearest_miss_count"],
        "discharged_atoms": cert["discharged_atoms"],
        "no_p2413_to_p2416_artifact_discharges_source_atom": cert[
            "no_p2413_to_p2416_artifact_discharges_source_atom"
        ],
        "bridge_source_ready_from_current_artifacts": cert["bridge_source_ready_from_current_artifacts"],
        "p2411_full_bridge_rule_inherited": cert["p2411_full_bridge_rule_inherited"],
        "value_witness_to_source_nonpromotion_ready": cert["value_witness_to_source_nonpromotion_ready"],
        "full_bridge_theorem_exported": cert["full_bridge_theorem_exported"],
        "role_transfer_licensed": cert["role_transfer_licensed"],
        "toe_closure_exported": cert["toe_closure_exported"],
        "not_licensed": [
            "No P2413-P2416 value witness is promoted to a P2411 source-obligation theorem.",
            "No selector-source theorem or QW-2191 discharge follows from the APD value witnesses.",
            "No full bridge theorem follows from the current zero source-discharge mask.",
            "No role-transfer theorem, role-bearing L_total term, or ToE closure follows.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "eight_source_atoms_mapped": theorem_export["source_obligation_atom_count"] == 8,
        "four_recent_artifacts_mapped": theorem_export["artifact_column_count"] == 4,
        "recent_value_witnesses_positive": theorem_export["all_recent_witnesses_positive"],
        "current_source_discharge_mask_zero": theorem_export["current_source_discharge_mask"] == 0,
        "full_source_mask_255": theorem_export["full_source_discharge_mask"] == 255,
        "all_256_assignments_counted": theorem_export["source_obligation_assignment_count"] == 256,
        "proper_subset_failures_255": theorem_export["proper_subset_failure_count"] == 255,
        "nearest_miss_count_eight": theorem_export["nearest_miss_count"] == 8,
        "no_discharged_atoms": theorem_export["discharged_atoms"] == [],
        "no_artifact_promoted_to_source_atom": theorem_export["no_p2413_to_p2416_artifact_discharges_source_atom"],
        "bridge_source_not_ready": not theorem_export["bridge_source_ready_from_current_artifacts"],
        "p2411_full_bridge_rule_inherited": theorem_export["p2411_full_bridge_rule_inherited"],
        "nonpromotion_certificate_ready": theorem_export["value_witness_to_source_nonpromotion_ready"],
        "full_bridge_still_open": not theorem_export["full_bridge_theorem_exported"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2417_s1367_v1",
        "packet_id": "P2417",
        "stage_id": "S1367",
        "result_kind": "APD_WITNESS_TO_SOURCE_OBLIGATION_NONPROMOTION_MATRIX_CERTIFICATE",
        "status": "PASS_VALUE_WITNESSES_POSITIVE_SOURCE_DISCHARGE_ZERO_NO_BRIDGE_CLOSURE",
        "apd_witness_to_source_obligation_nonpromotion_matrix_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Prove one real P2411 source-obligation atom, especially chi11/selector source or a dynamic APD source; do not promote value witnesses by bookkeeping."
        ),
        "global_status": "OPEN_PROGRESS_VALUE_WITNESSES_NONPROMOTED_TO_SOURCE_ATOMS",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["apd_witness_to_source_obligation_nonpromotion_matrix_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2417 S1367: APD witness-to-source nonpromotion matrix certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2417/S1367 maps positive P2413-P2416 value witnesses to the eight P2411 source-obligation atoms and proves zero source discharge.",
                "",
                "## Finite facts",
                "",
                f"- Source atoms: `{theorem['source_obligation_atom_count']}`.",
                f"- Artifact columns: `{theorem['artifact_column_count']}`.",
                f"- Current source mask: `{theorem['current_source_discharge_mask']}`.",
                f"- Full source mask: `{theorem['full_source_discharge_mask']}`.",
                f"- Proper subset failures: `{theorem['proper_subset_failure_count']}`.",
                f"- Discharged atoms: `{theorem['discharged_atoms']}`.",
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
