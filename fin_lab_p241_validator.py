#!/usr/bin/env python3
"""FIN Laboratory Program 241: blind-custody bundle validator.

The validator checks the eleven-field Program-228/241 transfer contract for
either an ordered twelve-state process record or an event-level double-slit
record.  It can validate hashes, schema, chronology, event columns, controls,
role separation, and a detached registrar signature.

It cannot create an independent empirical record, verify that named humans or
institutions are genuine, or turn a rendered interference image into event
data.  Production admission therefore additionally requires a detached GPG
signature verified against a separately supplied trusted registrar keyring.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys


SCHEMA_VERSION = "FIN-P241-1.0"
ACCEPTED_RESOURCE_TYPES = {"heat_process", "double_slit"}
DOUBLE_SLIT_CONFIGS = {"both_open", "left_only", "right_only", "both_closed"}
PROCESS_REQUIRED_COLUMNS = {
    "event_id",
    "timestamp_utc",
    "run_id",
    "subset",
    "preparation_id",
    "outcome_id",
    "evolution_multiple",
    "intervention",
}
DOUBLE_SLIT_REQUIRED_COLUMNS = {
    "event_id",
    "timestamp_utc",
    "run_id",
    "subset",
    "configuration",
    "x_detector",
    "y_detector",
    "intervention",
}
REQUIRED_FILES = {
    "events.csv",
    "calibration.json",
    "controls.json",
    "environment.json",
    "preregistration.json",
}


@dataclass
class Check:
    number: int
    name: str
    passed: bool
    detail: str

    def as_dict(self) -> dict:
        return {
            "number": self.number,
            "name": self.name,
            "passed": self.passed,
            "detail": self.detail,
        }


def canonical_digest(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_utc(value: str) -> datetime:
    if value.endswith("Z"):
        value = value[:-1] + "+00:00"
    parsed = datetime.fromisoformat(value)
    if parsed.tzinfo is None:
        raise ValueError("timestamp must contain a UTC offset")
    return parsed.astimezone(timezone.utc)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path.name} must contain a JSON object")
    return value


def safe_bundle_path(bundle: Path, relative: str) -> Path:
    candidate = (bundle / relative).resolve()
    root = bundle.resolve()
    if root != candidate and root not in candidate.parents:
        raise ValueError(f"path escapes bundle: {relative}")
    if candidate.is_symlink():
        raise ValueError(f"symlinks are forbidden in evidence bundles: {relative}")
    return candidate


def verify_gpg_signature(
    manifest_path: Path, signature_path: Path | None, keyring_path: Path | None
) -> tuple[bool, str]:
    if signature_path is None or keyring_path is None:
        return False, "production admission requires --signature and --trusted-keyring"
    gpgv = shutil.which("gpgv")
    if gpgv is None:
        return False, "gpgv is unavailable"
    if not signature_path.is_file() or not keyring_path.is_file():
        return False, "signature or trusted keyring file is missing"
    result = subprocess.run(
        [
            gpgv,
            "--keyring",
            str(keyring_path.resolve()),
            str(signature_path.resolve()),
            str(manifest_path.resolve()),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    detail = (result.stderr or result.stdout).strip()
    return result.returncode == 0, detail[-1200:]


def inspect_events(path: Path, resource_type: str) -> dict:
    required = (
        PROCESS_REQUIRED_COLUMNS
        if resource_type == "heat_process"
        else DOUBLE_SLIT_REQUIRED_COLUMNS
    )
    row_count = 0
    columns: set[str] = set()
    preparations: set[int] = set()
    outcomes: set[int] = set()
    multiples: set[int] = set()
    configurations: set[str] = set()
    subsets: set[str] = set()
    run_ids: set[str] = set()
    interventions: set[str] = set()
    minimum_timestamp: datetime | None = None
    maximum_timestamp: datetime | None = None
    count_cells: dict[str, int] = {}

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        missing = required - columns
        if missing:
            return {
                "schema_ok": False,
                "error": f"missing event columns: {sorted(missing)}",
                "columns": sorted(columns),
                "row_count": 0,
            }
        for row in reader:
            row_count += 1
            run_ids.add(row["run_id"].strip())
            subsets.add(row["subset"].strip())
            interventions.add(row["intervention"].strip())
            timestamp = parse_utc(row["timestamp_utc"].strip())
            minimum_timestamp = (
                timestamp
                if minimum_timestamp is None
                else min(minimum_timestamp, timestamp)
            )
            maximum_timestamp = (
                timestamp
                if maximum_timestamp is None
                else max(maximum_timestamp, timestamp)
            )
            if resource_type == "heat_process":
                preparation = int(row["preparation_id"])
                outcome = int(row["outcome_id"])
                multiple = int(row["evolution_multiple"])
                if not 0 <= preparation < 12 or not 0 <= outcome < 12:
                    raise ValueError("process preparation/outcome must lie in 0,...,11")
                if multiple not in {1, 2}:
                    raise ValueError("evolution_multiple must be 1 or 2")
                preparations.add(preparation)
                outcomes.add(outcome)
                multiples.add(multiple)
                key = f"{row['subset']}|{multiple}|{preparation}"
                count_cells[key] = count_cells.get(key, 0) + 1
            else:
                configuration = row["configuration"].strip()
                configurations.add(configuration)
                x_value = float(row["x_detector"])
                y_value = float(row["y_detector"])
                if not math.isfinite(x_value) or not math.isfinite(y_value):
                    raise ValueError("detector coordinates must be finite")

    return {
        "schema_ok": True,
        "columns": sorted(columns),
        "row_count": row_count,
        "preparations": sorted(preparations),
        "outcomes": sorted(outcomes),
        "evolution_multiples": sorted(multiples),
        "configurations": sorted(configurations),
        "subsets": sorted(subsets),
        "run_ids": sorted(run_ids),
        "interventions": sorted(interventions),
        "minimum_timestamp_utc": (
            minimum_timestamp.isoformat() if minimum_timestamp else None
        ),
        "maximum_timestamp_utc": (
            maximum_timestamp.isoformat() if maximum_timestamp else None
        ),
        "count_cells": count_cells,
    }


def bundle_digest(bundle: Path, file_hashes: dict[str, str]) -> str:
    normalized = "\n".join(
        f"{relative}\0{file_hashes[relative]}" for relative in sorted(file_hashes)
    ).encode("utf-8")
    return hashlib.sha256(normalized).hexdigest()


def validate_bundle(
    bundle: Path,
    signature_path: Path | None = None,
    keyring_path: Path | None = None,
) -> dict:
    bundle = bundle.resolve()
    manifest_path = bundle / "bundle_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError("bundle_manifest.json is missing")
    manifest = load_json(manifest_path)
    resource_type = manifest.get("resource_type")
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"schema_version must be {SCHEMA_VERSION}")
    if resource_type not in ACCEPTED_RESOURCE_TYPES:
        raise ValueError(f"resource_type must be one of {sorted(ACCEPTED_RESOURCE_TYPES)}")

    declared_files = manifest.get("files", {})
    if not isinstance(declared_files, dict):
        raise ValueError("manifest.files must be a path-to-SHA256 object")
    actual_hashes: dict[str, str] = {}
    file_hash_errors: list[str] = []
    for relative, expected in declared_files.items():
        path = safe_bundle_path(bundle, relative)
        if not path.is_file():
            file_hash_errors.append(f"{relative}: missing")
            continue
        actual = sha256_file(path)
        actual_hashes[relative] = actual
        if actual.lower() != str(expected).lower():
            file_hash_errors.append(f"{relative}: SHA256 mismatch")
    missing_required_files = REQUIRED_FILES - set(declared_files)

    events_path = bundle / "events.csv"
    event_info = (
        inspect_events(events_path, resource_type)
        if events_path.is_file()
        else {"schema_ok": False, "error": "events.csv missing", "row_count": 0}
    )
    calibration = (
        load_json(bundle / "calibration.json")
        if (bundle / "calibration.json").is_file()
        else {}
    )
    controls = (
        load_json(bundle / "controls.json")
        if (bundle / "controls.json").is_file()
        else {}
    )
    environment = (
        load_json(bundle / "environment.json")
        if (bundle / "environment.json").is_file()
        else {}
    )
    preregistration = (
        load_json(bundle / "preregistration.json")
        if (bundle / "preregistration.json").is_file()
        else {}
    )

    provider = manifest.get("provider", {})
    custody = manifest.get("custody", {})
    declaration = manifest.get("declaration", {})
    holdout = manifest.get("holdout", {})

    external_declaration = (
        manifest.get("evidence_status") == "external"
        and declaration.get("synthetic") is False
        and declaration.get("unit_test_fixture") is False
        and declaration.get("rendered_image_only") is False
    )
    provider_ok = all(
        str(provider.get(key, "")).strip()
        for key in ("name", "institution", "contact", "license")
    ) and external_declaration

    hash_ok = (
        not file_hash_errors
        and not missing_required_files
        and set(actual_hashes) == set(declared_files)
    )
    raw_record_ok = bool(event_info.get("schema_ok")) and event_info.get(
        "row_count", 0
    ) > 0

    if resource_type == "heat_process":
        preparation_ok = event_info.get("preparations") == list(range(12))
        intervention_ok = (
            set(event_info.get("evolution_multiples", [])) == {1, 2}
            and bool(set(event_info.get("interventions", [])) - {""})
        )
    else:
        preparation_ok = set(event_info.get("configurations", [])) == DOUBLE_SLIT_CONFIGS
        intervention_ok = bool(set(event_info.get("interventions", [])) - {""})

    clock_ok = all(
        key in calibration
        for key in (
            "clock_model",
            "clock_unit",
            "clock_uncertainty",
            "time_origin",
        )
    )
    detector_ok = all(
        key in calibration
        for key in (
            "detector_geometry",
            "efficiency_calibration",
            "dark_count_calibration",
            "blur_or_psf_calibration",
        )
    )
    environment_ok = all(
        key in environment
        for key in (
            "background_record",
            "temperature_or_not_applicable",
            "source_stability_record",
            "apparatus_state_record",
        )
    )
    controls_ok = bool(controls.get("negative_controls"))
    if resource_type == "double_slit":
        controls_ok = controls_ok and "both_closed" in controls.get(
            "negative_controls", []
        )

    chronology_ok = False
    chronology_detail = "missing preregistration/acquisition timestamps"
    try:
        preregistered_at = parse_utc(preregistration["registered_at_utc"])
        acquired_at = parse_utc(manifest["acquisition_started_at_utc"])
        chronology_ok = (
            preregistered_at < acquired_at
            and holdout.get("sealed_before_analysis") is True
            and holdout.get("committed_before_analysis") is True
            and bool(preregistration.get("held_out_target"))
            and bool(preregistration.get("frozen_analysis_hash"))
        )
        chronology_detail = (
            f"preregistered={preregistered_at.isoformat()}, "
            f"acquisition={acquired_at.isoformat()}"
        )
    except (KeyError, TypeError, ValueError) as error:
        chronology_detail = str(error)

    roles = [
        str(custody.get("provider_id", "")).strip(),
        str(custody.get("registrar_id", "")).strip(),
        str(custody.get("analyst_id", "")).strip(),
    ]
    role_metadata_ok = all(roles) and len(set(roles)) == 3
    signature_ok, signature_detail = verify_gpg_signature(
        manifest_path, signature_path, keyring_path
    )
    role_separation_ok = role_metadata_ok and signature_ok

    checks = [
        Check(1, "independent provider identity and license", provider_ok, str(provider)),
        Check(
            2,
            "immutable raw-event hash",
            hash_ok,
            "; ".join(file_hash_errors)
            or (
                f"all declared hashes match; bundle digest "
                f"{bundle_digest(bundle, actual_hashes)}"
            ),
        ),
        Check(3, "ordered events or detection coordinates", raw_record_ok, str(event_info)),
        Check(
            4,
            "preparation or slit-configuration labels",
            preparation_ok,
            (
                f"preparations={event_info.get('preparations')}"
                if resource_type == "heat_process"
                else f"configurations={event_info.get('configurations')}"
            ),
        ),
        Check(
            5,
            "intervention or shutter-control labels",
            intervention_ok,
            f"interventions={event_info.get('interventions')}",
        ),
        Check(6, "clock and dimensional calibration", clock_ok, str(calibration)),
        Check(7, "detector calibration", detector_ok, str(calibration)),
        Check(8, "environment and background record", environment_ok, str(environment)),
        Check(9, "negative-control record", controls_ok, str(controls)),
        Check(10, "held-out target committed before analysis", chronology_ok, chronology_detail),
        Check(
            11,
            "provider-registrar-analyst custody separation",
            role_separation_ok,
            f"roles={roles}; signature={signature_detail}",
        ),
    ]
    passed = sum(check.passed for check in checks)
    all_fields = passed == 11

    process_semigroup_ready = (
        resource_type == "heat_process"
        and all_fields
        and event_info.get("preparations") == list(range(12))
        and event_info.get("outcomes") == list(range(12))
        and set(event_info.get("evolution_multiples", [])) == {1, 2}
        and {"calibration", "holdout"}.issubset(event_info.get("subsets", []))
    )
    minimum_cell_count = (
        min(event_info.get("count_cells", {}).values())
        if event_info.get("count_cells")
        else 0
    )

    core = {
        "validator": "FIN-LAB-P241",
        "schema_version": SCHEMA_VERSION,
        "bundle_path_name": bundle.name,
        "resource_type": resource_type,
        "manifest_sha256": sha256_file(manifest_path),
        "bundle_digest": bundle_digest(bundle, actual_hashes),
        "checks": [check.as_dict() for check in checks],
        "passed_fields": passed,
        "all_eleven_fields_passed": all_fields,
        "registrar_signature_verified": signature_ok,
        "program_241_gate_passed": all_fields,
        "program_242_semigroup_ready": process_semigroup_ready,
        "minimum_process_cell_count": minimum_cell_count,
        "admission_status": (
            "ADMITTED_BY_SIGNED_REGISTRAR"
            if all_fields
            else "REJECTED_OR_AWAITING_EXTERNAL_CUSTODY"
        ),
        "software_limit": (
            "The validator checks declared metadata, hashes, event structure, "
            "chronology, role inequality, and a trusted-key signature. It "
            "does not create data or establish the real-world identity, "
            "competence, or independence of the named parties."
        ),
        "double_slit_boundary": (
            "An admitted double-slit event bundle is not automatically a "
            "twelve-state semigroup record and cannot unlock Program 242 "
            "without a separate typed process map."
        ),
    }
    return {"core": core, "canonical_core_sha256": canonical_digest(core)}


def emit_templates(target: Path) -> None:
    target.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "resource_type": "heat_process",
        "evidence_status": "EXTERNAL_REQUIRED",
        "acquisition_started_at_utc": "YYYY-MM-DDTHH:MM:SSZ",
        "provider": {
            "name": "",
            "institution": "",
            "contact": "",
            "license": "CC BY 4.0 or explicit data-use license",
        },
        "custody": {
            "provider_id": "",
            "registrar_id": "",
            "analyst_id": "",
        },
        "declaration": {
            "synthetic": False,
            "unit_test_fixture": False,
            "rendered_image_only": False,
        },
        "holdout": {
            "sealed_before_analysis": True,
            "committed_before_analysis": True,
        },
        "files": {
            name: "REPLACE_WITH_SHA256"
            for name in sorted(REQUIRED_FILES)
        },
    }
    (target / "bundle_manifest.template.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    calibration = {
        "clock_model": "",
        "clock_unit": "s",
        "clock_uncertainty": "",
        "time_origin": "",
        "detector_geometry": "",
        "efficiency_calibration": "",
        "dark_count_calibration": "",
        "blur_or_psf_calibration": "",
    }
    (target / "calibration.template.json").write_text(
        json.dumps(calibration, indent=2) + "\n", encoding="utf-8"
    )
    controls = {
        "negative_controls": [
            "time_order_permutation",
            "preparation_label_permutation",
            "nearest_neighbour_C12",
        ]
    }
    (target / "controls.template.json").write_text(
        json.dumps(controls, indent=2) + "\n", encoding="utf-8"
    )
    environment = {
        "background_record": "",
        "temperature_or_not_applicable": "",
        "source_stability_record": "",
        "apparatus_state_record": "",
    }
    (target / "environment.template.json").write_text(
        json.dumps(environment, indent=2) + "\n", encoding="utf-8"
    )
    preregistration = {
        "registered_at_utc": "YYYY-MM-DDTHH:MM:SSZ",
        "held_out_target": "P_2tau compared with P_tau@P_tau",
        "frozen_analysis_hash": "SHA256_OF_P242_ANALYSIS_LOCK",
        "one_execution_only": True,
        "failure_reported_without_model_repair": True,
    }
    (target / "preregistration.template.json").write_text(
        json.dumps(preregistration, indent=2) + "\n", encoding="utf-8"
    )
    process_header = ",".join(sorted(PROCESS_REQUIRED_COLUMNS)) + "\n"
    slit_header = ",".join(sorted(DOUBLE_SLIT_REQUIRED_COLUMNS)) + "\n"
    (target / "events_heat_process.header.csv").write_text(
        process_header, encoding="utf-8"
    )
    (target / "events_double_slit.header.csv").write_text(
        slit_header, encoding="utf-8"
    )
    readme = """FIN P241 TRANSFER TEMPLATE

These files define structure only. They contain no empirical events and
cannot pass the P241 gate.

Production sequence:
1. Provider records raw events and calibration/control/environment files.
2. Registrar replaces template names, computes every SHA-256, and completes
   bundle_manifest.json.
3. Registrar signs bundle_manifest.json with a detached GPG signature.
4. Analyst receives neither the sealed holdout nor registrar private key.
5. Run:
   python3 fin_lab_p241_validator.py BUNDLE \
     --signature bundle_manifest.json.asc \
     --trusted-keyring registrar-trustedkeys.gpg \
     --output admission_certificate.json
"""
    (target / "README.txt").write_text(readme, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", nargs="?", type=Path)
    parser.add_argument("--signature", type=Path)
    parser.add_argument("--trusted-keyring", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--emit-template", type=Path)
    args = parser.parse_args()
    if args.emit_template:
        emit_templates(args.emit_template)
        print(args.emit_template)
        return 0
    if args.bundle is None:
        parser.error("bundle is required unless --emit-template is used")
    result = validate_bundle(args.bundle, args.signature, args.trusted_keyring)
    text = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0 if result["core"]["program_241_gate_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
