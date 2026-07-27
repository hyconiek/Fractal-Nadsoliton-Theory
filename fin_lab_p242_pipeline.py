#!/usr/bin/env python3
"""FIN Laboratory Program 242: one-shot external semigroup pipeline.

The pipeline is deliberately fail-closed.  It re-runs Program 241, verifies
the frozen analysis-lock hash, requires a twelve-preparation heat-process
bundle, atomically consumes one execution token for the admitted bundle, and
then evaluates the registered semigroup holdout exactly once.

A vanilla double-slit record can pass the P241 event/custody contract but is
not a 12x12 transition record.  It therefore cannot unlock this semigroup
pipeline without a separately preregistered typed process map.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import sys

import numpy as np
from scipy.linalg import expm

import fin_lab_p240_optimal_tomography as p240
import fin_lab_p241_validator as p241


ROOT = Path(__file__).resolve().parent
ANALYSIS_LOCK = ROOT / "FIN_Lab_P242_Analysis_Lock.json"
ALPHA_TOTAL = 0.05
FINGERPRINT_THRESHOLD = 0.02
N = 12


def canonical_digest(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def analysis_lock_core() -> dict:
    return {
        "protocol_id": "FIN-LAB-P242-EXTERNAL-SEMIGROUP-HOLDOUT-V1",
        "depends_on": [
            "FIN-LAB-P240-OPTIMAL-HEAT-TOMOGRAPHY-V1",
            "FIN-P241-1.0 signed admission",
        ],
        "resource_type": "heat_process",
        "required_dimension": 12,
        "calibration_record": "subset=calibration and evolution_multiple=1",
        "held_out_record": "subset=holdout and evolution_multiple=2",
        "prediction": "P_2tau_predicted=P_tau@P_tau",
        "primary_statistic": (
            "maximum over preparation columns of total-variation distance "
            "between empirical P_2tau and P_tau@P_tau"
        ),
        "primary_alpha_total": ALPHA_TOTAL,
        "primary_finite_count_bound": (
            "Bretagnolle-Huber-Carol L1 multinomial bounds, unioned over "
            "twelve columns and propagated through the stochastic product"
        ),
        "primary_decision_labels": [
            "FALSIFIED_AT_REGISTERED_LEVEL",
            "NOT_FALSIFIED_AT_REGISTERED_LEVEL",
        ],
        "secondary_statistic": (
            "maximum absolute distance of the reconstructed projective "
            "positive spectrum from the frozen strict C12 fingerprint"
        ),
        "secondary_raw_threshold": FINGERPRINT_THRESHOLD,
        "negative_controls": [
            "cyclic preparation-label shift",
            "time reversal/misassignment",
            "best calibration-fit nearest-neighbour C12 heat model",
        ],
        "execution_rule": (
            "one atomic execution claim per admitted bundle digest; report "
            "failure and every control without model or threshold repair"
        ),
        "nonpromotion_rule": (
            "NOT_FALSIFIED is not proof that nature selected FIN; an emulator "
            "implementation tests the registered operator/instrument map only"
        ),
        "double_slit_rule": (
            "event-level double-slit data remain a P241 operational/custody "
            "record unless a separate 12-state process map was preregistered"
        ),
    }


def write_analysis_lock(path: Path = ANALYSIS_LOCK) -> dict:
    core = analysis_lock_core()
    record = {"core": core, "canonical_core_sha256": canonical_digest(core)}
    path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return record


def verify_analysis_lock(path: Path = ANALYSIS_LOCK) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        record = json.load(handle)
    expected = canonical_digest(record["core"])
    if record.get("canonical_core_sha256") != expected:
        raise ValueError("analysis lock canonical digest mismatch")
    if record["core"] != analysis_lock_core():
        raise ValueError("analysis lock content differs from executable specification")
    return record


def load_process_counts(events_path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    counts_tau = np.zeros((N, N), dtype=np.int64)
    counts_2tau = np.zeros((N, N), dtype=np.int64)
    ignored = 0
    with events_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            preparation = int(row["preparation_id"])
            outcome = int(row["outcome_id"])
            multiple = int(row["evolution_multiple"])
            subset = row["subset"].strip()
            if subset == "calibration" and multiple == 1:
                counts_tau[outcome, preparation] += 1
            elif subset == "holdout" and multiple == 2:
                counts_2tau[outcome, preparation] += 1
            else:
                ignored += 1
    totals_tau = counts_tau.sum(axis=0)
    totals_2tau = counts_2tau.sum(axis=0)
    if np.any(totals_tau <= 0) or np.any(totals_2tau <= 0):
        raise ValueError(
            "all twelve calibration P_tau and held-out P_2tau columns are required"
        )
    return counts_tau, counts_2tau, {
        "shots_tau_by_preparation": totals_tau.tolist(),
        "shots_2tau_by_preparation": totals_2tau.tolist(),
        "ignored_nonregistered_events": ignored,
    }


def transition_from_counts(counts: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=0)
    if np.any(totals <= 0):
        raise ValueError("empty preparation column")
    return counts / totals[None, :]


def max_column_tv(left: np.ndarray, right: np.ndarray) -> float:
    if left.shape != right.shape:
        raise ValueError("transition shapes differ")
    return float(0.5 * np.max(np.sum(np.abs(left - right), axis=0)))


def bhc_l1_radius(categories: int, shots: int, alpha: float) -> float:
    """Bretagnolle-Huber-Carol radius, union-ready for one column."""
    if categories <= 1 or shots <= 0 or not 0.0 < alpha < 1.0:
        raise ValueError("invalid BHC arguments")
    return min(
        2.0,
        math.sqrt(
            2.0
            * (categories * math.log(2.0) + math.log(1.0 / alpha))
            / shots
        ),
    )


def registered_primary_radius(
    shots_tau_min: int,
    shots_2tau_min: int,
    alpha_total: float = ALPHA_TOTAL,
) -> dict:
    # Half of alpha to each empirical time; union each half over N columns.
    alpha_per_tau_column = (alpha_total / 2.0) / N
    alpha_per_2tau_column = (alpha_total / 2.0) / N
    delta_tau_l1 = bhc_l1_radius(N, shots_tau_min, alpha_per_tau_column)
    delta_2tau_l1 = bhc_l1_radius(N, shots_2tau_min, alpha_per_2tau_column)
    # ||Phat^2-P^2||_1 <= 2||Phat-P||_1 for column-stochastic matrices.
    # TV contributes half the L1 error, hence delta_tau_l1.  Direct held-out
    # empirical TV contributes delta_2tau_l1/2.
    total_tv_radius = min(1.0, delta_tau_l1 + 0.5 * delta_2tau_l1)
    return {
        "alpha_total": alpha_total,
        "alpha_per_tau_column": alpha_per_tau_column,
        "alpha_per_2tau_column": alpha_per_2tau_column,
        "delta_tau_l1": delta_tau_l1,
        "delta_2tau_l1": delta_2tau_l1,
        "registered_tv_radius": total_tv_radius,
    }


def projective_fingerprint(transition: np.ndarray) -> tuple[np.ndarray | None, dict]:
    symmetric = (transition + transition.T) / 2.0
    eigenvalues = np.linalg.eigvalsh(symmetric)[::-1]
    diagnostic = {
        "transition_eigenvalues_descending": eigenvalues.tolist(),
        "minimum_transition_eigenvalue": float(eigenvalues[-1]),
    }
    if eigenvalues[-1] <= 0 or eigenvalues[0] <= 0:
        return None, diagnostic
    rates = -np.log(np.clip(eigenvalues / eigenvalues[0], 1e-15, 1.0))
    positive = np.sort(rates[1:])
    if positive[-1] <= 0:
        return None, diagnostic
    return positive / positive[-1], diagnostic


def nearest_neighbour_control(p_tau: np.ndarray, p_2tau: np.ndarray) -> dict:
    laplacian = np.zeros((N, N), dtype=float)
    for i in range(N):
        laplacian[i, i] = 2.0
        laplacian[(i - 1) % N, i] = -1.0
        laplacian[(i + 1) % N, i] = -1.0
    scale_grid = np.geomspace(0.003, 8.0, 1400)
    best_scale = None
    best_calibration = math.inf
    best_transition = None
    for scale in scale_grid:
        candidate = expm(-float(scale) * laplacian)
        distance = max_column_tv(candidate, p_tau)
        if distance < best_calibration:
            best_calibration = distance
            best_scale = float(scale)
            best_transition = candidate
    assert best_transition is not None and best_scale is not None
    prediction = expm(-2.0 * best_scale * laplacian)
    return {
        "fitted_scale": best_scale,
        "calibration_max_column_tv": best_calibration,
        "heldout_max_column_tv": max_column_tv(prediction, p_2tau),
    }


def analyze_counts(counts_tau: np.ndarray, counts_2tau: np.ndarray) -> dict:
    p_tau = transition_from_counts(counts_tau)
    p_2tau = transition_from_counts(counts_2tau)
    prediction = p_tau @ p_tau
    observed_primary = max_column_tv(prediction, p_2tau)
    shots_tau_min = int(np.min(counts_tau.sum(axis=0)))
    shots_2tau_min = int(np.min(counts_2tau.sum(axis=0)))
    primary_bound = registered_primary_radius(shots_tau_min, shots_2tau_min)
    primary_decision = (
        "FALSIFIED_AT_REGISTERED_LEVEL"
        if observed_primary > primary_bound["registered_tv_radius"]
        else "NOT_FALSIFIED_AT_REGISTERED_LEVEL"
    )

    a, _, _ = p240.strict_operator()
    target_signature = p240.projective_signature_from_generator(a)
    observed_signature, spectral_diagnostic = projective_fingerprint(p_tau)
    if observed_signature is None:
        fingerprint_distance = None
        fingerprint_decision = "FAIL_NONPOSITIVE_TRANSITION_SPECTRUM"
    else:
        fingerprint_distance = float(
            np.max(np.abs(observed_signature - target_signature))
        )
        fingerprint_decision = (
            "PASS_RAW_0P02_GATE"
            if fingerprint_distance <= FINGERPRINT_THRESHOLD
            else "FAIL_RAW_0P02_GATE"
        )

    shifted_calibration = np.roll(p_tau, 1, axis=1)
    shifted_prediction = p_tau @ shifted_calibration
    reverse_misassignment = p_2tau @ p_2tau
    controls = {
        "cyclic_preparation_label_shift_max_column_tv": max_column_tv(
            shifted_prediction, p_2tau
        ),
        "time_reversal_misassignment_max_column_tv": max_column_tv(
            reverse_misassignment, p_tau
        ),
        "nearest_neighbour_C12": nearest_neighbour_control(p_tau, p_2tau),
    }
    return {
        "primary": {
            "statistic": "maximum column total-variation distance",
            "observed": observed_primary,
            "finite_count_bound": primary_bound,
            "decision": primary_decision,
        },
        "secondary_fingerprint": {
            "observed_projective_signature": (
                observed_signature.tolist() if observed_signature is not None else None
            ),
            "strict_target_signature": target_signature.tolist(),
            "maximum_absolute_distance": fingerprint_distance,
            "raw_threshold": FINGERPRINT_THRESHOLD,
            "decision": fingerprint_decision,
            "diagnostic": spectral_diagnostic,
        },
        "negative_controls": controls,
        "nonpromotion": (
            "A non-rejection or raw fingerprint pass is not proof that nature "
            "selected FIN. Report model alternatives, calibration dependence, "
            "and whether the platform was engineered to realize the target."
        ),
    }


def consume_execution_token(ledger_dir: Path, bundle_digest: str) -> Path:
    ledger_dir.mkdir(parents=True, exist_ok=True)
    token = ledger_dir / f"FIN_P242_EXECUTED_{bundle_digest}.lock"
    descriptor = os.open(
        token,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL,
        0o444,
    )
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        handle.write(
            json.dumps(
                {
                    "bundle_digest": bundle_digest,
                    "protocol_id": analysis_lock_core()["protocol_id"],
                    "state": "EXECUTION_TOKEN_CONSUMED_BEFORE_UNBLINDED_ANALYSIS",
                },
                indent=2,
            )
            + "\n"
        )
    return token


def execute_once(
    bundle: Path,
    signature: Path,
    trusted_keyring: Path,
    ledger_dir: Path,
    output: Path,
) -> dict:
    lock = verify_analysis_lock()
    admission = p241.validate_bundle(bundle, signature, trusted_keyring)
    admission_core = admission["core"]
    if not admission_core["program_241_gate_passed"]:
        raise PermissionError("P241 signed 11/11 gate did not pass")
    if not admission_core["program_242_semigroup_ready"]:
        raise PermissionError(
            "bundle may be an admitted operational record but is not a "
            "twelve-state heat-process record ready for P242"
        )
    preregistration = p241.load_json(bundle / "preregistration.json")
    frozen_hash = preregistration.get("frozen_analysis_hash")
    actual_lock_hash = sha256_file(ANALYSIS_LOCK)
    if frozen_hash != actual_lock_hash:
        raise PermissionError(
            "preregistration does not commit the exact P242 analysis-lock hash"
        )

    token = consume_execution_token(
        ledger_dir, admission_core["bundle_digest"]
    )
    counts_tau, counts_2tau, count_diagnostic = load_process_counts(
        bundle / "events.csv"
    )
    analysis = analyze_counts(counts_tau, counts_2tau)
    core = {
        "protocol_id": analysis_lock_core()["protocol_id"],
        "bundle_digest": admission_core["bundle_digest"],
        "manifest_sha256": admission_core["manifest_sha256"],
        "analysis_lock_sha256": actual_lock_hash,
        "execution_token": str(token),
        "p241_admission": admission_core,
        "count_diagnostic": count_diagnostic,
        "analysis": analysis,
        "one_execution_completed": True,
        "posthoc_model_repair_permitted": False,
        "external_physical_validation_label": (
            "EXTERNAL_REGISTERED_TEST_EXECUTED; INTERPRET DECISIONS EXACTLY "
            "AS REPORTED, WITHOUT CALLING NONREJECTION PROOF"
        ),
    }
    record = {"core": core, "canonical_core_sha256": canonical_digest(core)}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-lock-only", action="store_true")
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--signature", type=Path)
    parser.add_argument("--trusted-keyring", type=Path)
    parser.add_argument("--ledger-dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.write_lock_only:
        write_analysis_lock()
        print(
            json.dumps(
                {
                    "analysis_lock": str(ANALYSIS_LOCK),
                    "sha256": sha256_file(ANALYSIS_LOCK),
                    "external_execution": False,
                },
                indent=2,
            )
        )
        return 0
    required = {
        "--bundle": args.bundle,
        "--signature": args.signature,
        "--trusted-keyring": args.trusted_keyring,
        "--ledger-dir": args.ledger_dir,
        "--output": args.output,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        parser.error(f"missing production arguments: {', '.join(missing)}")
    try:
        record = execute_once(
            args.bundle,
            args.signature,
            args.trusted_keyring,
            args.ledger_dir,
            args.output,
        )
    except FileExistsError:
        print(
            "P242 REFUSED: this bundle digest already consumed its one-shot "
            "execution token.",
            file=sys.stderr,
        )
        return 4
    except (PermissionError, ValueError, FileNotFoundError) as error:
        print(f"P242 LOCKED: {error}", file=sys.stderr)
        return 3
    print(json.dumps(record, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
