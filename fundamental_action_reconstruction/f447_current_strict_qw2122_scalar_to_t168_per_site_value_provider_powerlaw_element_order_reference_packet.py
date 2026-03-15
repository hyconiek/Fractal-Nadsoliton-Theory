#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

QW2122_JSON = REPO / "report_qw2122_psi_potential_diagonal_floor_gate.json"
ALPHA_GEO_JSON = GENERATED / "alpha_geo_strict_derived_v1.json"

R14_JSON = GENERATED / "r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"
R15_JSON = GENERATED / "r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"

# Persisted strict reference datum (no marked direction) + strict T169/T168 provider candidate.
OUT_REF_JSON = GENERATED / "r_ordpow_z12_v1_reference_distribution.json"
OUT_PROVIDER_JSON = GENERATED / "f447_current_strict_t169_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1.json"
OUT_PROVIDER_SUMMARY = GENERATED / "f447_current_strict_t169_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1_summary.json"

# Populate the designated harness inputs (keeps downstream probes/dashboard computable without manual copy).
P437_IN_JSON = GENERATED / "p437_input_vpsi_g4_g6_candidate.json"
P434_IN_JSON = GENERATED / "p434_input_sigma_opposite_pair_sum_values_candidate.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def parse_alpha_geo_numeric(obj: dict[str, Any]) -> float | None:
    # Accept either explicit numeric value or the canonical symbolic string "4 ln(2)".
    val = obj.get("value")
    if is_number(val):
        return float(val)
    if isinstance(val, str) and val.strip().replace(" ", "") in ("4ln(2)", "4*ln(2)", "4*log(2)"):
        return float(4.0 * math.log(2.0))
    # Fallback: use strict computation if definition says ln(16).
    if isinstance(obj.get("definition"), dict):
        d = obj["definition"]
        if isinstance(d.get("computation"), str) and "ln(16)" in d.get("computation", ""):
            return float(4.0 * math.log(2.0))
    return None


def ord_z12_by_x() -> list[int]:
    out: list[int] = []
    for x in range(12):
        if x == 0:
            out.append(1)
        else:
            out.append(int(12 // math.gcd(x, 12)))
    return out


def normalize_positive_weights(w: list[float]) -> list[float]:
    if not w or any((not math.isfinite(float(x)) or float(x) <= 0.0) for x in w):
        raise ValueError("weights must be finite positive")
    z = float(sum(float(x) for x in w))
    if not math.isfinite(z) or z <= 0.0:
        raise ValueError("weight sum must be finite positive")
    return [float(x) / z for x in w]


def extract_ktotal_from_r14(r14: dict[str, Any]) -> list[list[float]]:
    rows = r14.get("host_kernel_rows")
    if not (isinstance(rows, list) and len(rows) == 12 and all(isinstance(r, list) and len(r) == 12 for r in rows)):
        raise ValueError("R14.host_kernel_rows missing or not 12x12 list")
    k: list[list[float]] = []
    for i, row in enumerate(rows):
        rr: list[float] = []
        for j, v in enumerate(row):
            if not is_number(v):
                raise ValueError(f"R14.host_kernel_rows[{i}][{j}] not numeric")
            rr.append(float(v))
        k.append(rr)
    return k


def extract_m0_sq_from_r15(r15: dict[str, Any]) -> float:
    dd = r15.get("diagonal_decomposition")
    if not isinstance(dd, dict):
        raise ValueError("R15.diagonal_decomposition missing")
    m0 = dd.get("host_diagonal_floor_value")
    if not is_number(m0):
        raise ValueError("R15.diagonal_decomposition.host_diagonal_floor_value missing or non-numeric")
    return float(m0)


def mixing_energy_for_signed_vpsi(*, ktotal: list[list[float]], vpsi: list[float]) -> float:
    if len(vpsi) != 12:
        raise ValueError("vpsi must be length 12")
    e = 0.0
    for i in range(12):
        for j in range(i + 1, 12):
            e += float(ktotal[i][j]) * float(vpsi[i]) * float(vpsi[j])
    return float(e)


def choose_energy_minimizing_sign_vector(*, ktotal: list[list[float]], magnitudes: list[float]) -> dict[str, Any]:
    if len(magnitudes) != 12:
        raise ValueError("magnitudes must be length 12")
    mags = [float(abs(x)) for x in magnitudes]

    # Brute-force Ising minimization on 12 sites: 2^12 = 4096 states (deterministic, no hidden solver slot).
    best_e = float("inf")
    best: list[list[int]] = []
    tol = 1e-12
    for mask in range(1 << 12):
        s = [1 if ((mask >> i) & 1) == 0 else -1 for i in range(12)]
        v = [float(s[i]) * mags[i] for i in range(12)]
        e = mixing_energy_for_signed_vpsi(ktotal=ktotal, vpsi=v)
        if e < best_e - tol:
            best_e = e
            best = [s]
        elif abs(e - best_e) <= tol:
            best.append(s)

    if not best:
        raise RuntimeError("unexpected: no sign vectors scanned")

    # Fix the unavoidable global Z2 sign ambiguity deterministically by requiring s0 = +1.
    # Since the reference distribution already distinguishes the identity orbit {0}, this does not introduce a new slot.
    candidates = [s for s in best if s[0] == 1]
    if not candidates:
        candidates = [[-x for x in best[0]]]

    chosen = min(tuple(s) for s in candidates)
    chosen_list = [int(x) for x in chosen]

    return {
        "rule": "choose s ∈ {±1}^12 minimizing E_mix(s)=Σ_{i<j} K_total[i,j] (s_i|v_i|)(s_j|v_j|); fix global Z2 by s0=+1",
        "min_energy": float(best_e),
        "minimizer_count": int(len(best)),
        "chosen_sign_vector": chosen_list,
    }


def compute_d_local_residual_profile(*, ktotal: list[list[float]], m0_sq: float, vpsi: list[float], g4: list[float], g6: list[float]) -> list[float]:
    # Mirror P437/N477 computation (no additional assumptions beyond vpsi_k != 0 for division).
    if not (len(vpsi) == len(g4) == len(g6) == 12):
        raise ValueError("vpsi,g4,g6 must be length 12")
    d: list[float] = []
    for k in range(12):
        denom = float(vpsi[k])
        if denom == 0.0:
            raise ValueError("vpsi entries must be nonzero for N477 division premise")
        mix = 0.0
        for j in range(12):
            if j == k:
                continue
            mix += float(ktotal[k][j]) * (float(vpsi[j]) / denom)
        val = -mix + 2.0 * float(g4[k]) * (denom**2) + 4.0 * float(g6[k]) * (denom**4) - float(m0_sq)
        d.append(float(val))
    return d


def compute_sigma_opposite_pair_sums(d: list[float]) -> dict[str, float]:
    if len(d) != 12:
        raise ValueError("d must be length 12")
    sigmas: dict[str, float] = {}
    for k in range(6):
        sigmas[f"Sigma_psi{k}_psi{k+6}"] = float(d[k] + d[k + 6])
    return sigmas


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing_files: list[str] = []
    for p in (QW2122_JSON, ALPHA_GEO_JSON, R14_JSON, R15_JSON):
        if not p.exists():
            missing_files.append(str(p.relative_to(REPO)))
    if missing_files:
        summary = {
            "stage": "F447",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FILES",
            "missing_files": missing_files,
            "no_false_pass": True,
        }
        OUT_PROVIDER_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_PROVIDER_SUMMARY)
        return

    r2122 = load_json(QW2122_JSON)
    alpha_geo_obj = load_json(ALPHA_GEO_JSON)
    r14 = load_json(R14_JSON)
    r15 = load_json(R15_JSON)

    lambda_psi_strict = (r2122.get("inputs") or {}).get("lambda_psi_strict")
    rho_star_sq = ((r2122.get("branch_results") or {}).get("broken_branch_strict") or {}).get("rho_star_sq")

    alpha_geo = parse_alpha_geo_numeric(alpha_geo_obj)

    missing_fields: list[str] = []
    if not is_number(lambda_psi_strict):
        missing_fields.append("QW-2122.inputs.lambda_psi_strict (finite number)")
    if not is_number(rho_star_sq):
        missing_fields.append("QW-2122.branch_results.broken_branch_strict.rho_star_sq (finite number)")
    if alpha_geo is None or not math.isfinite(float(alpha_geo)) or float(alpha_geo) <= 0.0:
        missing_fields.append("alpha_geo_strict_derived_v1.value (positive numeric; e.g. 4 ln 2)")

    if missing_fields:
        summary = {
            "stage": "F447",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FIELDS",
            "missing_fields": missing_fields,
            "no_false_pass": True,
        }
        OUT_PROVIDER_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_PROVIDER_SUMMARY)
        return

    lambda_psi = float(lambda_psi_strict)
    rho_sq = float(rho_star_sq)
    alpha = float(alpha_geo)

    ords = ord_z12_by_x()

    # Strict reference distribution r_ordpow(x) ∝ ord(x)^(-alpha_geo).
    # This remains Aut(Z12)-invariant (no marked direction), but is not translation-invariant (marks identity orbit).
    w = [float(math.exp(-alpha * math.log(float(o)))) for o in ords]
    rprob = normalize_positive_weights(w)

    ref_artifact = {
        "object_id": "r_ordpow_z12_v1_reference_distribution",
        "type": "z12_reference_distribution_v1",
        "definition": "r_ordpow(x) ∝ ord_Z12(x)^(-alpha_geo)",
        "alpha_geo": {"source_object_id": "alpha_geo_strict_derived_v1", "value_numeric": alpha, "value_symbolic": "4 ln 2"},
        "carrier": {"group_object_id": "z_12_v1_group", "index_set_object_id": "i_12_v1_index_set"},
        "ord_z12_by_x": ords,
        "reference_prob": rprob,
        "invariance_notes": [
            "Aut(Z_12)-invariant reference shape (depends only on ord_Z12(x)): no marked generator/direction slot.",
            "Not translation-invariant on the regular action: distinguishes the identity orbit {0}.",
        ],
        "provenance": {"packet": "F447", "inputs": ["F329", "alpha_geo_strict_derived_v1"]},
        "no_false_pass": True,
    }

    # T169 lift: vpsi_i^2 := rho_star_sq * r_ordpow(i) (unique minimizer of KL(q||r) under sum constraint).
    vpsi = [float(math.sqrt(rho_sq * float(p))) for p in rprob]

    # Scalar-to-per-site quartic lift: enforce the QW-2122 quartic coefficient along the vacuum direction with minimal uniform g4_i.
    # If n_i := vpsi_i / rho_* then sum_i n_i^4 = sum_i q_i^2 where q_i = vpsi_i^2 / rho_*^2 = rprob_i.
    sum_q2 = float(sum(float(q) ** 2 for q in rprob))
    g4_val = float(lambda_psi / sum_q2)
    g4 = [g4_val for _ in range(12)]
    g6 = [0.0 for _ in range(12)]

    # Compute the six sigma sums via the N477/P437 formula (uses R14/R15 strict kernel+floor).
    ktotal = extract_ktotal_from_r14(r14)
    m0_sq = extract_m0_sq_from_r15(r15)

    sign_selection = choose_energy_minimizing_sign_vector(ktotal=ktotal, magnitudes=vpsi)
    svec = sign_selection["chosen_sign_vector"]
    vpsi = [float(svec[i]) * float(vpsi[i]) for i in range(12)]

    d_local = compute_d_local_residual_profile(ktotal=ktotal, m0_sq=m0_sq, vpsi=vpsi, g4=g4, g6=g6)
    sigmas = compute_sigma_opposite_pair_sums(d_local)

    provider = {
        "stage": "F447",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "STRICT_DERIVED_T169_CONSTRAINED_LIFT_VALUE_PROVIDER_PACKET",
        "classification": "strict_derived",
        "theorems": [
            "N479_CURRENT_FIRST_STRICT_Z12_ELEMENT_ORDER_REFERENCE_AUT_Z12_INVARIANCE_NO_MARKED_DIRECTION_THEOREM",
            "N483_CURRENT_FIRST_STRICT_T169_POWERLAW_ELEMENT_ORDER_CONSTRAINED_LIFT_EXISTENCE_UNIQUENESS_THEOREM",
        ],
        "object_id": "Delta_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1",
        "goal": "export_one_strict_derived_T169_constrained_lift_value_provider_for_(vpsi,g4,g6)_and_populate_designated_harness_inputs_without_hidden_direction_slots",
        "strict_inputs": {
            "QW-2122": str(QW2122_JSON.relative_to(REPO)),
            "alpha_geo_strict_derived_v1": str(ALPHA_GEO_JSON.relative_to(REPO)),
            "R14": str(R14_JSON.relative_to(REPO)),
            "R15": str(R15_JSON.relative_to(REPO)),
        },
        "reference": {"r_ordpow": str(OUT_REF_JSON.relative_to(REPO))},
        "lift_definition": {
            "vpsi_magnitudes": "vpsi_i := sqrt(rho_star_sq * r_ordpow(i))",
            "vpsi_signs": "choose s ∈ {±1}^12 minimizing E_mix(s)=Σ_{i<j} K_total[i,j] (s_i|v_i|)(s_j|v_j|); fix global Z2 by s0=+1; set vpsi_i := s_i * vpsi_i",
            "g6": "g6_i := 0 (scalar QW-2122 model has no sextic term)",
            "g4_uniform_minimal_matching": "g4_i := lambda_psi_strict / sum_i r_ordpow(i)^2 (uniform g4 matching scalar quartic along the lifted vacuum direction)",
        },
        "values": {
            "rho_star_sq": rho_sq,
            "lambda_psi_strict": lambda_psi,
            "alpha_geo_numeric": alpha,
            "sum_q2": sum_q2,
            "g4_uniform_value": g4_val,
            "m0_sq": m0_sq,
        },
        # Compatibility: P437/P444 expect numeric lists at top-level keys.
        "vpsi": vpsi,
        "g4": g4,
        "g6": g6,
        # Also persist the induced sigma-six-sums at top-level keys for direct P434 compute (T167 interface).
        **sigmas,
        "computed": {
            "d_local_residual_profile": d_local,
            "sigma_opposite_pair_sums": sigmas,
            "diagnostics": {
                "vpsi_abs_min": float(min(abs(x) for x in vpsi)),
                "vpsi_abs_ratio_max_over_min": float(max(abs(x) for x in vpsi) / min(abs(x) for x in vpsi)),
                "mixing_energy_sign_selection": sign_selection,
            },
        },
        "hard_limits": [
            "This packet exports a strict lift/value-provider for the diagonal-lane harness; it does not claim ToE closure.",
            "Does not claim global QW-2191 discharge.",
            "Does not claim strict theta export beyond the pair1 diagonal/local axis computed from the resulting F2(d).",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F447",
        "status": "PASS_STRICT_DERIVED_T169_PROVIDER_EXPORTED",
        "export_provider": str(OUT_PROVIDER_JSON.relative_to(REPO)),
        "export_reference": str(OUT_REF_JSON.relative_to(REPO)),
        "populated_harness_inputs": {
            "P437_in": str(P437_IN_JSON.relative_to(REPO)),
            "P434_in": str(P434_IN_JSON.relative_to(REPO)),
        },
        "quick_flags": {
            "vpsi_nonzero_all": bool(all(float(x) != 0.0 and math.isfinite(float(x)) for x in vpsi)),
            "sigma_keys_present": bool(all(k in sigmas for k in [f"Sigma_psi{i}_psi{i+6}" for i in range(6)])),
        },
        "no_false_pass": True,
    }

    # Populate designated harness inputs (so P432/P434/P437/P438 dashboards become decision-ready without manual copy).
    p437_in = {
        "stage": "P437_INPUT",
        "status": "STRICT_DERIVED_INPUT_FROM_F447",
        "source_provider": str(OUT_PROVIDER_JSON.relative_to(REPO)),
        "vpsi": vpsi,
        "g4": g4,
        "g6": g6,
        "no_false_pass": True,
    }

    p434_in = {
        "stage": "P434_INPUT",
        "status": "STRICT_DERIVED_INPUT_FROM_F447",
        "source_provider": str(OUT_PROVIDER_JSON.relative_to(REPO)),
        **{k: float(sigmas[k]) for k in sorted(sigmas.keys())},
        "no_false_pass": True,
    }

    OUT_REF_JSON.write_text(json.dumps(ref_artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_PROVIDER_JSON.write_text(json.dumps(provider, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_PROVIDER_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    P437_IN_JSON.write_text(json.dumps(p437_in, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    P434_IN_JSON.write_text(json.dumps(p434_in, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    print(OUT_PROVIDER_JSON)


if __name__ == "__main__":
    main()
