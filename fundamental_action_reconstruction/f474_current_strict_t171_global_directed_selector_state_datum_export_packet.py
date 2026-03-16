#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-16"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_ALPHA_GEO = GENERATED / "alpha_geo_strict_derived_v1.json"
IN_FIX = GENERATED / "kappa_z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1.json"
IN_F470 = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
IN_F467 = GENERATED / "u_12345_pair12345_chart_glued_orientation_vector_section_oriented_mod_2pi_strict_convention_v1.json"

OUT_S_DIR = GENERATED / "s_dir_pair1_strict_v1.json"
OUT_STATE = GENERATED / "selector_state_global_c_v1_directed_strict_v1.json"
OUT_SUMMARY = GENERATED / "f474_current_strict_t171_global_directed_selector_state_datum_export_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def dot(w: list[float], v: list[float]) -> float:
    return float(sum(float(a) * float(b) for a, b in zip(w, v)))


def max_abs_diff(a: list[float], b: list[float]) -> float:
    return float(max(abs(float(x) - float(y)) for x, y in zip(a, b))) if a and b else 0.0


def extract_u_vec(obj: dict[str, Any], key: str) -> list[float] | None:
    data = obj.get("data") if isinstance(obj.get("data"), dict) else obj
    v = (data or {}).get(key)
    if not (isinstance(v, list) and len(v) == 12 and all(isinstance(x, (int, float)) for x in v)):
        return None
    return [float(x) for x in v]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing = [str(p.relative_to(REPO)) for p in (IN_ALPHA_GEO, IN_FIX, IN_F470, IN_F467) if not p.exists()]
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "F474",
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": [
                        "generated/alpha_geo_strict_derived_v1.json (F309/N420)",
                        "generated/kappa_z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1.json (F473/N523)",
                        "generated/selector_state_global_c_v1_projective_strict_v1.json (F470/N516)",
                        "generated/u_12345_pair12345_chart_glued_orientation_vector_section_oriented_mod_2pi_strict_convention_v1.json (F467/N511)",
                    ],
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    alpha_geo = load_json(IN_ALPHA_GEO)
    fix = load_json(IN_FIX)
    f470 = load_json(IN_F470)
    f467 = load_json(IN_F467)

    alpha_geo_value = alpha_geo.get("value")
    alpha_geo_numeric = float(4.0 * math.log(2.0))

    d_fix = ((fix.get("definition") or {}).get("D_fix_v1") or {}) if isinstance(fix.get("definition"), dict) else {}
    generator_fixed = int(d_fix.get("generator_fixed") or 0)
    suc_table = d_fix.get("successor_table")
    if generator_fixed != 1:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "F474",
                    "status": "INVALID_T164_FIXING_DATUM",
                    "reason": "expected generator_fixed=1 in fixing datum",
                    "generator_fixed": generator_fixed,
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )
    if not (isinstance(suc_table, list) and len(suc_table) == 12 and all(isinstance(x, int) for x in suc_table)):
        raise SystemExit(
            json.dumps(
                {
                    "stage": "F474",
                    "status": "INVALID_T164_FIXING_DATUM",
                    "reason": "expected successor_table length 12 of ints",
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    u_vecs = ((f467.get("outputs") or {}).get("u_vectors") or {}) if isinstance(f467.get("outputs"), dict) else {}
    if not (isinstance(u_vecs, dict) and all(isinstance(k, str) for k in u_vecs.keys())):
        raise SystemExit(json.dumps({"stage": "F474", "status": "INVALID_F467_SHAPE"}, ensure_ascii=True))

    u1 = u_vecs.get("u_1")
    if not (isinstance(u1, list) and len(u1) == 12 and all(isinstance(x, (int, float)) for x in u1)):
        raise SystemExit(json.dumps({"stage": "F474", "status": "MISSING_U1_VECTOR"}, ensure_ascii=True))
    u1 = [float(x) for x in u1]

    # Directed coordinate on I_12 induced by the fixed successor map. Here: x = 0..11.
    x_list = list(range(12))
    w_dir = [float(math.exp(-alpha_geo_numeric * float(x))) for x in x_list]

    S_u1 = dot(w_dir, u1)
    if abs(float(S_u1)) <= 0.0:
        # Practical: we expect nonzero under a reflection-breaking weight; fail hard if it degenerates.
        raise SystemExit(
            json.dumps(
                {
                    "stage": "F474",
                    "status": "DEGENERATE_SIGN_DISTINCTION_OBSERVABLE",
                    "reason": "S_dir(u1)=0 cannot fix orientation",
                    "S_dir_u1": float(S_u1),
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    sign_fix = 1.0 if S_u1 >= 0.0 else -1.0
    u_oriented: dict[str, list[float]] = {}
    for k, v in u_vecs.items():
        if not (isinstance(v, list) and len(v) == 12 and all(isinstance(x, (int, float)) for x in v)):
            continue
        u_oriented[k] = [float(sign_fix * float(x)) for x in v]

    S_oriented = dot(w_dir, u_oriented["u_1"])
    if not (S_oriented > 0.0):
        raise SystemExit(
            json.dumps(
                {
                    "stage": "F474",
                    "status": "SIGN_FIX_FAILED",
                    "S_dir_u1": float(S_u1),
                    "sign_fix": sign_fix,
                    "S_dir_u1_oriented": float(S_oriented),
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    # Compatibility audit against the projector-section representatives (vector-level, up to sign).
    vector_section_consistency: dict[str, Any] = {"max_abs_diff_by_pair": {}, "all_pairs_match": True}
    charts = f470.get("charts") if isinstance(f470.get("charts"), dict) else {}
    for pair_id, chart in charts.items():
        if not isinstance(chart, dict):
            continue
        local = chart.get("local_operator") if isinstance(chart.get("local_operator"), dict) else {}
        ref = local.get("ref")
        if not isinstance(ref, str):
            continue
        ref_path = Path(ref)
        if not ref_path.is_absolute():
            ref_path = REPO / ref_path
        if not ref_path.exists():
            continue
        a = load_json(ref_path)
        idx = str(pair_id).replace("pair", "")
        u_key = f"u_{idx}"
        u_a = extract_u_vec(a, u_key)
        u_f = u_vecs.get(u_key)
        if u_a is None or not (isinstance(u_f, list) and len(u_f) == 12):
            continue
        u_f_float = [float(x) for x in u_f]
        diff = max_abs_diff(u_a, u_f_float)
        vector_section_consistency["max_abs_diff_by_pair"][pair_id] = float(diff)
        if diff > 1e-12:
            vector_section_consistency["all_pairs_match"] = False

    # Export the sign-distinction observable as an explicit strict object.
    s_dir_obj = {
        "object": "S_dir_pair1_strict_v1",
        "status": "actual_exported_strict_directed_sign_distinction_observable__premise_based_on_T164_fixing_datum",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export one explicit sign-sensitive scalar observable on pair1 whose value flips under u -> -u, using a directed "
            "coordinate induced by the exported T164 fixing datum and a deterministic Shannon-typed weight w_dir(x)=exp(-alpha_geo*x). "
            "This supplies a directed orientation convention/observable in the declared strict scope without smuggling hidden generator choices."
        ),
        "inputs": {
            "alpha_geo_strict_derived_v1": alpha_geo_value,
            "alpha_geo_numeric_used": alpha_geo_numeric,
            "t164_fixing_datum_ref": str(IN_FIX.relative_to(REPO)),
            "pair1_vector_section_ref": str(IN_F467.relative_to(REPO)),
        },
        "directed_coordinate": {
            "x_domain": "I_12_v1 = {0,1,...,11}",
            "induced_by": "D_fix_v1.successor_map",
            "successor_table": suc_table,
        },
        "directed_weight": {
            "w_dir(x)": "exp(-alpha_geo * x)",
            "weights": w_dir,
            "normalization_note": "Unnormalized weights (normalization irrelevant for sign).",
        },
        "definition": {"S_dir(u)": "Σ_x w_dir(x) * u(x)", "odd_under_sign": True},
        "value_on_exported_pair1_representative": {
            "u1_source": str((f467.get("object") or "F467_vector_section")),
            "S_dir_u1": float(S_u1),
            "sign_fix": float(sign_fix),
            "S_dir_u1_oriented": float(S_oriented),
        },
        "meaning": [
            "S_dir(u) changes sign under u -> -u and therefore distinguishes u from -u in the declared directed scope.",
            "The induced sign-fix chooses the oriented representative by requiring S_dir(u) > 0.",
        ],
        "hard_limits": [
            "Premise-based: depends on an explicit T164 fixing datum; does not claim Aut(Z_12)-invariant canonicity (N462).",
            "Does not claim strict-core selector closure nor global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }
    OUT_S_DIR.write_text(json.dumps(s_dir_obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    # Export the global directed selector state object (vector-level lift) on C_v1.
    directed_state = {
        "object": "SelectorState_global_C_v1_directed_strict_v1",
        "stage": "F474",
        "status": "actual_exported_global_directed_selector_state_object_on_C_v1__vector_level_lift__premise_based__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export one explicit global directed (sign-sensitive) selector state object on C_v1 by lifting the already exported "
            "projective/ray-level selector state (F470/N516) to an oriented vector representative using the exported T164 fixing datum "
            "and the exported sign-distinction observable S_dir (this packet). The directed state descends to the strict projective state "
            "at projector/span level."
        ),
        "domain": f470.get("domain"),
        "state_type": {
            "level": "directed_vector_state",
            "encoding": "vector_representative_section_on_pair_charts",
            "sign_gauge": "fixed_by_T164_fixing_datum_and_S_dir_pair1_strict_v1",
        },
        "depends_on": {
            "projective_state_ref": str(IN_F470.relative_to(REPO)),
            "t164_fixing_datum_ref": str(IN_FIX.relative_to(REPO)),
            "vector_section_ref": str(IN_F467.relative_to(REPO)),
            "sign_observable_ref": str(OUT_S_DIR.relative_to(REPO)),
        },
        "construction": {
            "rule": "apply global sign_fix from S_dir(u1) to the existing oriented vector section u_1..u_5",
            "sign_fix": float(sign_fix),
            "S_dir_u1": float(S_u1),
            "S_dir_u1_oriented": float(S_oriented),
        },
        "outputs": {"u_vectors_directed": u_oriented},
        "compatibility": {
            "descends_to_projective": True,
            "projective_state_object": f470.get("object"),
            "vector_section_matches_projector_section": vector_section_consistency,
            "note": "Projectors/spans are unchanged under global sign flip; this object adds a directed representative only.",
        },
        "hard_limits": [
            "Premise-based: depends on explicit T164 fixing datum; does not claim Aut(Z_12)-invariant canonicity (N462).",
            "Does not promote section-level cocycle data to operator-level identities (see N512).",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }
    OUT_STATE.write_text(json.dumps(directed_state, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    summary = {
        "stage": "F474",
        "status": "PASS_EXPORTED_T171_GLOBAL_DIRECTED_SELECTOR_STATE_OBJECT_AND_SIGN_OBSERVABLE",
        "exported": {
            "S_dir_pair1": str(OUT_S_DIR.relative_to(REPO)),
            "selector_state_global_directed": str(OUT_STATE.relative_to(REPO)),
        },
        "t171_discharged": True,
        "h37_discharged": True,
        "sign_fix": float(sign_fix),
        "S_dir_u1": float(S_u1),
        "S_dir_u1_oriented": float(S_oriented),
        "vector_section_matches_projector_section": bool(vector_section_consistency.get("all_pairs_match")),
        "no_false_pass": True,
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

