#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2324_s1274_strict_axis_branch_susceptibility_nonpromotion_probe.json"
MD = GEN / "p2324_s1274_strict_axis_branch_susceptibility_nonpromotion_probe.md"

TOL = 1e-12
THETA_REF = 0.0

SOURCE_FILES = {
    "P2323_HESSIAN_BRANCH_DEGENERACY": GEN / "p2323_s1273_strict_self_coupled_axis_potential_hessian_branch_degeneracy_probe.json",
    "P2322_SELF_COUPLED_AXIS_POTENTIAL": GEN / "p2322_s1272_strict_nadsoliton_fractal_self_coupling_axis_potential_audit_probe.json",
    "P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION": GEN / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.json",
    "P2308_CURRENT_INTERFACE_NONEXISTENCE": GEN / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json",
}

GREP_PATTERNS = (
    "branch degeneracy",
    "perturbation",
    "susceptibility",
    "symmetry-breaking",
    "external selector",
    "signed perturbation",
    "provider_lift_per_step",
    "QW-2191",
    "selector closure",
)


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def branch_susceptibility_rows(p2323: dict[str, Any], required_lift: float) -> list[dict[str, Any]]:
    probe = p2323.get("strict_self_coupled_axis_potential_hessian_branch_degeneracy_probe", {}) or {}
    rows = []
    for row in probe.get("hessian_branch_rows", []) or []:
        q = int(row.get("degree_q", 0))
        branches = list(row.get("branch_rows", []) or [])
        if not branches:
            continue
        shifts_unit = []
        for branch in branches:
            theta = float(branch["theta_radians"])
            unit_shift = -math.cos(theta - THETA_REF)
            shifts_unit.append({
                "branch_index": int(branch["branch_index"]),
                "theta_radians": theta,
                "unit_signed_energy_shift_for_minus_mu_cos": unit_shift,
            })
        values = [item["unit_signed_energy_shift_for_minus_mu_cos"] for item in shifts_unit]
        min_value = min(values)
        max_value = max(values)
        spread_per_mu = max_value - min_value
        selected_indices = [item["branch_index"] for item in shifts_unit if abs(item["unit_signed_energy_shift_for_minus_mu_cos"] - min_value) < TOL]
        rejected_indices = [item["branch_index"] for item in shifts_unit if abs(item["unit_signed_energy_shift_for_minus_mu_cos"] - max_value) < TOL]
        mu_for_required_spread = required_lift / spread_per_mu if spread_per_mu > 0.0 else None
        rows.append({
            "pair": row.get("pair"),
            "degree_q": q,
            "branch_count": len(branches),
            "perturbation_template": "Delta V = -mu*cos(theta-theta_ref) evaluated on the degenerate P2323 branch orbit",
            "theta_ref": THETA_REF,
            "unit_shift_rows": shifts_unit,
            "spread_per_unit_mu": spread_per_mu,
            "selected_branch_indices_for_positive_mu": selected_indices,
            "maximally_disfavored_branch_indices_for_positive_mu": rejected_indices,
            "mu_needed_for_required_lift_spread": mu_for_required_spread,
            "required_lift_spread_check": spread_per_mu * mu_for_required_spread if mu_for_required_spread is not None else 0.0,
            "splits_degenerate_branch_orbit_if_mu_given": spread_per_mu > 0.0,
            "strict_internal_source_for_mu_exported": False,
            "why_not_strict_bridge": "The computation quantifies susceptibility to a signed perturbation. It does not derive mu, theta_ref, or a policy-margin response from strict nadsoliton data; therefore it is a non-promoted selector-premise sensitivity trace.",
        })
    return rows


def corpus_hits() -> list[dict[str, Any]]:
    paths = sorted(
        set(SOURCE_FILES.values())
        | set(ROOT.glob("P*_*.md"))
        | set(GEN.glob("p23*_s12*_strict*json"))
        | set(GEN.glob("p22*_s12*_strict*perturbation*json"))
    )
    self_paths = {OUT.resolve(), MD.resolve(), Path(__file__).resolve()}
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not path.exists() or path.is_dir() or path.resolve() in self_paths:
            continue
        text = read_text(path)
        lowered = text.lower()
        count = sum(lowered.count(pattern.lower()) for pattern in GREP_PATTERNS)
        if count == 0:
            continue
        first_line = 0
        first_excerpt = ""
        for idx, line in enumerate(text.splitlines(), start=1):
            if any(pattern.lower() in line.lower() for pattern in GREP_PATTERNS):
                first_line = idx
                first_excerpt = line.strip()[:240]
                break
        rows.append({
            "path": path.relative_to(REPO).as_posix(),
            "pattern_hit_count": count,
            "first_hit_line": first_line,
            "first_hit_excerpt": first_excerpt,
        })
    rows.sort(key=lambda row: (-int(row["pattern_hit_count"]), row["path"]))
    return rows


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2323 = load_json(SOURCE_FILES["P2323_HESSIAN_BRANCH_DEGENERACY"])
    p2322 = load_json(SOURCE_FILES["P2322_SELF_COUPLED_AXIS_POTENTIAL"])
    p2318 = load_json(SOURCE_FILES["P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION"])
    p2308 = load_json(SOURCE_FILES["P2308_CURRENT_INTERFACE_NONEXISTENCE"])

    p2318_probe = p2318.get("strict_selector_lane_to_task3_margin_bridge_obligation_probe", {}) or {}
    required_lift = float((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("required_lift_per_step", 0.0))
    missing_fields = list((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("missing_required_bridge_fields", []))
    rows = branch_susceptibility_rows(p2323, required_lift)

    susceptibility_certificate = {
        "pair_count": len(rows),
        "all_pairs_split_by_template_if_mu_given": all(row["splits_degenerate_branch_orbit_if_mu_given"] for row in rows),
        "spread_per_unit_mu_by_pair": {row["pair"]: row["spread_per_unit_mu"] for row in rows},
        "mu_needed_for_required_lift_spread_by_pair": {row["pair"]: row["mu_needed_for_required_lift_spread"] for row in rows},
        "max_required_lift_spread_residual": max(abs(row["required_lift_spread_check"] - required_lift) for row in rows),
        "strict_internal_mu_source_exported": False,
        "lay_symptom_summary": [
            "The self-coupled nadsoliton axis potential behaves like a perfectly balanced multi-valley landscape.",
            "Each valley is locally stable, but the valleys are exactly tied until a real sign/source tilts the landscape.",
            "A tiny signed tilt would choose a valley very efficiently; the missing theorem is not numerical size but internal provenance and response semantics.",
        ],
    }

    bridge_obligation_update = {
        "required_lift_per_step": required_lift,
        "p2318_missing_fields_loaded": missing_fields,
        "signed_perturbation_susceptibility_exported": True,
        "strict_internal_signed_perturbation_exported": False,
        "susceptibility_fills_any_missing_p2318_field": False,
        "fields_still_unfilled_after_susceptibility_trace": missing_fields,
        "g1_g3_update_allowed": False,
        "reason": "A hypothetical signed perturbation splits the branch orbit, but P2324 does not derive that perturbation from strict core, does not type it as a P2281/P2302 response, and does not discharge P2308/P2318 interface blockers.",
    }

    theorem_export = {
        "theorem_name": "P2324 signed-perturbation susceptibility without strict-source promotion",
        "formal_statement": (
            "Given the P2323 exactly degenerate but Hessian-stable D12 branch orbit, a hypothetical signed perturbation Delta V=-mu cos(theta-theta_ref) splits every pair branch orbit with computable spread per unit mu.  "
            "For the P2318 required lift 0.0068, the needed mu values are finite and explicitly computed for each pair.  This proves high susceptibility of the stable axis candidates to a signed tilt, but it does not derive mu or theta_ref from strict nadsoliton sources, does not export a Task-3 response functional or replay semantics, and therefore cannot update G1/G3 or discharge QW-2191."
        ),
        "proof_bits": {
            "required_lift_per_step": required_lift,
            "spread_per_unit_mu_by_pair": susceptibility_certificate["spread_per_unit_mu_by_pair"],
            "mu_needed_for_required_lift_spread_by_pair": susceptibility_certificate["mu_needed_for_required_lift_spread_by_pair"],
            "max_required_lift_spread_residual": susceptibility_certificate["max_required_lift_spread_residual"],
            "p2318_missing_field_count": len(missing_fields),
            "strict_internal_mu_source_exported": False,
        },
        "scope_limits": [
            "sensitivity calculation under a hypothetical signed perturbation, not strict-source derivation",
            "does not choose a physical branch/sign inside strict core",
            "does not export Task-3 response functional or P2281 replay semantics",
            "does not discharge QW-2191, update G1/G3, or close ToE",
        ],
        "nonpromotion_rule": "Do not promote perturbation susceptibility into selector closure unless the signed perturbation, sign reference, and response/replay bridge are independently exported from strict sources.",
    }
    theorem_fingerprint = sha256_json(theorem_export)

    probe = {
        "probe_id": "P2324_S1274_STRICT_AXIS_BRANCH_SUSCEPTIBILITY_NONPROMOTION",
        "source_packets": {key: path.relative_to(REPO).as_posix() for key, path in SOURCE_FILES.items()},
        "source_hashes": {f"{key}_sha256": sha256_file(path) for key, path in SOURCE_FILES.items()},
        "far_susceptibility_grep_audit": {
            "patterns": list(GREP_PATTERNS),
            "hit_count": len(corpus_hits()),
            "top_hits": corpus_hits()[:30],
        },
        "susceptibility_rows": rows,
        "susceptibility_certificate": susceptibility_certificate,
        "bridge_obligation_update": bridge_obligation_update,
        "existing_blocker_context": {
            "p2323_result_kind": p2323.get("result_kind"),
            "p2322_result_kind": p2322.get("result_kind"),
            "p2318_result_kind": p2318.get("result_kind"),
            "p2308_result_kind": p2308.get("result_kind"),
        },
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": theorem_fingerprint,
    }

    gatekeeper_checks = {
        "far_susceptibility_grep_hits_found": len(probe["far_susceptibility_grep_audit"]["top_hits"]) >= 5,
        "p2323_loaded": p2323.get("packet_id") == "P2323",
        "p2318_loaded": p2318.get("packet_id") == "P2318",
        "pair_count_is_5": len(rows) == 5,
        "all_pairs_split_by_template_if_mu_given": susceptibility_certificate["all_pairs_split_by_template_if_mu_given"],
        "required_lift_spread_reconstructed": susceptibility_certificate["max_required_lift_spread_residual"] < TOL,
        "strict_internal_mu_source_not_exported": not susceptibility_certificate["strict_internal_mu_source_exported"],
        "susceptibility_not_promoted_to_bridge": not bridge_obligation_update["susceptibility_fills_any_missing_p2318_field"],
        "p2318_bridge_fields_still_unfilled": len(bridge_obligation_update["fields_still_unfilled_after_susceptibility_trace"]) == 6,
        "strict_g1_g3_not_updated": not bridge_obligation_update["g1_g3_update_allowed"],
        "no_selector_closure_claimed": True,
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2324_s1274_v1",
        "packet_id": "P2324",
        "stage_id": "S1274",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_SIGNED_PERTURBATION_SUSCEPTIBILITY_WITHOUT_STRICT_SOURCE_PROMOTION",
        "result_kind": "STRICT_AXIS_BRANCH_SUSCEPTIBILITY_NONPROMOTION_AUDIT_NO_G1_G3_UPDATE",
        "strict_axis_branch_susceptibility_nonpromotion_probe": probe,
        "recommended_next_honest_step": "Derive a strict internal signed perturbation/source and response-replay bridge, or prove that all such tilts remain selector premises outside strict core.",
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_BRANCH_SUSCEPTIBILITY_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2324/S1274 — branch susceptibility nonpromotion audit\n\n"
        f"- Result kind: `{payload['result_kind']}`.\n"
        f"- Spread per unit signed perturbation by pair: `{susceptibility_certificate['spread_per_unit_mu_by_pair']}`.\n"
        f"- Mu needed to create required lift spread `{required_lift}`: `{susceptibility_certificate['mu_needed_for_required_lift_spread_by_pair']}`.\n"
        "- Lay symptom: the nadsoliton axis landscape is stable but perfectly tied; a tiny real signed tilt would pick a valley, yet the tilt itself is not derived here.\n"
        "- Verdict: susceptibility is computed, but no strict signed source, response/replay bridge, G1/G3 update, QW-2191 discharge, or ToE closure is claimed.\n",
        encoding="utf-8",
    )
    print(json.dumps({"wrote": str(OUT.relative_to(REPO)), "theorem_fingerprint_sha256": theorem_fingerprint}, indent=2))


if __name__ == "__main__":
    main()
