#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2372_s1322_bridge_kernel_direct_band_polarity_audit.json"
MD = GEN / "p2372_s1322_bridge_kernel_direct_band_polarity_audit.md"

SOURCE_FILES = {
    "P2363_BRIDGE_MOMENTS": GEN / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.json",
    "P2371_AUT_BANDPASS_OBSTRUCTION": GEN / "p2371_s1321_aut_invariant_unit_bandpass_obstruction_theorem.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
CHAMBER_THRESHOLD = 1.0 / 3.0
STRICT_PARAMS = {"omega": 743.0 / 4000.0, "phi": 13.0 / 80.0, "beta": 1.0, "eta": 9.0 / 5.0}
LEGACY_PARAMS = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 1.0 / 100.0}


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


def k_strict(d: int) -> float:
    return math.cos(STRICT_PARAMS["omega"] * d + STRICT_PARAMS["phi"]) / (
        1.0 + STRICT_PARAMS["beta"] * d ** STRICT_PARAMS["eta"]
    )


def k_legacy_amplitude_normalized(d: int) -> float:
    return math.cos(LEGACY_PARAMS["omega"] * d + LEGACY_PARAMS["phi"]) / (1.0 + LEGACY_PARAMS["beta_tors"] * d)


def internal_edges(support: tuple[int, ...], step: int) -> int:
    support_set = set(support)
    edges = set()
    for node in support_set:
        for neighbor in ((node + step) % Z12_NODE_COUNT, (node - step) % Z12_NODE_COUNT):
            if neighbor in support_set:
                edges.add(tuple(sorted((node, neighbor))))
    return len(edges)


def support_rows() -> list[dict[str, Any]]:
    rows = []
    for support in combinations(range(Z12_NODE_COUNT), SUPPORT_SIZE):
        h1 = internal_edges(support, 1)
        h5 = internal_edges(support, 5)
        rows.append({"support": list(support), "h1": h1, "h5": h5})
    return rows


def score_maximizers(rows: list[dict[str, Any]], a: float, b: float) -> dict[str, Any]:
    scored = [(a * row["h1"] + b * row["h5"], row) for row in rows]
    maximum = max(score for score, _ in scored)
    tolerance = 1e-12
    maximizers = [row for score, row in scored if abs(score - maximum) <= tolerance]
    pair_distribution: dict[str, int] = {}
    for row in maximizers:
        key = f"{row['h1']},{row['h5']}"
        pair_distribution[key] = pair_distribution.get(key, 0) + 1
    return {
        "weights": {"a_h1": a, "b_h5": b},
        "a_over_b": None if abs(b) < tolerance else a / b,
        "d5_chamber_a_over_b_lt_1_over_3": b > 0 and a >= 0 and (a / b) < CHAMBER_THRESHOLD,
        "maximum_score": maximum,
        "maximizer_count": len(maximizers),
        "maximizer_h1_h5_pair_distribution": dict(sorted(pair_distribution.items())),
        "maximizers_sample": maximizers[:12],
    }


def polarity_rows() -> dict[str, Any]:
    rows = support_rows()
    strict_k1 = k_strict(1)
    strict_k5 = k_strict(5)
    legacy_k1 = k_legacy_amplitude_normalized(1)
    legacy_k5 = k_legacy_amplitude_normalized(5)
    candidates = {
        "strict_direct_kernel_pair_weight": (strict_k1, strict_k5),
        "strict_absolute_kernel_pair_weight": (abs(strict_k1), abs(strict_k5)),
        "apd_completed_direct_pair_weight_equals_strict": (strict_k1, strict_k5),
        "legacy_amplitude_normalized_signed_pair_weight": (legacy_k1, legacy_k5),
        "legacy_amplitude_normalized_absolute_pair_weight": (abs(legacy_k1), abs(legacy_k5)),
    }
    scored = {name: score_maximizers(rows, a, b) for name, (a, b) in candidates.items()}
    return {
        "strict_kernel_values": {"K1": strict_k1, "K5": strict_k5, "K1_over_K5": strict_k1 / strict_k5},
        "legacy_amplitude_normalized_kernel_values": {"K1": legacy_k1, "K5": legacy_k5, "K1_over_K5": legacy_k1 / legacy_k5},
        "required_chamber": "For nonnegative score a*h1+b*h5, P2371 requires a/b < 1/3 for unique d5 path selection.",
        "candidate_score_audits": scored,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    polarity = polarity_rows()

    strict_direct = polarity["candidate_score_audits"]["strict_direct_kernel_pair_weight"]
    legacy_abs = polarity["candidate_score_audits"]["legacy_amplitude_normalized_absolute_pair_weight"]

    theorem_export = {
        "theorem_name": "P2372 bridge-kernel direct band-polarity audit",
        "claim": (
            "The bridge-completed strict kernel does not, by its direct distance-1/distance-5 pair weights, "
            "export the P2371 band-polarity inequality a/b<1/3. Numerically K_strict(1)/K_strict(5) is about "
            f"{polarity['strict_kernel_values']['K1_over_K5']:.6g}, so the direct nonnegative pair score lies far outside "
            "the d5 chamber and selects distance-1 path supports rather than the P2370 distance-5 paths. Raw or amplitude-normalized "
            "legacy pair weights also fail to provide a licensed strict selector source. A separate dynamical polarity theorem remains required."
        ),
        "positive_content": [
            "P2363 bridge-completed strict kernel parameters are replayed for K_strict(d).",
            "P2371 chamber inequality a/b<1/3 is applied to direct kernel pair weights at d=1 and d=5.",
            "Strict direct and strict absolute pair weights have K1/K5 far above 1/3 and select h1 path supports.",
            "Legacy amplitude-normalized signed weights have K5<0 and are not a nonnegative d5 chamber witness.",
            "Legacy absolute pair weights still fail the a/b<1/3 chamber and are not licensed for strict role transfer.",
        ],
        "not_licensed": [
            "strict derivation of band-polarity source",
            "strict derivation of distance-5 band-pass action",
            "promotion of direct K(d) pair weights to selector action",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "QW-2191 discharge",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2372_S1322_BRIDGE_KERNEL_DIRECT_BAND_POLARITY_AUDIT",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "p2363_replay": {
            "packet_id": artifacts["P2363_BRIDGE_MOMENTS"].get("packet_id"),
            "completed_moments_match_strict": artifacts["P2363_BRIDGE_MOMENTS"].get("gatekeeper_checks", {}).get("completed_moments_match_strict"),
            "no_selector_or_role_transfer_claimed": artifacts["P2363_BRIDGE_MOMENTS"].get("gatekeeper_checks", {}).get("no_selector_or_role_transfer_claimed"),
        },
        "p2371_replay": {
            "packet_id": artifacts["P2371_AUT_BANDPASS_OBSTRUCTION"].get("packet_id"),
            "full_aut_does_not_select_d5": artifacts["P2371_AUT_BANDPASS_OBSTRUCTION"].get("gatekeeper_checks", {}).get("full_aut_does_not_select_d5"),
        },
        "polarity_audit": polarity,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2363_loaded": probe["p2363_replay"]["packet_id"] == "P2363",
        "p2371_loaded": probe["p2371_replay"]["packet_id"] == "P2371",
        "strict_direct_fails_d5_chamber": strict_direct["d5_chamber_a_over_b_lt_1_over_3"] is False,
        "strict_direct_selects_h1_paths": strict_direct["maximizer_h1_h5_pair_distribution"] == {"4,0": 12},
        "strict_ratio_above_threshold": polarity["strict_kernel_values"]["K1_over_K5"] > CHAMBER_THRESHOLD,
        "legacy_signed_not_nonnegative_d5_witness": polarity["legacy_amplitude_normalized_kernel_values"]["K5"] < 0,
        "legacy_absolute_fails_d5_chamber": legacy_abs["d5_chamber_a_over_b_lt_1_over_3"] is False,
        "docs_updated_with_p2372_audit": all("P2372/S1322" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2372_s1322_v1",
        "packet_id": "P2372",
        "stage_id": "S1322",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_BRIDGE_KERNEL_DIRECT_POLARITY_AUDIT_NO_QW2191_DISCHARGE",
        "result_kind": "BRIDGE_KERNEL_DIRECT_BAND_POLARITY_AUDIT",
        "bridge_kernel_direct_band_polarity_audit": probe,
        "recommended_next_honest_step": (
            "Do not use direct K(d) pair weights as the d5 selector source. Search for a separate "
            "bridge-completed dynamical term whose effective h1/h5 coefficient ratio satisfies a/b<1/3, "
            "or classify the band-polarity premise as non-strict."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2372 S1322: bridge-kernel direct band-polarity audit",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2372/S1322 checks whether the bridge-completed strict kernel itself supplies the P2371 band-polarity inequality for the d5 selector. It does not.",
                "",
                "## Certificate",
                "",
                f"- Strict K1: `{polarity['strict_kernel_values']['K1']}`.",
                f"- Strict K5: `{polarity['strict_kernel_values']['K5']}`.",
                f"- Strict K1/K5: `{polarity['strict_kernel_values']['K1_over_K5']}`; required `< 1/3`.",
                f"- Strict direct maximizer pair distribution: `{strict_direct['maximizer_h1_h5_pair_distribution']}`.",
                f"- Legacy normalized K5 is negative: `{polarity['legacy_amplitude_normalized_kernel_values']['K5']}`.",
                "",
                "## Hard limits",
                "",
                "- This is a direct-pair-weight audit, not a theorem that K(d) is the selector action.",
                "- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
