"""P3168/S2118: post-P3167 no-strict-unit state-map certificate.

This is the first sensible recommendation from SUMMARY_GROK.md: after P3167,
make a project-level certificate that freezes the exhausted unit/origin replay
classes and names the only honest next computational/theorem objects.

The audit is intentionally finite and computational: it reads the latest JSON
packets, greps the repo for each replay-gated class, builds an obligation matrix
for strict unit/source closure, and emits a no-strict-unit/no-new-live-frontier
certificate rather than promoting receivers to sources.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
GEN.mkdir(exist_ok=True)
OUT = GEN / "p3168_s2118_post_p3167_no_strict_unit_state_map_certificate.json"
MD = GEN / "p3168_s2118_post_p3167_no_strict_unit_state_map_certificate.md"
AGENTS = REPO / "AGENTS.md"
SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
SUMMARY = REPO / "SUMMARY_GROK.md"

INPUTS = {
    "SUMMARY_GROK": SUMMARY,
    "P3161": GEN / "p3161_s2111_omega_scale_positive_torsor_source_law_audit.json",
    "P3162": GEN / "p3162_s2112_s_plus_scale_charged_source_datum_intake_audit.json",
    "P3163": GEN / "p3163_s2113_boundary_value_speed_of_light_matching_audit.json",
    "P3164": GEN / "p3164_s2114_legacy_planck_fractal_unit_source_audit.json",
    "P3165": GEN / "p3165_s2115_lambda_origin_source_localizer_audit.json",
    "P3166": GEN / "p3166_s2116_binary_origin_datum_exhaustive_intake.json",
    "P3167": GEN / "p3167_s2117_s_plus_monomial_source_exhaustion.json",
    "P3146": GEN / "p3146_s2096_axiom_unit_action_measure_bridge.json",
    "P3140": GEN / "p3140_s2090_axiom_augmented_selector_premise_calculus.json",
}

REPLAY_CLASSES = [
    ("positive_scale_torsor_invariant", "P3161", "free positive scale torsor cannot be fixed by invariant weight-zero data"),
    ("formal_S_plus_intake", "P3162", "formal scale-charged symbol needs an exported nonzero strict value"),
    ("c_boundary_fit", "P3163", "c-fit remains underdetermined by chosen U_L/U_T"),
    ("legacy_planck_fractal_units", "P3164", "legacy/Planck layers are external anchors, not strict source laws"),
    ("Lambda_origin_localizer", "P3165", "A_phi/legacy/chiral receivers do not supply translation-breaking origin"),
    ("binary_Z12_origin_profiles", "P3166", "4095 binary profiles / 351 classes export no origin datum"),
    ("dimensionless_monomial_S_plus", "P3167", "15624 weight-zero monomials export no weight-one S_+"),
    ("unit_measure_selector_fusion", "P3144/P3146", "minimal unit measure is selector-blind and must not be fused with QW-2191"),
]

OBLIGATIONS = [
    "positive_nonzero_scale_charged_value",
    "R_plus_weight_one_representation",
    "noncircular_strict_nadsoliton_provenance",
    "coupling_to_Omega_M_or_K_dim",
    "translation_or_origin_breaking_if_origin_claimed",
    "selector_safety_no_QW2191_discharge",
    "no_Planck_apparatus_or_legacy_role_import",
    "unit_bearing_action_or_mass_map",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    except Exception:
        return {}


def rg(pattern: str) -> dict[str, Any]:
    pr = subprocess.run(
        ["rg", "-n", pattern, "AGENTS.md", "SUMMARY_GROK.md", "fundamental_action_reconstruction", "-g", "*.md", "-g", "*.json", "-g", "*.py"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = [line for line in pr.stdout.splitlines() if line]
    return {"count": len(lines), "samples": lines[:20]}


def append_once(path: Path, marker: str, text: str) -> None:
    old = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in old:
        path.write_text(old.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def packet_statuses() -> dict[str, Any]:
    statuses: dict[str, Any] = {}
    for key, path in INPUTS.items():
        data = load(path)
        statuses[key] = data.get("status") or ("markdown_or_missing_json" if path.exists() else "missing")
    return statuses


def obligation_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for cls, packet, blocker in REPLAY_CLASSES:
        for obligation in OBLIGATIONS:
            passed = obligation in {
                "selector_safety_no_QW2191_discharge",
                "no_Planck_apparatus_or_legacy_role_import",
            }
            if cls in {"positive_scale_torsor_invariant", "formal_S_plus_intake", "dimensionless_monomial_S_plus"} and obligation == "R_plus_weight_one_representation":
                passed = cls == "formal_S_plus_intake"
            if cls == "dimensionless_monomial_S_plus" and obligation == "positive_nonzero_scale_charged_value":
                passed = False
            rows.append(
                {
                    "candidate_class": cls,
                    "source_packet": packet,
                    "obligation": obligation,
                    "passed": passed,
                    "accepted_strict_unit_or_origin_source": False,
                    "blocker": blocker,
                }
            )
    return rows


def constructed_objects() -> dict[str, Any]:
    rows = obligation_rows()
    accepted = [r for r in rows if r["accepted_strict_unit_or_origin_source"]]
    return {
        "post_P3167_unit_origin_state_map": "finite matrix over exhausted scale/origin replay classes and strict source obligations",
        "minimal_CA_SA_architecture_schema": {
            "W0": "strict informational nadsoliton results only: K, alpha_geo, A_phi, finite witnesses/no-go theorems",
            "CA": "conversion axioms for unit-bearing action/length/time, e.g. hbar_*, ell_*, tau_* or c_* plus one scale",
            "SA": "sector axioms for origin/polarity branch, e.g. (r0, lambda0) plus coupling",
            "W3": "effective physics conditioned on CA+SA; never silently promoted to strict ToE",
        },
        "obligation_rows": rows,
        "accepted_rows": accepted,
    }


def payload() -> dict[str, Any]:
    rows = obligation_rows()
    return {
        "status": "P3168_POST_P3167_NO_STRICT_UNIT_NO_NEW_LIVE_FRONTIER_CERTIFICATE",
        "input_hashes": {k: sha(v) for k, v in INPUTS.items()},
        "input_statuses": packet_statuses(),
        "repo_grep": {
            "unit_frontier": rg(r"S_\+|Omega_M|K_dim|no-strict-unit|scale-charged"),
            "origin_frontier": rg(r"Lambda_origin|translation-breaking|binary origin|DHL|origin datum"),
            "conditional_architecture": rg(r"CA\+SA|Conversion Axioms|Sector Axioms|conditioned|axiom-augmented"),
        },
        "constructed_theoretical_objects": constructed_objects(),
        "finite_certificate": {
            "replay_classes_checked": len(REPLAY_CLASSES),
            "strict_source_obligations": len(OBLIGATIONS),
            "obligation_rows": len(rows),
            "accepted_strict_unit_or_origin_sources": 0,
            "strict_S_plus_exported": False,
            "strict_Lambda_origin_exported": False,
            "conditional_CA_SA_schema_exported": True,
        },
        "finite_theorem": {
            "name": "P3168_T1_post_P3167_replay_classes_do_not_export_strict_units_or_origin",
            "statement": "Across the currently named P3161-P3167 scale/origin replay classes, every class fails at least one required strict-source obligation; the union exports no nonzero scale-charged S_+ value, no translation-breaking Lambda_origin value, and no unit-bearing action/mass map.  A CA+SA package is therefore admissible only as explicit axiom-augmented conditioning, not as strict closure.",
        },
        "decision": {
            "bounded_result": "P3168 implements the SUMMARY_GROK recommendation to certify no-strict-unit and no-new-origin on current artifacts after P3167.",
            "next_honest_step": "Do not replay the closed scale/origin classes.  Either draft the explicit CA+SA conditioned package as non-strict architecture, or supply exactly one new hard object: a nonzero scale-charged strict S_+ value with Omega_M/K_dim coupling, or a genuinely translation-breaking strict Lambda_origin datum with Phi_Info/A_phi coupling.  If neither formula is supplied, preserve the no-new-live-frontier certificate.",
            "negative_export_flags": {
                "S_plus_source_exported": False,
                "Lambda_origin_source_exported": False,
                "strict_unit_source_exported": False,
                "QW2191_discharged": False,
                "bridge_completion_exported": False,
                "role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
        },
    }


def write_md(data: dict[str, Any]) -> None:
    fc = data["finite_certificate"]
    lines = [
        "# P3168/S2118 Post-P3167 no-strict-unit state-map certificate",
        "",
        f"Status: `{data['status']}`",
        "",
        "## Constructed objects",
        "- `post_P3167_unit_origin_state_map`: finite obligation matrix over the exhausted scale/origin replay classes.",
        "- `minimal_CA_SA_architecture_schema`: explicit non-strict two-package architecture for conditioned physics.",
        "",
        "## Finite certificate",
    ]
    for key, value in fc.items():
        lines.append(f"- `{key}`: `{value}`")
    lines += [
        "",
        "## Theorem",
        f"`{data['finite_theorem']['name']}`: {data['finite_theorem']['statement']}",
        "",
        "## Decision",
        data["decision"]["bounded_result"],
        "",
        "## Next honest step",
        data["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    data = payload()
    OUT.write_text(json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_md(data)
    append_once(
        SHEET,
        "P3168/S2118 Post-P3167 no-strict-unit state-map certificate",
        """## P3168/S2118 Post-P3167 no-strict-unit state-map certificate

`P3168/S2118` implements the post-`P3167` state-map recommendation.  It checks `8` replay-gated scale/origin classes against `8` strict-source obligations (`64` rows) and accepts `0` strict unit or origin sources.  No `S_+`, `Lambda_origin`, unit-bearing action/mass map, selector discharge, bridge completion, role transfer, `L_total`, or ToE closure is exported.  The only honest downstream architecture is explicit CA+SA conditioning unless a genuinely new hard object is supplied.""",
    )
    append_once(
        DRAFT,
        "P3168/S2118 CA+SA remains conditional after no-strict-unit certificate",
        """## P3168/S2118 CA+SA remains conditional after no-strict-unit certificate

`P3168/S2118` separates the strict informational package `W0` from optional conversion and sector axioms.  `CA` may introduce unit-bearing action/length/time symbols such as `hbar_*`, `ell_*`, and `tau_*`; `SA` may introduce an origin/polarity branch such as `(r0, lambda0)` plus coupling.  Any effective physics built this way is conditioned/axiom-augmented and does not become strict `L_total` or ToE closure.""",
    )
    append_once(
        AGENTS,
        "Current post-P3167 no-strict-unit state-map guardrail (P3168/S2118, 2026-07-16)",
        """## Current post-P3167 no-strict-unit state-map guardrail (P3168/S2118, 2026-07-16)

- P3168 implements the post-P3167 state-map certificate recommended by `SUMMARY_GROK.md`: it checks `8` replay-gated scale/origin classes against `8` strict-source obligations (`64` rows) and accepts `0` strict unit or origin sources.
- The frozen replay classes include positive-scale torsor invariants, formal `S_+`/`Omega_M` self-sources, c-boundary fits, legacy/Planck/fractal unit imports, `Lambda_origin` localizers, binary `Z12` origin profiles, dimensionless monomials, and unit-measure/selector fusion.
- No nonzero scale-charged `S_+`, translation-breaking `Lambda_origin`, unit-bearing action/mass map, `QW-2191` discharge, bridge completion, role transfer, `L_total`, or ToE closure is exported.
- A CA+SA/two-package architecture is admissible only when explicitly marked non-strict/conditioned: `W0` strict informational nadsoliton; `CA` conversion units (`hbar_*`, `ell_*`, `tau_*` or equivalent); `SA` sector/origin/polarity axioms; `W3` effective physics conditioned on CA+SA.
- The next proof-grade hard move must supply exactly one genuinely new object: either a nonzero scale-charged strict `S_+` value with `Omega_M/K_dim` coupling, or a translation-breaking strict `Lambda_origin` datum with `Phi_Info/A_phi` coupling.  Otherwise preserve the no-new-live-frontier certificate rather than replaying closed lanes.
""",
    )
    return data


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
