"""P3169/S2119: ternary origin-datum exhaustive audit.

This is the next hard-object test after P3168.  Instead of another summary
certificate, it constructs a genuinely new finite candidate class for the
missing Lambda_origin object: all nonzero ternary signed Z12 profiles with
values {-1,0,+1}.  The class is outside the already closed binary profile class
(P3166) and is audited as a possible translation-breaking origin/polarity datum
coupled to Phi_Info/A_phi.

The finite result is deliberately conservative: many ternary receivers have
trivial translation stabilizer, nonzero resultant, and unique extrema, but the
profile class itself supplies no strict nadsoliton provenance law selecting one
absolute representative.  Translation quotienting still gives only orbits, and
canonical representatives import an ordering convention.  Therefore no strict
Lambda_origin source is exported.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
GEN.mkdir(exist_ok=True)
OUT = GEN / "p3169_s2119_ternary_origin_datum_exhaustive_audit.json"
MD = GEN / "p3169_s2119_ternary_origin_datum_exhaustive_audit.md"
AGENTS = REPO / "AGENTS.md"
SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

INPUTS = {
    "P3165": GEN / "p3165_s2115_lambda_origin_source_localizer_audit.json",
    "P3166": GEN / "p3166_s2116_binary_origin_datum_exhaustive_intake.json",
    "P3168": GEN / "p3168_s2118_post_p3167_no_strict_unit_state_map_certificate.json",
    "SUMMARY_GROK": REPO / "SUMMARY_GROK.md",
}
VALUES = (-1, 0, 1)
N = 12
ROOTS = [complex(math.cos(2 * math.pi * k / N), math.sin(2 * math.pi * k / N)) for k in range(N)]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    except Exception:
        return {}


def rotate(profile: tuple[int, ...], shift: int) -> tuple[int, ...]:
    shift %= N
    return profile[shift:] + profile[:shift]


def affine(profile: tuple[int, ...], unit: int, shift: int) -> tuple[int, ...]:
    out = [0] * N
    for i, value in enumerate(profile):
        out[(unit * i + shift) % N] = value
    return tuple(out)


def translation_orbit(profile: tuple[int, ...]) -> list[tuple[int, ...]]:
    return [rotate(profile, s) for s in range(N)]


def canonical_translation(profile: tuple[int, ...]) -> tuple[int, ...]:
    return min(translation_orbit(profile))


def canonical_affine(profile: tuple[int, ...]) -> tuple[int, ...]:
    return min(affine(profile, u, s) for u in (1, 5, 7, 11) for s in range(N))


def stabilizer_size(profile: tuple[int, ...]) -> int:
    return sum(1 for s in range(N) if rotate(profile, s) == profile)


def resultant(profile: tuple[int, ...]) -> complex:
    return sum(value * ROOTS[i] for i, value in enumerate(profile))


def unique_extrema(profile: tuple[int, ...]) -> dict[str, Any]:
    max_value = max(profile)
    min_value = min(profile)
    max_positions = [i for i, value in enumerate(profile) if value == max_value]
    min_positions = [i for i, value in enumerate(profile) if value == min_value]
    return {
        "unique_max": len(max_positions) == 1,
        "unique_min": len(min_positions) == 1,
        "max_positions": max_positions[:4],
        "min_positions": min_positions[:4],
    }


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


def finite_scan() -> dict[str, Any]:
    translation_classes: set[tuple[int, ...]] = set()
    affine_classes: set[tuple[int, ...]] = set()
    stabilizers = Counter()
    sum_counts = Counter()
    support_counts = Counter()
    unique_max = 0
    unique_min = 0
    trivial_stabilizer = 0
    nonzero_resultant = 0
    chiral_pairs = 0
    representative_samples: list[dict[str, Any]] = []

    total = 0
    for profile in itertools.product(VALUES, repeat=N):
        if all(value == 0 for value in profile):
            continue
        total += 1
        translation_classes.add(canonical_translation(profile))
        affine_classes.add(canonical_affine(profile))
        stab = stabilizer_size(profile)
        stabilizers[stab] += 1
        if stab == 1:
            trivial_stabilizer += 1
        s = sum(profile)
        sum_counts[s] += 1
        supp = sum(1 for value in profile if value != 0)
        support_counts[supp] += 1
        ext = unique_extrema(profile)
        if ext["unique_max"]:
            unique_max += 1
        if ext["unique_min"]:
            unique_min += 1
        res = resultant(profile)
        if abs(res) > 1e-12:
            nonzero_resultant += 1
        neg = tuple(-value for value in profile)
        if profile < neg:
            chiral_pairs += 1
        if len(representative_samples) < 12 and stab == 1 and abs(res) > 1e-12 and (ext["unique_max"] or ext["unique_min"]):
            representative_samples.append(
                {
                    "profile": profile,
                    "sum": s,
                    "support": supp,
                    "resultant_re": round(res.real, 12),
                    "resultant_im": round(res.imag, 12),
                    "translation_orbit_size": N // stab,
                    "unique_max": ext["unique_max"],
                    "unique_min": ext["unique_min"],
                    "accepted_Lambda_origin_source": False,
                    "blocker": "receiver has no strict provenance law selecting this absolute representative",
                }
            )
    return {
        "ternary_profiles_nonzero": total,
        "translation_classes": len(translation_classes),
        "affine_translation_classes": len(affine_classes),
        "stabilizer_histogram": dict(sorted(stabilizers.items())),
        "support_histogram": dict(sorted(support_counts.items())),
        "sum_histogram_excerpt": {str(k): sum_counts[k] for k in sorted(sum_counts)[:5] + sorted(sum_counts)[-5:]},
        "trivial_translation_stabilizer_profiles": trivial_stabilizer,
        "nonzero_first_resultant_profiles": nonzero_resultant,
        "unique_max_profiles": unique_max,
        "unique_min_profiles": unique_min,
        "negation_pairs": chiral_pairs,
        "representative_receiver_samples": representative_samples,
    }


def gate_rows(scan: dict[str, Any]) -> list[dict[str, Any]]:
    gates = {
        "new_nonbinary_candidate_class_constructed": True,
        "finite_exhaustion_completed": True,
        "translation_breaking_receivers_exist": scan["trivial_translation_stabilizer_profiles"] > 0,
        "phase_resultant_receivers_exist": scan["nonzero_first_resultant_profiles"] > 0,
        "unique_extremum_receivers_exist": scan["unique_max_profiles"] > 0 or scan["unique_min_profiles"] > 0,
        "strict_nadsoliton_provenance_law_exported": False,
        "absolute_origin_representative_without_convention": False,
        "Phi_Info_A_phi_coupling_theorem_exported": False,
        "selector_safety_no_QW2191_discharge": True,
        "accepted_Lambda_origin_source": False,
    }
    return [{"candidate_class": "ternary_signed_Z12_profiles", "gate": gate, "passed": passed} for gate, passed in gates.items()]


def payload() -> dict[str, Any]:
    scan = finite_scan()
    gates = gate_rows(scan)
    return {
        "status": "P3169_TERNARY_ORIGIN_DATUM_EXHAUSTIVE_AUDIT_BOUNDED_NO_GO",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: load(path).get("status") or ("markdown_or_missing_json" if path.exists() else "missing") for key, path in INPUTS.items()},
        "repo_grep": {
            "origin_frontier": rg(r"Lambda_origin|translation-breaking|origin datum|Phi_Info|A_phi"),
            "closed_binary_or_DHL": rg(r"binary origin|binary `Z12`|DHL|receiver-to-source|P3166|P3139"),
            "scale_alternative": rg(r"S_\+|Omega_M|K_dim|scale-charged|no-strict-unit"),
        },
        "constructed_theoretical_objects": {
            "ternary_signed_Z12_profile_family": "all nonzero profiles f: Z12 -> {-1,0,+1}",
            "translation_quotient": "orbits under cyclic shifts of source labels",
            "affine_quotient": "orbits under translations and Aut(Z12)=U(12)",
            "first_phase_resultant_receiver": "R1(f)=sum_i f_i exp(2*pi*i*i/12), a receiver coupled in shape to Phi_Info/A_phi but not a source law",
            "gate_rows": gates,
            "representative_receiver_samples": scan["representative_receiver_samples"],
        },
        "finite_certificate": {
            **{key: value for key, value in scan.items() if key != "representative_receiver_samples"},
            "gate_rows": len(gates),
            "passed_receiver_gates": sum(1 for row in gates if row["passed"]),
            "accepted_Lambda_origin_sources": 0,
        },
        "finite_theorem": {
            "name": "P3169_T1_ternary_receivers_do_not_supply_absolute_strict_origin",
            "statement": "The exhaustive ternary signed Z12 profile family contains many translation-breaking and phase-resultant receivers, but translation/affine quotienting selects only orbits and no current artifact exports a strict nadsoliton provenance law or Phi_Info/A_phi coupling theorem selecting one absolute representative.  Therefore the family exports no strict Lambda_origin source and does not discharge QW-2191.",
        },
        "decision": {
            "bounded_result": "P3169 tests one genuinely new non-binary origin candidate class after P3168 and closes it as receiver-only on current artifacts.",
            "next_honest_step": "Do not escalate from binary to ternary/k-ary profile inventories as origin closure.  The next hard move should pivot to the other P3168 branch: construct one nonzero scale-charged strict S_+ value with Omega_M/K_dim coupling.  If no such formula is supplied, preserve the P3168-P3169 no-new-live-frontier certificate or draft CA+SA only as explicit non-strict conditioning.",
            "negative_export_flags": {
                "Lambda_origin_source_exported": False,
                "Phi_Info_A_phi_coupling_exported": False,
                "S_plus_source_exported": False,
                "unit_source_exported": False,
                "QW2191_discharged": False,
                "selector_closure_exported": False,
                "bridge_completion_exported": False,
                "role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
        },
    }


def write_md(data: dict[str, Any]) -> None:
    cert = data["finite_certificate"]
    lines = [
        "# P3169/S2119 Ternary origin-datum exhaustive audit",
        "",
        f"Status: `{data['status']}`",
        "",
        "## Constructed objects",
        "- `ternary_signed_Z12_profile_family`: all nonzero maps `Z12 -> {-1,0,+1}`.",
        "- Translation and affine quotient certificates for this new non-binary origin candidate class.",
        "- `R1(f)=sum_i f_i exp(2*pi*i*i/12)` phase-resultant receiver, treated only as a receiver coupled in shape to `Phi_Info/A_phi`.",
        "",
        "## Finite certificate",
    ]
    for key in [
        "ternary_profiles_nonzero",
        "translation_classes",
        "affine_translation_classes",
        "trivial_translation_stabilizer_profiles",
        "nonzero_first_resultant_profiles",
        "unique_max_profiles",
        "unique_min_profiles",
        "gate_rows",
        "accepted_Lambda_origin_sources",
    ]:
        lines.append(f"- `{key}`: `{cert[key]}`")
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
        "P3169/S2119 ternary origin-datum exhaustive audit",
        """## P3169/S2119 ternary origin-datum exhaustive audit

`P3169/S2119` constructs the non-binary hard-object branch left by `P3168`: all nonzero ternary signed `Z12 -> {-1,0,+1}` profiles as candidate `Lambda_origin` data.  The exhaustive scan finds `531440` profiles, `44367` translation classes, `12768` affine classes, many translation-breaking/phase-resultant/unique-extremum receivers, but `0` accepted strict `Lambda_origin` sources because no strict nadsoliton provenance law or `Phi_Info/A_phi` coupling theorem selects an absolute representative without convention.  No selector discharge, bridge completion, role transfer, `L_total`, or ToE closure is exported.""",
    )
    append_once(
        DRAFT,
        "P3169/S2119 ternary origin receivers do not add Lagrangian source",
        """## P3169/S2119 ternary origin receivers do not add Lagrangian source

`P3169/S2119` proves that the new ternary signed origin-profile family is receiver-rich but source-empty on current artifacts.  Since no strict `Lambda_origin` law or `Phi_Info/A_phi` coupling theorem is exported, no origin-localized action density, unit-bearing measure, EOM term, selector closure, or `L_total` term is added.""",
    )
    append_once(
        AGENTS,
        "Current ternary origin-datum exhaustive audit guardrail (P3169/S2119, 2026-07-16)",
        """## Current ternary origin-datum exhaustive audit guardrail (P3169/S2119, 2026-07-16)

- P3169 tests the non-binary hard-origin branch left by P3168: all `531440` nonzero ternary signed profiles `Z12 -> {-1,0,+1}` as candidate `Lambda_origin` data.
- The scan finds `44367` translation classes and `12768` affine classes, with many translation-breaking, nonzero phase-resultant, and unique-extremum receivers, but `0` accepted strict `Lambda_origin` sources.
- Translation/affine quotienting selects only orbits, and representative choices import conventions unless a strict nadsoliton provenance law plus `Phi_Info/A_phi` coupling theorem is exported.
- Do not escalate binary profiles into ternary/k-ary profile inventories as origin, selector, unit, bridge, role-transfer, `L_total`, or ToE closure.
- The next hard move should pivot to the other P3168 branch: one nonzero scale-charged strict `S_+` value with `Omega_M/K_dim` coupling.  Without such a formula, preserve the P3168-P3169 no-new-live-frontier certificate or draft CA+SA only as explicit non-strict conditioning.
""",
    )
    return data


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
