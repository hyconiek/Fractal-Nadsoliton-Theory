#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from collections.abc import Callable
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2646_s1596_frozen_kernel_compression_signature_preregistration_certificate.json"
MD = GEN / "p2646_s1596_frozen_kernel_compression_signature_preregistration_certificate.md"

BETA_TORS = 0.01
STRICT_BETA = 1.0
ETA = 9.0 / 5.0
TEST_PAIRS = [(1, 7), (1, 12), (2, 8), (3, 9)]

SOURCE_FILES = {
    "P2635_EMPIRICAL_INTERFACE": GEN / "p2635_s1585_toe_neural_universe_empirical_signature_audit.json",
    "P2643_BETA_THRESHOLD": GEN / "p2643_s1593_inverse_hierarchy_beta_threshold_role_rejection_certificate.json",
    "P2644_COMPRESSED_SUCCESSOR": GEN / "p2644_s1594_modified_compressed_inverse_hierarchy_successor_theorem.json",
    "P2645_ROLE_MATRIX": GEN / "p2645_s1595_role_transfer_matrix_and_closure_route_rerun.json",
    "STRICT_EQUATION_SHEET": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "STRICT_LAGRANGIAN_DRAFT": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "empirical_confirmation_claimed",
    "data_fit_performed",
    "kernel_retuning_allowed",
    "beta_source_exported",
    "unchanged_inverse_hierarchy_restored",
    "legacy_role_transfer_revalidated",
    "bridge_completion_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "*.lean",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:45]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "blind_frozen_empirical_content": (
            "blind frozen-kernel|frozen-kernel empirical|holdout|preregistration|pre-registered|"
            "kernel retune|without retuning|adversarial controls"
        ),
        "compression_signature_content": (
            "compression signature|compression/locality-bias|locality-bias|strict denominator|"
            "heavy-tailed compression|attention suppression|suppression factor|tail ratio"
        ),
        "phase_amplitude_invariant_content": (
            "phase-blind|phase masked|amplitude.*cancel|denominator-only|envelope|tail slope|"
            "cosine.*factor|positional encoding"
        ),
        "legacy_strict_discriminator_content": (
            "legacy hyperbolic|beta_tors|inverse hierarchy|distant octave|unchanged inverse-hierarchy|"
            "strict beta=1|eta=9/5"
        ),
        "guardrail_nonclosure_content": (
            "beta source|target-independent|role-transfer|role-bearing L_total|QW-2191|ToE closure|"
            "empirical confirmation"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for frozen-kernel compression-signature preregistration, not packet-name search",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def latest_commit_audit() -> list[dict[str, Any]]:
    proc = subprocess.run(["git", "log", "-n", "10", "--oneline", "--name-only"], cwd=REPO, text=True, stdout=subprocess.PIPE, check=True)
    rows: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        if re.match(r"^[0-9a-f]{7,12} ", line):
            if current:
                rows.append(current)
            sha, subject = line.split(" ", 1)
            current = {"sha": sha, "subject": subject, "files": []}
        elif current is not None:
            current["files"].append(line)
    if current:
        rows.append(current)
    return rows


def strict_attention(d: float) -> float:
    return 1.0 / (1.0 + STRICT_BETA * d**ETA)


def legacy_attention(d: float) -> float:
    return 1.0 / (1.0 + BETA_TORS * d)


def tail_ratio(attention: Callable[[float], float], near: float, far: float) -> float:
    return attention(far) / attention(near)


def log_tail_slope(ratio: float, near: float, far: float) -> float:
    return -math.log(ratio) / math.log(far / near)


def frozen_discriminator_rows() -> list[dict[str, Any]]:
    rows = []
    for near, far in TEST_PAIRS:
        strict_ratio = tail_ratio(strict_attention, near, far)
        legacy_ratio = tail_ratio(legacy_attention, near, far)
        strict_slope = log_tail_slope(strict_ratio, near, far)
        legacy_slope = log_tail_slope(legacy_ratio, near, far)
        rows.append({
            "near": near,
            "far": far,
            "strict_denominator_tail_ratio": strict_ratio,
            "legacy_hyperbolic_tail_ratio": legacy_ratio,
            "strict_over_legacy_tail_ratio": strict_ratio / legacy_ratio,
            "midpoint_ratio_threshold": 0.5 * (strict_ratio + legacy_ratio),
            "strict_log_tail_slope": strict_slope,
            "legacy_log_tail_slope": legacy_slope,
            "midpoint_slope_threshold": 0.5 * (strict_slope + legacy_slope),
            "ratio_gap": legacy_ratio - strict_ratio,
            "slope_gap": strict_slope - legacy_slope,
        })
    return rows


def preregistered_tests(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "locked_parameters": {
            "legacy_beta_tors": BETA_TORS,
            "strict_beta": STRICT_BETA,
            "strict_eta": ETA,
            "phase_and_amplitude_policy": (
                "Use denominator/envelope observables or phase-masked bins so alpha and cos(omega*d+phi) do not decide the test; "
                "omega, phi, amplitude, beta, and eta are not fitted on the holdout."
            ),
        },
        "observable_family": "phase/amplitude-invariant denominator tail ratios R(a,b)=A(b)/A(a) and log-tail slopes sigma(a,b)=-log R/log(b/a)",
        "pass_rule_for_strict_compression_signature": (
            "On blind holdout data, the measured denominator/envelope tail ratio must lie below the preregistered midpoint threshold "
            "and the measured log-tail slope must lie above the preregistered midpoint threshold for the audited pairs, with no kernel retuning."
        ),
        "failure_rule": (
            "Fail the strict compression-signature claim if ratios are legacy-like, if thresholds pass only after refitting beta/eta/phase per dataset, "
            "or if a matched exponential/spline/control baseline removes the frozen-kernel advantage."
        ),
        "rows": rows,
        "all_ratio_gaps_positive": all(row["ratio_gap"] > 0 for row in rows),
        "all_slope_gaps_positive": all(row["slope_gap"] > 0 for row in rows),
        "minimum_ratio_gap": min(row["ratio_gap"] for row in rows),
        "minimum_slope_gap": min(row["slope_gap"] for row in rows),
    }


def upstream_consistency() -> dict[str, Any]:
    p2644 = load_json(SOURCE_FILES["P2644_COMPRESSED_SUCCESSOR"])
    p2643 = load_json(SOURCE_FILES["P2643_BETA_THRESHOLD"])
    p2635 = load_json(SOURCE_FILES["P2635_EMPIRICAL_INTERFACE"])
    p2645 = load_json(SOURCE_FILES["P2645_ROLE_MATRIX"])
    return {
        "p2644_modified_successor_not_full_transfer": "UNCHANGED_INVERSE_HIERARCHY_REJECTED" in json.dumps(p2644, sort_keys=True),
        "p2643_threshold_rejection_available": "beta_crit" in json.dumps(p2643, sort_keys=True) or "threshold" in json.dumps(p2643, sort_keys=True),
        "p2635_empirical_interface_available": "holdout" in json.dumps(p2635, sort_keys=True).lower(),
        "p2645_role_matrix_points_to_holdout": "blind frozen-kernel empirical compression test" in json.dumps(p2645, sort_keys=True),
    }


def closure_decision(tests: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "PREREGISTERED_FROZEN_COMPRESSION_SIGNATURE_INTERFACE_ONLY__NO_EMPIRICAL_CONFIRMATION_NO_SOURCE_EXPORT",
        "professorial_verdict": (
            "P2646 converts the P2644 compression/locality-bias reading and P2645 role-matrix route into a falsifiable, "
            "phase/amplitude-invariant holdout protocol. The exact strict-vs-legacy denominator gaps are large on the audited pairs, "
            "so the next empirical check is no longer a vague analogy: it is a locked tail-ratio/log-slope discriminator. But this is a preregistration certificate only; "
            "without blind data it does not prove beta sourcehood, role transfer, L_total, or ToE closure."
        ),
        "professorial_closure_path": [
            "Run the blind frozen-kernel compression-signature holdout exactly as locked here, including phase masks and control baselines.",
            "In parallel, continue the target-independent beta-source theorem; empirical success cannot replace the source theorem but can prioritize it.",
            "If the holdout fails or requires retuning, demote strict compression to a useful model-level successor and do not reopen role-bearing L_total.",
            "If the holdout passes, rerun the role-transfer matrix only for the modified/compressed successor semantics, not for unchanged inverse-hierarchy.",
            "Only after beta source, selector/source QW-2191, role-transfer matrix, and empirical controls pass can frozen-kernel L_total be reconsidered.",
        ],
        "next_honest_step": (
            "Do not add more node/offset/stride or amplitude rescues. The next honest step is either an actual blind holdout execution of the P2646 tail-ratio/log-slope compression signature, "
            "or a target-independent beta-source theorem; the strongest route is to run the holdout while separately proving the beta source."
        ),
        "can_update_role_transfer_table": tests["all_ratio_gaps_positive"] and tests["all_slope_gaps_positive"],
        "full_kernel_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    tests = payload["preregistered_frozen_kernel_tests"]
    decision = payload["closure_decision"]
    lines = [
        "# P2646/S1596 frozen-kernel compression-signature preregistration certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps blind/frozen empirical content, compression signatures, phase/amplitude-invariant observables, legacy-vs-strict discriminators, and nonclosure guardrails before adding the certificate.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Locked discriminator",
        "",
        f"Observable family: `{tests['observable_family']}`.",
        tests["pass_rule_for_strict_compression_signature"],
        tests["failure_rule"],
        "",
        "| near | far | strict tail ratio | legacy tail ratio | midpoint ratio threshold | strict slope | legacy slope | midpoint slope threshold |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in tests["rows"]:
        lines.append(
            f"| {row['near']} | {row['far']} | `{row['strict_denominator_tail_ratio']:.12f}` | "
            f"`{row['legacy_hyperbolic_tail_ratio']:.12f}` | `{row['midpoint_ratio_threshold']:.12f}` | "
            f"`{row['strict_log_tail_slope']:.12f}` | `{row['legacy_log_tail_slope']:.12f}` | `{row['midpoint_slope_threshold']:.12f}` |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Full kernel now? `{decision['full_kernel_now']}`.",
        f"ToE closure now? `{decision['toe_closure_now']}`.",
        "",
        "## Professorial closure path",
        "",
    ])
    for item in decision["professorial_closure_path"]:
        lines.append(f"- {item}")
    lines.extend([
        "",
        "## Next honest step",
        "",
        decision["next_honest_step"],
        "",
        "## Negative exports",
        "",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    rows = frozen_discriminator_rows()
    tests = preregistered_tests(rows)
    decision = closure_decision(tests)
    payload: dict[str, Any] = {
        "status": "P2646_FROZEN_KERNEL_COMPRESSION_SIGNATURE_PREREGISTRATION_CERTIFICATE_NO_EMPIRICAL_CONFIRMATION_NO_FALSE_PASS",
        "script_correction_note": (
            "The user-supplied draft was mathematically usable, but it was pasted twice and reused P2645/S1595, "
            "which is already occupied by the role-transfer matrix. This implementation keeps the proof content and promotes it to the next packet P2646/S1596."
        ),
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {name: sha256_file(path) for name, path in SOURCE_FILES.items()},
        "upstream_consistency": upstream_consistency(),
        "preregistered_frozen_kernel_tests": tests,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        SOURCE_FILES["STRICT_EQUATION_SHEET"],
        "P2646/S1596 frozen-kernel compression-signature preregistration guard",
        "## P2646/S1596 frozen-kernel compression-signature preregistration guard\n\n"
        "`P2646/S1596` turns the P2644 compression/locality-bias successor and P2645 role-matrix route into a locked empirical discriminator rather than a new fit: phase/amplitude-invariant denominator tail ratios and log-tail slopes are preregistered for audited pairs such as `(1,7)` and `(1,12)`, with strict `beta=1, eta=9/5` predictions separated far below the legacy `beta_tors=0.01` hyperbolic tail.  This is only a blind-holdout interface; it exports no empirical confirmation, beta source, bridge completion, role transfer, `QW-2191` discharge, role-bearing `L_total`, or ToE closure.\n",
    )
    append_once(
        SOURCE_FILES["STRICT_LAGRANGIAN_DRAFT"],
        "P2646/S1596 frozen-kernel compression-signature Ltotal guard",
        "## P2646/S1596 frozen-kernel compression-signature Ltotal guard\n\n"
        "`P2646/S1596` does not re-enable `L_total`: it preregisters a no-retuning empirical tail-ratio/log-slope test for the modified compression successor, but a role-bearing variational term still needs blind data plus a target-independent beta/source theorem and the later role-transfer matrix.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
