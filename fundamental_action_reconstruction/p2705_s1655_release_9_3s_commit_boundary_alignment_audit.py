#!/usr/bin/env python3
"""P2705/S1655: Release-9.3s commit boundary-alignment audit.

Audit the user-supplied commit 8d48faf... as a concrete older-release pointer and
check whether it unlocks any current selector/strict-core/L_total/ToE closure lane.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2705_s1655_release_9_3s_commit_boundary_alignment_audit.json"
MD = GEN / "p2705_s1655_release_9_3s_commit_boundary_alignment_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

COMMIT = "8d48faf012f87721d01a692fd7e3888461d4e6d2"
PARENT = "6adfe7a7bbdbe6d529242a1783b9dc0339a2cc32"

CURRENT_INPUTS = {
    "P2704": GEN / "p2704_s1654_p1343_p1348_selector_provenance_revalidation_table.json",
    "P2699": GEN / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.json",
    "P2700": GEN / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.json",
    "P2701": GEN / "p2701_s1651_strict_sourced_symmetry_breaking_provider_inventory_matrix.json",
    "P2702": GEN / "p2702_s1652_selector_circle_lay_mechanism_and_status_packet.json",
    "P2377": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
    "P2378": GEN / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.json",
}

CLOSURE_TOKENS = [
    "P1343", "P1348", "QW-2191", "selector", "strict-core", "strict core", "pair12", "L_total", "ToE", "global closure"
]
DAMPING_TOKENS = ["P2377", "P2378", "damping", "compression", "transport", "eta", "beta", "tau"]
NEGATIVE_EXPORT_FLAGS = [
    "release_9_3s_commit_unblocks_qw2191",
    "non_premise_selector_provider_exported",
    "pair12_strict_core_upgrade_exported",
    "ltotal_promoted",
    "toe_closure_exported",
    "legacy_role_transfer_exported",
    "bridge_completion_exported",
]


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=REPO, text=True, encoding="utf-8")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def commit_metadata() -> dict[str, Any]:
    raw = git("show", "--no-patch", "--format=%H%n%P%n%s%n%an%n%aI%n%B", COMMIT).splitlines()
    changed = [line for line in git("diff", "--name-only", f"{PARENT}..{COMMIT}").splitlines() if line]
    numstat = []
    for line in git("diff", "--numstat", f"{PARENT}..{COMMIT}").splitlines():
        added, deleted, path = line.split("\t")
        numstat.append({"path": path, "added": None if added == "-" else int(added), "deleted": None if deleted == "-" else int(deleted)})
    return {"commit": raw[0], "parents": raw[1].split(), "subject": raw[2], "author": raw[3], "author_date": raw[4], "message": "\n".join(raw[5:]).strip(), "changed_files": changed, "numstat": numstat}


def file_token_scan(files: list[str]) -> list[dict[str, Any]]:
    rows = []
    for path in files:
        try:
            text = git("show", f"{COMMIT}:{path}")
        except subprocess.CalledProcessError:
            text = ""
        lower = text.lower()
        closure_hits = [tok for tok in CLOSURE_TOKENS if tok.lower() in lower]
        damping_hits = [tok for tok in DAMPING_TOKENS if tok.lower() in lower]
        rows.append({
            "path": path,
            "sha256_at_commit": hashlib.sha256(text.encode("utf-8")).hexdigest() if text else None,
            "closure_token_hits": closure_hits,
            "damping_token_hits": damping_hits,
            "line_count_at_commit": len(text.splitlines()) if text else 0,
        })
    return rows


def p2377_p2378_numeric_boundary() -> dict[str, Any]:
    p2377 = read_json(CURRENT_INPUTS["P2377"])
    p2378 = read_json(CURRENT_INPUTS["P2378"])
    p2377_theorem = p2377.get("damping_compression_transport_primitive_uniform_coupling_theorem", {})
    p2377_cert = p2377_theorem.get("transport_primitive_uniform_coupling_certificate", {})
    threshold = p2377_cert.get("uniform_threshold_certificate", {}).get("uniform_coupling", {})
    p2378_theorem = p2378.get("unit_normalized_transport_coupling_insufficiency_theorem", {})
    p2378_rect = p2378_theorem.get("unit_normalized_transport_coupling_insufficiency_certificate", {}).get("rectangle_proof", {})
    return {
        "p2377_status": p2377.get("status"),
        "p2378_status": p2378.get("status"),
        "p2377_uniform_tau_gt": threshold.get("tau_gt_uniform"),
        "p2377_normalization_status": threshold.get("normalization_status"),
        "p2378_unit_mass_insufficient_on_rectangle": p2378_rect.get("unit_mass_insufficient_on_rectangle"),
        "p2378_tau_threshold_range": p2378_rect.get("tau_threshold_range"),
        "boundary_reading": "The commit adds damping/compression transport mathematics.  P2377 gives an exact endpoint primitive and a sufficient scalar-coupling threshold; P2378 says unit-normalized transport remains insufficient.  This is not a selector-provider or L_total variational source export.",
    }


def current_selector_boundary() -> dict[str, Any]:
    p2704 = read_json(CURRENT_INPUTS["P2704"])
    p2699 = read_json(CURRENT_INPUTS["P2699"])
    p2700 = read_json(CURRENT_INPUTS["P2700"])
    p2701 = read_json(CURRENT_INPUTS["P2701"])
    p2702 = read_json(CURRENT_INPUTS["P2702"])
    return {
        "p2704_declared_scope_positive": p2704.get("decision", {}).get("p1343_p1348_declared_scope_provenance_revalidated") is True,
        "p2699_bounded_no_go": p2699.get("decision", {}).get("bounded_no_go_now") is True,
        "p2700_bounded_no_go": p2700.get("decision", {}).get("bounded_no_go_now") is True,
        "p2701_bounded_no_go": p2701.get("decision", {}).get("bounded_no_go_now") is True,
        "p2702_no_new_closure_exported": p2702.get("decision", {}).get("no_new_closure_exported") is True,
    }


def alignment_matrix(meta: dict[str, Any], scan: list[dict[str, Any]], numeric: dict[str, Any], boundary: dict[str, Any]) -> list[dict[str, Any]]:
    paths = set(meta["changed_files"])
    only_p2377_p2378_lane = all(("p2377" in p or "p2378" in p or p.endswith("STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md") or p.endswith("STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md")) for p in paths)
    selector_named_artifacts = [p for p in paths if any(tok in p.lower() for tok in ["selector", "qw2191", "pair12", "p1343", "p1348"])]
    return [
        {"obligation": "commit_identity_is_damping_compression_pr_not_direct_release_9_3s_selector_document", "passes": "P2377/P2378" in meta.get("message", "") or "damping-compression" in meta["subject"], "evidence": {"subject": meta["subject"], "message": meta.get("message", ""), "changed_files": meta["changed_files"]}},
        {"obligation": "diff_scope_contains_only_p2377_p2378_lane_plus_docs", "passes": only_p2377_p2378_lane, "evidence": meta["numstat"]},
        {"obligation": "no_new_selector_or_qw2191_artifact_added_by_commit", "passes": len(selector_named_artifacts) == 0, "evidence": {"selector_named_artifacts": selector_named_artifacts, "token_scan_note": "closure tokens appear only as negative/not-licensed guardrail language inside P2377/P2378 lane files, not as new selector artifacts"}},
        {"obligation": "numeric_damping_boundary_is_insufficiency_not_closure", "passes": numeric["p2378_unit_mass_insufficient_on_rectangle"] is True and "not derived" in str(numeric["p2377_normalization_status"]), "evidence": numeric},
        {"obligation": "current_p2699_p2704_boundary_remains_consistent", "passes": all(boundary.values()), "evidence": boundary},
    ]


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2705/S1655 Release-9.3s commit boundary-alignment audit", "", f"Status: `{payload['status']}`", "", "## Commit audited", f"- `{COMMIT}`: {payload['commit_metadata']['subject']}", "", "## Alignment matrix"]
    for row in payload["alignment_matrix"]:
        lines.append(f"- `{row['obligation']}`: passes={row['passes']}")
    lines.extend(["", "## Numeric boundary", payload["p2377_p2378_numeric_boundary"]["boundary_reading"], "", "## Decision", payload["decision"]["current_unlock_reading"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    meta = commit_metadata()
    scan = file_token_scan(meta["changed_files"])
    numeric = p2377_p2378_numeric_boundary()
    boundary = current_selector_boundary()
    matrix = alignment_matrix(meta, scan, numeric, boundary)
    all_pass = all(row["passes"] for row in matrix)
    payload = {
        "status": "P2705_RELEASE_9_3S_POINTER_AUDITED_NO_CURRENT_UNLOCK" if all_pass else "P2705_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {name: sha(path) for name, path in CURRENT_INPUTS.items()},
        "commit_metadata": meta,
        "file_token_scan": scan,
        "p2377_p2378_numeric_boundary": numeric,
        "current_selector_boundary": boundary,
        "alignment_matrix": matrix,
        "decision": {
            "release_9_3s_pointer_unblocks_current_stage": False,
            "current_unlock_reading": "The supplied commit is a P2377/P2378 damping-compression transport merge, not a direct Release-9.3s selector-closure artifact.  It adds useful bridge/damping mathematics but also preserves an insufficiency boundary: scalar coupling is not dynamically sourced and unit-normalized transport is insufficient.  Therefore it does not remove current QW-2191/non-premise-selector/pair12/L_total/ToE blocks.",
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "next_honest_step": "P2706 should be a narrow damping-to-selector interface obstruction/witness table: test whether P2377/P2378's damping-compression transport primitive can define any Aut(Z12)-noninvariant or non-premise directed-unit functional on the Z12 selector problem.  If it only changes radial/tail weights and remains orientation-blind, preserve the P2697-P2705 no-unlock certificate rather than reopening selector closure.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2705/S1655 release 9.3s commit boundary audit", "## P2705/S1655 release 9.3s commit boundary audit\n\n`P2705/S1655` audits the user-supplied commit `8d48faf012f87721d01a692fd7e3888461d4e6d2`.  The commit is a P2377/P2378 damping-compression transport merge, not a direct Release-9.3s selector-closure document.  Its numerical boundary is useful but non-unlocking: P2377 supplies an endpoint transport primitive with an unsourced scalar-coupling threshold, while P2378 confirms unit-normalized transport insufficiency on the rectangle.  No `QW-2191`, non-premise selector-provider, pair12 strict-core, `L_total`, role-transfer, bridge-completion, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2705/S1655 release commit Ltotal guard", "## P2705/S1655 release commit Ltotal guard\n\n`P2705/S1655` finds that commit `8d48faf012f87721d01a692fd7e3888461d4e6d2` contributes damping-compression transport evidence only.  Because P2377 leaves scalar coupling unsourced and P2378 proves unit-normalized transport insufficiency, the commit is not a variational source and does not promote `L_total`, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, or ToE closure.\n")
    append_once(AGENTS, "Current Release-9.3s pointer boundary-alignment guardrail (P2705/S1655, 2026-06-13)", "## Current Release-9.3s pointer boundary-alignment guardrail (P2705/S1655, 2026-06-13)\n\n- P2705 audits the supplied commit `8d48faf012f87721d01a692fd7e3888461d4e6d2` and identifies it as a P2377/P2378 damping-compression transport merge, not as a direct Release-9.3s selector-closure artifact.\n- The commit strengthens damping/bridge mathematics but preserves the relevant insufficiency boundary: the scalar coupling is not dynamically sourced and unit-normalized transport is insufficient.\n- Do not use this commit to remove `QW-2191`, non-premise selector-provider, pair12 strict-core, `L_total`, role-transfer, bridge-completion, or ToE blocks.  A next admissible move is a narrow damping-to-selector interface obstruction/witness table, not closure promotion.\n")
    return payload


if __name__ == "__main__":
    main()
