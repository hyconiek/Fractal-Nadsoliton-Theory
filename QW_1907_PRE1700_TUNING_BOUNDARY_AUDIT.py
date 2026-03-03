#!/usr/bin/env python3
"""
QW-1907: PRE-1700 TUNING BOUNDARY AUDIT
---------------------------------------
Goal:
1) Separate "parameter inference / fitting" from "kernel core retuning".
2) Audit QW-700..1699 files (py/md) with static markers.
3) Produce a reproducible report with explicit caveats.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1907_pre1700_tuning_boundary_audit.json"
OUT_MD = ROOT / "RAPORT_QW1907_PRE1700_TUNING_BOUNDARY_AUDIT.md"

QW_MIN = 700
QW_MAX = 1699
EXTENSIONS = {".py", ".md"}

RX_QW = re.compile(r"(?:QW[_-]?|qw[_-]?)(\d{3,4})")
RX_INFERENCE = re.compile(
    r"\b(minimize|minimizer|curve_fit|least_squares|optuna|objective|"
    r"log_likelihood|likelihood|posterior|bayes|mcmc|sampler|inference)\b",
    re.IGNORECASE,
)
RX_SIMULATION = re.compile(
    r"\b(simulat|simulation|mock|injection|synthetic|generate_catalog|proof of concept)\b",
    re.IGNORECASE,
)
RX_EXTERNAL = re.compile(
    r"\b(raw_strain|nanograv|nano15|gwosc|hdf5|\.h5|ligo|virgo|read_csv|loadtxt)\b",
    re.IGNORECASE,
)
RX_EXTERNAL_STRICT = re.compile(
    r"(raw_strain|nano15|nanograv|gwosc|hdf5|\.h5|/kaggle/working/raw_strain)",
    re.IGNORECASE,
)
RX_IO = re.compile(
    r"(h5py|read_csv|loadtxt|np\.load|load_residuals|open\s*\([^)]*\.h5)",
    re.IGNORECASE,
)
RX_KERNEL = re.compile(r"\b(kernel|k\(d\)|alpha_geo|beta_tors|omega|phi)\b", re.IGNORECASE)
RX_KERNEL_CONTEXT = re.compile(
    r"(def\s+K\s*\(\s*d\s*\)|K\(d\)|coupling kernel|kernel K\(d\))",
    re.IGNORECASE,
)
RX_FREEZE = re.compile(
    r"(frozen parameter set|kernel frozen|no fitting|bez fittingu|frozen)",
    re.IGNORECASE,
)

RX_ASSIGNMENTS = {
    "alpha_geo": re.compile(r"\balpha_geo\s*=\s*([^\n#]+)"),
    "beta_tors": re.compile(r"\bbeta_tors\s*=\s*([^\n#]+)"),
    "omega": re.compile(r"\bomega\s*=\s*([^\n#]+)"),
    "phi": re.compile(r"\bphi\s*=\s*([^\n#]+)"),
}


def parse_qw_ids(path: Path) -> Set[int]:
    ids: Set[int] = set()
    for m in RX_QW.finditer(path.name):
        try:
            ids.add(int(m.group(1)))
        except ValueError:
            pass
    return {x for x in ids if QW_MIN <= x <= QW_MAX}


def is_probably_canonical(param: str, expr: str) -> bool:
    raw = expr.strip()
    # Symbolic indirection (e.g. ALPHA_GEO = ... elsewhere) is not direct evidence of retuning.
    if re.fullmatch(r"[A-Z_][A-Z0-9_]*", raw):
        return True

    s = expr.strip().lower().replace(" ", "")

    if param == "alpha_geo":
        return (
            "4*np.log(2)" in s
            or "4*numpy.log(2)" in s
            or "2.7726" in s
            or "2.77" in s
        )
    if param == "beta_tors":
        return s.startswith("0.01") or s == "1e-2"
    if param == "omega":
        return "pi/4" in s or "np.pi/4" in s or "0.7854" in s
    if param == "phi":
        return "pi/6" in s or "np.pi/6" in s or "0.5236" in s
    return True


@dataclass
class Row:
    path: str
    qws: List[int]
    has_inference: bool
    has_simulation: bool
    has_external: bool
    has_external_data: bool
    has_kernel: bool
    has_freeze_claim: bool
    noncanonical_assignments: List[Tuple[str, str]]
    category: str


def classify(text: str, qws: List[int], path: Path) -> Row:
    ext = path.suffix.lower()
    is_code = ext == ".py"
    has_inference = bool(RX_INFERENCE.search(text))
    has_simulation = bool(RX_SIMULATION.search(text))
    has_external = bool(RX_EXTERNAL.search(text))
    has_external_strict = bool(RX_EXTERNAL_STRICT.search(text))
    has_io = bool(RX_IO.search(text))
    has_external_data = has_external_strict and has_io
    has_kernel = bool(RX_KERNEL.search(text))
    has_kernel_ctx = bool(RX_KERNEL_CONTEXT.search(text))
    has_freeze = bool(RX_FREEZE.search(text))

    noncanonical: List[Tuple[str, str]] = []
    for name, rx in RX_ASSIGNMENTS.items():
        for m in rx.finditer(text):
            expr = m.group(1).strip()
            if not is_probably_canonical(name, expr):
                noncanonical.append((name, expr[:120]))

    core_mentions = sum(int(tok in text.lower()) for tok in ["alpha_geo", "beta_tors", "omega", "phi"])
    has_kernel_core_context = has_kernel_ctx and core_mentions >= 3
    assigned_params = {name for name, _ in noncanonical}
    has_multicore_noncanonical = len(assigned_params) >= 2

    if is_code and has_inference and has_simulation and not has_external_data:
        category = "INFERENCE_SIMULATION"
    elif is_code and has_inference and has_external_data:
        category = "INFERENCE_EXTERNAL"
    elif is_code and has_external_data and not has_inference:
        category = "EXTERNAL_ANALYSIS_NO_INFERENCE"
    elif is_code and noncanonical and has_kernel_core_context and has_multicore_noncanonical:
        category = "KERNEL_PARAM_SWEEP_CANDIDATE"
    else:
        category = "OTHER"

    return Row(
        path=str(path.relative_to(ROOT)),
        qws=sorted(qws),
        has_inference=has_inference,
        has_simulation=has_simulation,
        has_external=has_external,
        has_external_data=has_external_data,
        has_kernel=has_kernel,
        has_freeze_claim=has_freeze,
        noncanonical_assignments=noncanonical,
        category=category,
    )


def min_qw(rows: List[Row], predicate) -> Optional[int]:
    vals: List[int] = []
    for r in rows:
        if predicate(r):
            vals.extend(r.qws)
    return min(vals) if vals else None


def head_paths(rows: List[Row], n: int = 20) -> List[str]:
    return [r.path for r in rows[:n]]


def main() -> None:
    rows: List[Row] = []
    scanned = 0

    for p in ROOT.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in EXTENSIONS:
            continue

        qws = parse_qw_ids(p)
        if not qws:
            continue

        scanned += 1
        try:
            txt = p.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue

        rows.append(classify(txt, sorted(qws), p))

    inf_rows = [r for r in rows if r.has_inference]
    inf_sim_rows = [r for r in rows if r.category == "INFERENCE_SIMULATION"]
    inf_ext_rows = [r for r in rows if r.category == "INFERENCE_EXTERNAL"]
    ext_rows = [r for r in rows if r.has_external]
    kvar_rows = [r for r in rows if r.category == "KERNEL_PARAM_SWEEP_CANDIDATE"]
    kvar_external_rows = [r for r in kvar_rows if r.has_external_data]
    freeze_rows = [r for r in rows if r.has_freeze_claim]

    out: Dict[str, object] = {
        "scope": {
            "qw_min": QW_MIN,
            "qw_max": QW_MAX,
            "extensions": sorted(EXTENSIONS),
            "files_scanned": scanned,
        },
        "counts": {
            "rows_total": len(rows),
            "rows_with_inference": len(inf_rows),
            "rows_inference_simulation": len(inf_sim_rows),
            "rows_inference_external": len(inf_ext_rows),
            "rows_with_external_markers": len(ext_rows),
            "rows_with_freeze_claim": len(freeze_rows),
            "rows_kernel_param_sweep_candidates": len(kvar_rows),
            "rows_kernel_param_sweep_candidates_with_external_data": len(kvar_external_rows),
        },
        "first_occurrence_qw": {
            "inference_any": min_qw(rows, lambda r: r.has_inference),
            "inference_simulation": min_qw(rows, lambda r: r.category == "INFERENCE_SIMULATION"),
            "inference_external": min_qw(rows, lambda r: r.category == "INFERENCE_EXTERNAL"),
            "external_any": min_qw(rows, lambda r: r.has_external),
            "kernel_param_sweep_candidate": min_qw(
                rows, lambda r: r.category == "KERNEL_PARAM_SWEEP_CANDIDATE"
            ),
        },
        "evidence_heads": {
            "inference_simulation": head_paths(inf_sim_rows),
            "inference_external": head_paths(inf_ext_rows),
            "external_no_inference": head_paths(
                [r for r in rows if r.category == "EXTERNAL_ANALYSIS_NO_INFERENCE"]
            ),
            "kernel_param_sweep_candidates": head_paths(kvar_rows),
            "freeze_claims": head_paths(freeze_rows),
        },
        "caveat": (
            "Static lexical audit. Detects markers, not full causal proof. "
            "Use as screening/audit evidence, not as final semantic adjudication."
        ),
    }

    inf_any = out["counts"]["rows_with_inference"] > 0  # type: ignore[index]
    kretune_external_any = (
        out["counts"]["rows_kernel_param_sweep_candidates_with_external_data"] > 0  # type: ignore[index]
    )

    verdict = {
        "analysis_parameter_tuning_pre1700": (
            "DETECTED" if inf_any else "NOT_DETECTED"
        ),
        "kernel_core_retuning_pre1700_static_signal": (
            "CANDIDATES_ON_EXTERNAL_DATA_DETECTED"
            if kretune_external_any
            else "NO_EXTERNAL_DATA_RETUNING_SIGNAL"
        ),
        "overall": (
            "PRE1700_HAS_INFERENCE_BUT_NO_EXTERNAL_KERNEL_RETUNING_SIGNAL"
            if inf_any and not kretune_external_any
            else "PRE1700_HAS_INFERENCE_AND_EXTERNAL_KERNEL_RETUNING_CANDIDATES"
            if inf_any and kretune_external_any
            else "NO_PRE1700_TUNING_SIGNAL_DETECTED"
        ),
    }
    out["verdict"] = verdict

    OUT_JSON.write_text(json.dumps(out, indent=2, ensure_ascii=True), encoding="utf-8")

    md: List[str] = [
        "# RAPORT QW-1907: PRE-1700 TUNING BOUNDARY AUDIT (QW-700..1699)",
        "",
        "## Wynik",
        f"- Analysis-parameter tuning pre-1700: **{verdict['analysis_parameter_tuning_pre1700']}**",
        f"- Kernel-core retuning pre-1700 (static signal): **{verdict['kernel_core_retuning_pre1700_static_signal']}**",
        f"- Overall: **{verdict['overall']}**",
        "",
        "## Zakres i statystyki",
        f"- files_scanned: {out['scope']['files_scanned']}",
        f"- rows_total: {out['counts']['rows_total']}",
        f"- rows_with_inference: {out['counts']['rows_with_inference']}",
        f"- rows_inference_simulation: {out['counts']['rows_inference_simulation']}",
        f"- rows_inference_external: {out['counts']['rows_inference_external']}",
        f"- rows_with_external_markers: {out['counts']['rows_with_external_markers']}",
        f"- rows_with_freeze_claim: {out['counts']['rows_with_freeze_claim']}",
        f"- rows_kernel_param_sweep_candidates: {out['counts']['rows_kernel_param_sweep_candidates']}",
        (
            "- rows_kernel_param_sweep_candidates_with_external_data: "
            f"{out['counts']['rows_kernel_param_sweep_candidates_with_external_data']}"
        ),
        "",
        "## Pierwsze wystapienia (QW id)",
        f"- inference_any: {out['first_occurrence_qw']['inference_any']}",
        f"- inference_simulation: {out['first_occurrence_qw']['inference_simulation']}",
        f"- inference_external: {out['first_occurrence_qw']['inference_external']}",
        f"- external_any: {out['first_occurrence_qw']['external_any']}",
        f"- kernel_param_sweep_candidate: {out['first_occurrence_qw']['kernel_param_sweep_candidate']}",
        "",
        "## Przykladowe pliki (head)",
        "- inference_simulation:",
    ]
    md.extend([f"  - {x}" for x in out["evidence_heads"]["inference_simulation"][:10]])  # type: ignore[index]
    md.append("- inference_external:")
    md.extend([f"  - {x}" for x in out["evidence_heads"]["inference_external"][:10]])  # type: ignore[index]
    md.append("- kernel_param_sweep_candidates:")
    md.extend([f"  - {x}" for x in out["evidence_heads"]["kernel_param_sweep_candidates"][:10]])  # type: ignore[index]
    md.extend(
        [
            "",
            "## Ograniczenie",
            f"- {out['caveat']}",
            "",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"[OK] wrote: {OUT_JSON.name}")
    print(f"[OK] wrote: {OUT_MD.name}")
    print(f"[VERDICT] {verdict['overall']}")


if __name__ == "__main__":
    main()
