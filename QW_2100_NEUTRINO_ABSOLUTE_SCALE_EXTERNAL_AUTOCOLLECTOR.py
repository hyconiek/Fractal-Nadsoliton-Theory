#!/usr/bin/env python3
"""
QW-2100: external neutrino absolute-scale autocollector for QW-2091.

Builds a strict-input JSON with full source metadata and integrity hash.
This script does not claim external independence by itself; it only freezes
and validates the provided external snapshot inputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict


ROOT = Path(__file__).resolve().parent
DEFAULT_SOURCE = ROOT / "external_neutrino_absolute_scale_qw2100.json"
DEFAULT_OUT_INPUT = ROOT / "neutrino_absolute_scale_input_qw2091.json"
OUT_JSON = ROOT / "report_qw2100_neutrino_absolute_scale_external_autocollector.json"
OUT_MD = ROOT / "RAPORT_QW2100_NEUTRINO_ABSOLUTE_SCALE_EXTERNAL_AUTOCOLLECTOR.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def pick(cli: Any, src: Dict[str, Any], key: str, default: Any = None) -> Any:
    if cli is not None:
        return cli
    if key in src and src[key] is not None:
        return src[key]
    return default


def is_finite_positive(x: Any) -> bool:
    try:
        v = float(x)
    except Exception:  # noqa: BLE001
        return False
    return math.isfinite(v) and v > 0.0


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2100 neutrino absolute-scale external autocollector")
    p.add_argument("--source-file", default=str(DEFAULT_SOURCE), help="Source snapshot JSON path.")
    p.add_argument("--output-input", default=str(DEFAULT_OUT_INPUT), help="Output input JSON for QW-2091.")
    p.add_argument("--source", default="planck_bao_neutrino_sum_snapshot", help="Short source label.")
    p.add_argument("--citation", required=True, help="Human-readable citation.")
    p.add_argument("--reference-url", required=True, help="Reference URL or DOI resolver.")
    p.add_argument("--source-version", required=True, help="Version/release tag for source snapshot.")

    p.add_argument("--mode", choices=["sum_mnu", "m_beta"], default=None, help="Absolute-scale observable mode.")
    p.add_argument("--sum-mnu-ev", type=float, default=None, help="Absolute sum neutrino mass in eV.")
    p.add_argument("--m-beta-eff-ev", type=float, default=None, help="Effective beta-decay mass in eV.")
    p.add_argument("--sigma-ev", type=float, default=None, help="One-sigma uncertainty for selected observable.")
    p.add_argument("--sum-mnu-upper-bound-ev", type=float, default=None, help="Optional upper bound for sum(mnu).")
    p.add_argument("--hierarchy", default=None, help="Hierarchy mode (currently QW-2091 supports normal).")
    p.add_argument("--dm21-ev2", type=float, default=None, help="Delta m21^2 in eV^2.")
    p.add_argument("--dm31-ev2", type=float, default=None, help="Delta m31^2 in eV^2.")
    args = p.parse_args()

    src_path = Path(args.source_file).resolve()
    out_input = Path(args.output_input).resolve()
    if not src_path.exists():
        raise FileNotFoundError(f"Source file not found: {src_path}")
    src = load_json(src_path)
    src_sha = sha256_file(src_path)

    mode = str(pick(args.mode, src, "mode", "sum_mnu")).strip().lower()
    hierarchy = str(pick(args.hierarchy, src, "hierarchy", "normal")).strip().lower()
    sigma_ev = float(pick(args.sigma_ev, src, "sigma_ev", 0.01))
    dm21_ev2 = float(pick(args.dm21_ev2, src, "dm21_ev2", 7.42e-5))
    dm31_ev2 = float(pick(args.dm31_ev2, src, "dm31_ev2", 2.517e-3))
    sum_mnu_upper_bound_ev = pick(args.sum_mnu_upper_bound_ev, src, "sum_mnu_upper_bound_ev", None)

    if not is_finite_positive(sigma_ev):
        raise ValueError("sigma_ev must be finite and > 0.")
    if not is_finite_positive(dm21_ev2) or not is_finite_positive(dm31_ev2):
        raise ValueError("dm21_ev2 and dm31_ev2 must be finite and > 0.")
    if mode not in {"sum_mnu", "m_beta"}:
        raise ValueError(f"Unsupported mode: {mode}")

    payload: Dict[str, Any] = {
        "source": str(args.source),
        "citation": str(args.citation),
        "reference_url": str(args.reference_url),
        "source_sha256": src_sha,
        "source_version": str(args.source_version),
        "provenance_anchor_free": True,
        "hierarchy": hierarchy,
        "mode": mode,
        "sigma_ev": float(sigma_ev),
        "dm21_ev2": float(dm21_ev2),
        "dm31_ev2": float(dm31_ev2),
    }

    if mode == "sum_mnu":
        sum_mnu_ev = pick(args.sum_mnu_ev, src, "sum_mnu_ev", None)
        if not is_finite_positive(sum_mnu_ev):
            raise ValueError("mode=sum_mnu requires finite positive sum_mnu_ev.")
        payload["sum_mnu_ev"] = float(sum_mnu_ev)
    else:
        m_beta_eff_ev = pick(args.m_beta_eff_ev, src, "m_beta_eff_ev", None)
        if not is_finite_positive(m_beta_eff_ev):
            raise ValueError("mode=m_beta requires finite positive m_beta_eff_ev.")
        payload["m_beta_eff_ev"] = float(m_beta_eff_ev)

    if sum_mnu_upper_bound_ev is not None:
        if not is_finite_positive(sum_mnu_upper_bound_ev):
            raise ValueError("sum_mnu_upper_bound_ev must be finite and > 0 when provided.")
        payload["sum_mnu_upper_bound_ev"] = float(sum_mnu_upper_bound_ev)

    out_input.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    flags = {
        "mode_supported": bool(mode in {"sum_mnu", "m_beta"}),
        "hierarchy_present": bool(hierarchy),
        "observable_present": bool(("sum_mnu_ev" in payload) or ("m_beta_eff_ev" in payload)),
        "sigma_positive": bool(is_finite_positive(payload["sigma_ev"])),
        "mass_splittings_positive": bool(is_finite_positive(payload["dm21_ev2"]) and is_finite_positive(payload["dm31_ev2"])),
        "source_hash_present": bool(bool(src_sha)),
    }

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": "NEUTRINO_ABSOLUTE_SCALE_EXTERNAL_AUTOCOLLECTED",
        "source_file": str(src_path),
        "source_file_sha256": src_sha,
        "output_input_json": str(out_input),
        "mode": mode,
        "hierarchy": hierarchy,
        "payload_preview": {
            "sum_mnu_ev": payload.get("sum_mnu_ev"),
            "m_beta_eff_ev": payload.get("m_beta_eff_ev"),
            "sigma_ev": payload.get("sigma_ev"),
            "sum_mnu_upper_bound_ev": payload.get("sum_mnu_upper_bound_ev"),
        },
        "flags": flags,
        "pass_count": int(sum(1 for v in flags.values() if bool(v))),
        "total_flags": int(len(flags)),
        "required_next_step": f"RUN_QW2091_WITH:{out_input.name}",
    }
    OUT_JSON.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2100: NEUTRINO ABSOLUTE SCALE EXTERNAL AUTOCOLLECTOR",
        "",
        f"- Date UTC: {report['generated_utc']}",
        f"- Verdict: **{report['verdict']}**",
        f"- source_file: `{src_path}`",
        f"- source_file_sha256: `{src_sha}`",
        f"- mode: `{mode}`",
        f"- hierarchy: `{hierarchy}`",
        f"- pass_count: `{report['pass_count']}/{report['total_flags']}`",
        "",
        "## Output",
        f"- input JSON: `{out_input}`",
        f"- report JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2100] Saved input JSON: {out_input}")
    print(f"[QW-2100] Saved report JSON: {OUT_JSON.name}")
    print(f"[QW-2100] Saved report MD:   {OUT_MD.name}")
    print(
        f"[QW-2100] verdict={report['verdict']} pass_count={report['pass_count']}/{report['total_flags']}"
    )


if __name__ == "__main__":
    main()

