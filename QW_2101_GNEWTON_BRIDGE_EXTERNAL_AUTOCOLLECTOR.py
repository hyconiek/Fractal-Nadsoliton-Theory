#!/usr/bin/env python3
"""
QW-2101: G_newton bridge external autocollector for QW-2092.

Builds a strict-input JSON with full source metadata and integrity hash.
This script does not claim external independence by itself; it freezes
and validates a provided external bridge snapshot.
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
DEFAULT_SOURCE = ROOT / "external_gnewton_bridge_qw2101.json"
DEFAULT_OUT_INPUT = ROOT / "gnewton_si_bridge_input_qw2092.json"
OUT_JSON = ROOT / "report_qw2101_gnewton_bridge_external_autocollector.json"
OUT_MD = ROOT / "RAPORT_QW2101_GNEWTON_BRIDGE_EXTERNAL_AUTOCOLLECTOR.md"


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


def g_si_from_dimless(g_dimless: float, mu_ref_gev: float, hbar_si: float, c_si: float, gev_to_j: float) -> float:
    g_nat_gev_m2 = g_dimless / (mu_ref_gev * mu_ref_gev)
    return float(g_nat_gev_m2 * hbar_si * (c_si**5) / (gev_to_j**2))


def g_dimless_from_g_si(g_si: float, mu_ref_gev: float, hbar_si: float, c_si: float, gev_to_j: float) -> float:
    g_nat_gev_m2 = g_si * (gev_to_j**2) / (hbar_si * (c_si**5))
    return float(g_nat_gev_m2 * (mu_ref_gev * mu_ref_gev))


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2101 G_newton bridge external autocollector")
    p.add_argument("--source-file", default=str(DEFAULT_SOURCE), help="Source snapshot JSON path.")
    p.add_argument("--output-input", default=str(DEFAULT_OUT_INPUT), help="Output input JSON for QW-2092.")
    p.add_argument("--source", default="codata_gnewton_bridge_snapshot", help="Short source label.")
    p.add_argument("--citation", required=True, help="Human-readable citation.")
    p.add_argument("--reference-url", required=True, help="Reference URL or DOI resolver.")
    p.add_argument("--source-version", required=True, help="Version/release tag for source snapshot.")

    p.add_argument("--mu-ref-gev", type=float, default=None, help="Reference scale in GeV.")
    p.add_argument("--g-si", type=float, default=None, help="Measured G in SI units.")
    p.add_argument(
        "--g-dimensionless-mu-ref",
        type=float,
        default=None,
        help="Dimensionless bridge observable at mu_ref.",
    )
    p.add_argument("--hbar-si", type=float, default=None, help="Planck constant reduced in SI.")
    p.add_argument("--c-si", type=float, default=None, help="Speed of light in SI.")
    p.add_argument("--gev-to-j", type=float, default=None, help="GeV -> Joule conversion.")
    p.add_argument(
        "--strict-dimensionless-only",
        action="store_true",
        help="Reject bridge construction unless a direct external dimensionless observable is provided.",
    )
    p.add_argument(
        "--omit-g-si-optional",
        action="store_true",
        help="Force g_si_input_optional=null in generated QW-2092 input to avoid SI-primary ambiguity.",
    )
    p.add_argument(
        "--require-strict-ready",
        action="store_true",
        help="Return non-zero exit code unless strict provenance conditions are met.",
    )
    args = p.parse_args()

    src_path = Path(args.source_file).resolve()
    out_input = Path(args.output_input).resolve()
    if not src_path.exists():
        raise FileNotFoundError(f"Source file not found: {src_path}")
    src = load_json(src_path)
    src_sha = sha256_file(src_path)

    mu_ref_gev = float(pick(args.mu_ref_gev, src, "mu_ref_gev", 1.0))
    hbar_si = float(pick(args.hbar_si, src, "hbar_si", 1.054571817e-34))
    c_si = float(pick(args.c_si, src, "c_si", 299792458.0))
    gev_to_j = float(pick(args.gev_to_j, src, "gev_to_j", 1.602176634e-10))
    g_si = pick(args.g_si, src, "g_si", None)
    g_dimless = pick(args.g_dimensionless_mu_ref, src, "g_dimensionless_mu_ref", None)
    g_dimless_direct = g_dimless is not None

    if not all(is_finite_positive(x) for x in [mu_ref_gev, hbar_si, c_si, gev_to_j]):
        raise ValueError("mu_ref_gev/hbar_si/c_si/gev_to_j must be finite and > 0.")
    if g_si is None and g_dimless is None:
        raise ValueError("Provide at least one of g_si or g_dimensionless_mu_ref.")

    if g_dimless is None:
        if args.strict_dimensionless_only:
            raise ValueError(
                "strict-dimensionless-only requires direct g_dimensionless_mu_ref (no backsolve from g_si)."
            )
        if not is_finite_positive(g_si):
            raise ValueError("g_si must be finite and > 0.")
        g_dimless = g_dimless_from_g_si(float(g_si), mu_ref_gev, hbar_si, c_si, gev_to_j)
        bridge_observable_origin = "backsolved_from_g_si"
    elif not is_finite_positive(g_dimless):
        raise ValueError("g_dimensionless_mu_ref must be finite and > 0.")
    else:
        src_origin = str(src.get("bridge_observable_origin", "")).strip()
        bridge_observable_origin = src_origin or "external_dimensionless_observable"
        if args.strict_dimensionless_only and bridge_observable_origin != "external_dimensionless_observable":
            raise ValueError(
                "strict-dimensionless-only requires bridge_observable_origin=external_dimensionless_observable."
            )

    g_dimless = float(g_dimless)
    g_si_from_dim = g_si_from_dimless(g_dimless, mu_ref_gev, hbar_si, c_si, gev_to_j)

    rel_delta_g_si = None
    if g_si is not None and is_finite_positive(g_si):
        g_si = float(g_si)
        rel_delta_g_si = abs(g_si_from_dim - g_si) / abs(g_si) * 100.0

    seeded_from_registry = bool(src.get("seeded_from_registry", False))
    provenance_anchor_free = bool(src.get("provenance_anchor_free", True))
    if bridge_observable_origin == "backsolved_from_g_si":
        # Backsolving dimensionless bridge from G_SI is tautological for QW-2092 strict claim.
        provenance_anchor_free = False

    g_si_optional_payload = None if args.omit_g_si_optional else (float(g_si) if g_si is not None else None)
    strict_provenance_ready = bool(
        bridge_observable_origin == "external_dimensionless_observable"
        and provenance_anchor_free
        and not seeded_from_registry
        and g_si_optional_payload is None
    )

    payload = {
        "source": str(args.source),
        "citation": str(args.citation),
        "reference_url": str(args.reference_url),
        "source_sha256": src_sha,
        "source_version": str(args.source_version),
        "provenance_anchor_free": provenance_anchor_free,
        "seeded_from_registry": seeded_from_registry,
        "bridge_observable_origin": bridge_observable_origin,
        "mu_ref_gev": float(mu_ref_gev),
        "g_dimensionless_mu_ref": float(g_dimless),
        "g_si_input_optional": g_si_optional_payload,
        "hbar_si": float(hbar_si),
        "c_si": float(c_si),
        "gev_to_j": float(gev_to_j),
    }
    out_input.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    flags = {
        "bridge_positive": bool(is_finite_positive(g_dimless)),
        "mu_ref_positive": bool(is_finite_positive(mu_ref_gev)),
        "constants_positive": bool(all(is_finite_positive(x) for x in [hbar_si, c_si, gev_to_j])),
        "source_hash_present": bool(bool(src_sha)),
        "self_consistency_if_g_si_provided": bool(rel_delta_g_si is None or rel_delta_g_si <= 1e-6),
        "dimensionless_directly_provided": bool(g_dimless_direct),
        "backsolve_from_g_si_detected": bool(bridge_observable_origin == "backsolved_from_g_si"),
        "strict_provenance_ready": strict_provenance_ready,
    }
    if strict_provenance_ready:
        verdict = "GNEWTON_BRIDGE_EXTERNAL_AUTOCOLLECTED_STRICT_READY"
    elif bridge_observable_origin == "backsolved_from_g_si":
        verdict = "GNEWTON_BRIDGE_EXTERNAL_AUTOCOLLECTED_BACKSOLVED_NONSTRICT"
    else:
        verdict = "GNEWTON_BRIDGE_EXTERNAL_AUTOCOLLECTED_NONSTRICT"

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": verdict,
        "source_file": str(src_path),
        "source_file_sha256": src_sha,
        "output_input_json": str(out_input),
        "bridge": {
            "mu_ref_gev": float(mu_ref_gev),
            "g_dimensionless_mu_ref": float(g_dimless),
            "g_si_from_dimensionless": float(g_si_from_dim),
            "g_si_input_optional": g_si_optional_payload,
            "rel_delta_g_si_pct": rel_delta_g_si,
            "bridge_observable_origin": bridge_observable_origin,
            "strict_provenance_ready": strict_provenance_ready,
        },
        "flags": flags,
        "pass_count": int(sum(1 for v in flags.values() if bool(v))),
        "total_flags": int(len(flags)),
        "required_next_step": (
            f"RUN_QW2103_AND_QW2092_WITH:{out_input.name}"
            if strict_provenance_ready
            else "PROVIDE_EXTERNAL_DIMENSIONLESS_BRIDGE_WITH_ANCHOR_FREE_PROVENANCE_AND_RERUN_QW2101"
        ),
    }
    OUT_JSON.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2101: GNEWTON BRIDGE EXTERNAL AUTOCOLLECTOR",
        "",
        f"- Date UTC: {report['generated_utc']}",
        f"- Verdict: **{report['verdict']}**",
        f"- source_file: `{src_path}`",
        f"- source_file_sha256: `{src_sha}`",
        f"- pass_count: `{report['pass_count']}/{report['total_flags']}`",
        "",
        "## Bridge",
        f"- mu_ref_gev: `{mu_ref_gev:.12g}`",
        f"- g_dimensionless_mu_ref: `{g_dimless:.12e}`",
        f"- g_si_from_dimensionless: `{g_si_from_dim:.12e}`",
        f"- bridge_observable_origin: `{bridge_observable_origin}`",
        f"- strict_provenance_ready: `{strict_provenance_ready}`",
        "",
        "## Output",
        f"- input JSON: `{out_input}`",
        f"- report JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2101] Saved input JSON: {out_input}")
    print(f"[QW-2101] Saved report JSON: {OUT_JSON.name}")
    print(f"[QW-2101] Saved report MD:   {OUT_MD.name}")
    print(
        f"[QW-2101] verdict={report['verdict']} pass_count={report['pass_count']}/{report['total_flags']}"
    )
    if args.require_strict_ready and not strict_provenance_ready:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
