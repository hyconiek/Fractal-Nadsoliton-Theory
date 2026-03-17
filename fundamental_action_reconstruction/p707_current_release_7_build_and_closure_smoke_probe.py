#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from shutil import copyfile
from shutil import which
from typing import Any


AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_TEX = REPO / "TOE_FINAL_DOCUMENTATION_RELEASE_7_STRICT_FULL.tex"
OUT_PDF = Path("/tmp/TOE_FINAL_DOCUMENTATION_RELEASE_7_STRICT_FULL.pdf")
OUT_PDF_REPO_ROOT = REPO / "TOE_FINAL_DOCUMENTATION_RELEASE_7_STRICT_FULL.pdf"

P706 = ROOT / "p706_current_release_7_strict_projective_operational_toe_os_closure_dashboard_probe.py"
P706_SUMMARY = GENERATED / "p706_current_release_7_strict_projective_operational_toe_os_closure_dashboard_probe_summary.json"

P441 = ROOT / "p441_current_strict_global_closure_next_move_dashboard_probe.py"
P441_SUMMARY = GENERATED / "p441_current_strict_global_closure_next_move_dashboard_probe_summary.json"

OUT_JSON = GENERATED / "p707_current_release_7_build_and_closure_smoke_probe.json"
OUT_SUMMARY = GENERATED / "p707_current_release_7_build_and_closure_smoke_probe_summary.json"


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def run_cmd(cmd: list[str], cwd: Path) -> dict[str, Any]:
    completed = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    stdout_tail = (completed.stdout or "")[-2000:]
    stderr_tail = (completed.stderr or "")[-2000:]
    return {
        "cmd": cmd,
        "returncode": int(completed.returncode),
        "stdout_tail": stdout_tail,
        "stderr_tail": stderr_tail,
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    if not IN_TEX.exists():
        missing.append(str(IN_TEX.relative_to(REPO)))
    if not P706.exists():
        missing.append(str(P706.relative_to(REPO)))
    if not P441.exists():
        missing.append(str(P441.relative_to(REPO)))

    if missing:
        artifact = {
            "stage": "P707",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    if which("pdflatex") is None:
        artifact = {
            "stage": "P707",
            "status": "NOT_COMPUTABLE_MISSING_PDFLATEX",
            "as_of": AS_OF,
            "error": "pdflatex not found in PATH",
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    tex_build_1 = run_cmd(
        [
            "pdflatex",
            "-interaction=nonstopmode",
            "-halt-on-error",
            "-output-directory=/tmp",
            str(IN_TEX),
        ],
        cwd=REPO,
    )
    tex_build_2 = run_cmd(
        [
            "pdflatex",
            "-interaction=nonstopmode",
            "-halt-on-error",
            "-output-directory=/tmp",
            str(IN_TEX),
        ],
        cwd=REPO,
    )

    pdf_ok = (tex_build_1["returncode"] == 0) and (tex_build_2["returncode"] == 0) and OUT_PDF.exists()
    repo_pdf_ok = False
    repo_pdf_error = None
    if pdf_ok:
        try:
            copyfile(OUT_PDF, OUT_PDF_REPO_ROOT)
            repo_pdf_ok = OUT_PDF_REPO_ROOT.exists()
        except Exception as exc:  # pragma: no cover
            repo_pdf_error = str(exc)

    p706_run = run_cmd([sys.executable, str(P706)], cwd=REPO)
    p706_ok = p706_run["returncode"] == 0 and P706_SUMMARY.exists()
    p706_status = None
    if p706_ok:
        p706_status = (load_json(P706_SUMMARY) or {}).get("status")

    p441_run = run_cmd([sys.executable, str(P441)], cwd=REPO)
    p441_ok = p441_run["returncode"] == 0 and P441_SUMMARY.exists()
    p441_recommended = None
    if p441_ok:
        p441_recommended = (load_json(P441_SUMMARY) or {}).get("recommended_next_strict_target")

    checks = [
        {"id": "pdf_build", "pass": bool(pdf_ok), "pdf_path": str(OUT_PDF)},
        {
            "id": "pdf_copied_to_repo_root",
            "pass": bool(repo_pdf_ok) if pdf_ok else True,
            "pdf_path": str(OUT_PDF_REPO_ROOT),
            "error": repo_pdf_error,
        },
        {"id": "p706_dashboard", "pass": bool(p706_ok), "status": p706_status},
        {"id": "p441_dashboard", "pass": bool(p441_ok), "recommended_next_strict_target": p441_recommended},
    ]

    all_ok = all(bool(c.get("pass")) for c in checks)
    status = "PASS_RELEASE_7_BUILD_AND_CLOSURE_SMOKE_READY" if all_ok else "REQUIRES_REVIEW_RELEASE_7_BUILD_AND_CLOSURE_SMOKE"

    artifact = {
        "stage": "P707",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "operational_release_7_build_and_strict_dashboard_smoke_only",
        "inputs": {
            "tex": str(IN_TEX.relative_to(REPO)),
            "P706": str(P706.relative_to(REPO)),
            "P441": str(P441.relative_to(REPO)),
        },
        "outputs": {
            "pdf": str(OUT_PDF),
            "P706_summary": str(P706_SUMMARY.relative_to(REPO)),
            "P441_summary": str(P441_SUMMARY.relative_to(REPO)),
        },
        "runs": {
            "pdflatex_pass1": tex_build_1,
            "pdflatex_pass2": tex_build_2,
            "P706": p706_run,
            "P441": p441_run,
        },
        "checks": checks,
        "hard_limits": [
            "Operational smoke probe only; does not claim any physics identification.",
            "No kernel-alone/global QW-2191 discharge.",
            "No directed/sign-sensitive physical orientation datum promotion.",
            "No Standard Model host matching claim.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P707",
        "status": status,
        "as_of": AS_OF,
        "pdf_build_ok": bool(pdf_ok),
        "pdf_path": str(OUT_PDF),
        "p706_status": p706_status,
        "recommended_next_strict_target": p441_recommended,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
