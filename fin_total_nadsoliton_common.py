#!/usr/bin/env python3
"""Shared exact/numerical helpers for FIN ST1002--ST1091.

This module does not assume new physics.  It supplies finite-dimensional
state-space examples used to distinguish persistence of the total state from
decay of an internal pattern.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from fin_st402_st416_research import independent_strict_matrix_float


ROOT = Path(__file__).resolve().parent
STRICT_A = independent_strict_matrix_float()
N = int(STRICT_A.shape[0])
ONES = np.ones(N)
UNIFORM = ONES / N


def write_packet(program: int, name: str, status: str, boundary: str,
                 payload: dict[str, Any]) -> dict[str, Any]:
    path = ROOT / f"FIN_ST{program}_{name}.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return {
        "program": f"ST{program}",
        "object": name,
        "packet_file": path.name,
        "packet_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        **payload,
        "status": status,
        "boundary": boundary,
    }


def write_round(lo: int, hi: int, results: dict[str, Any]) -> None:
    assert hi - lo + 1 == 15
    assert list(results) == [f"ST{k}" for k in range(lo, hi + 1)]
    out = ROOT / f"FIN_ST{lo}_ST{hi}_Results.json"
    summary = ROOT / f"FIN_ST{lo}_ST{hi}_Summary.csv"
    out.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    with summary.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["program", "status", "object", "boundary"])
        for key, value in results.items():
            writer.writerow([key, value["status"], value["object"], value["boundary"]])


def spectral_facts() -> dict[str, float]:
    eig = np.linalg.eigvalsh(STRICT_A)
    return {
        "dimension": N,
        "symmetry_residual_inf": float(np.linalg.norm(STRICT_A - STRICT_A.T, ord=np.inf)),
        "row_sum_residual_inf": float(np.linalg.norm(STRICT_A @ ONES, ord=np.inf)),
        "lambda_min": float(eig[0]),
        "spectral_gap": float(eig[1]),
        "lambda_max": float(eig[-1]),
    }


def normalized(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm == 0.0:
        raise ValueError("the zero vector has no normalized representative")
    return vector / norm


def simplex_mass(vector: np.ndarray) -> float:
    return float(np.sum(vector))

