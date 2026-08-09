#!/usr/bin/env python3
"""FIN P485/P487 exact campaign with Singular as the Groebner engine."""

from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import traceback
import zipfile

import sympy as sp

import fin_phase_exact_algebra as algebra


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(algebra.json_ready(payload), indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def singular_text(expression: str) -> str:
    return expression.replace("**", "^")


def checkpoint(root: Path, label: str) -> Path:
    path = root / f"FIN_Singular_Checkpoint_{label}.zip"
    temporary = path.with_suffix(".zip.tmp")
    patterns = ("FIN_Singular_*.json", "FIN_Singular_*.log", "FIN_Singular_*.sing",
                "FIN_Singular_*.txt", "FIN_Singular_*.txt.gz")
    files = sorted({p for pattern in patterns for p in root.glob(pattern)})
    with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED) as archive:
        for source in files:
            if source not in (path, temporary):
                archive.write(source, source.name)
    temporary.replace(path)
    return path


def prepare(root: Path) -> dict[str, object]:
    status = root / "FIN_Singular_Prepare_Status.json"
    atomic_json(status, {"stage": "constructing_exact_reduced_system", "started": time.time()})
    context = algebra.exact_context()
    reduced = algebra.reduced_branch_system(context)
    tangent = algebra.tangent_consistency(context)
    branch = reduced["branch_substitution"]
    targets = [sp.expand(value.subs(branch)) for value in tangent["consistency_polynomials"]]
    localizer = sp.symbols("pivot_inverse", real=True)
    variables = reduced["variables"][:-1] + (localizer, reduced["variables"][-1])
    generators = list(reduced["equations"]) + [reduced["alpha_polynomial"],
        sp.expand(localizer*reduced["pivot_determinant"]-1)]
    # Mandatory closure audit before any remote CAS work.
    for expression in generators + targets:
        if expression.has(sp.sqrt(2)):
            raise AssertionError(f"unclosed sqrt(2) coefficient: {expression}")
        sp.Poly(expression, *variables, domain=sp.QQ)

    system = algebra.polynomial_system()
    alpha = context["alpha"]
    s1, s2, s3 = system["sines"]
    substitution = {s1:alpha, s2:1-2*alpha**2, s3:3*alpha-4*alpha**3}
    A,B,u,L,a,b,c,d,e,f,g,h,i = system["variables"]
    full_variables = (A,B,u,a,b,c,d,e,f,g,h,i,alpha,L)
    full_generators = [sp.expand(value.subs(substitution)) for value in system["equations"]]
    full_generators.append(8*alpha**4-8*alpha**2+1)
    for expression in full_generators:
        sp.Poly(expression, *full_variables, domain=sp.QQ)

    payload: dict[str, object] = {
        "status": "prepared_exact_QQ_system",
        "coefficient_closure": "sqrt(2)=2-4*alpha^2",
        "reduced_variables": [str(value) for value in variables],
        "reduced_generators": [str(value) for value in generators],
        "p485_targets": [str(value) for value in targets],
        "active_pivot_rows": list(reduced["pivot"]["rows"]),
        "rho_reference_orbit": tangent["reference"],
        "full_variables": [str(value) for value in full_variables],
        "full_generators": [str(value) for value in full_generators],
        "prepared_unix_time": time.time(),
    }
    atomic_json(root / "FIN_Singular_Exact_Input.json", payload)
    atomic_json(status, {"stage": "completed", "finished": time.time(),
                         "generators": len(generators), "targets": len(targets)})
    checkpoint(root, "after_PREPARE")
    return payload


def load_input(root: Path) -> dict[str, object]:
    path = root / "FIN_Singular_Exact_Input.json"
    if not path.exists():
        return prepare(root)
    return json.loads(path.read_text(encoding="utf-8"))


def run_singular(root: Path, program: str, script_text: str) -> tuple[int, str]:
    executable = shutil.which("Singular") or shutil.which("singular")
    if not executable:
        raise RuntimeError("Singular executable not found")
    script = root / f"FIN_Singular_{program}.sing"
    log = root / f"FIN_Singular_{program}.log"
    status = root / f"FIN_Singular_{program}_Status.json"
    script.write_text(script_text, encoding="utf-8")
    atomic_json(status, {"program": program, "stage": "singular_running", "started": time.time(),
                         "executable": executable})
    with log.open("w", encoding="utf-8") as handle:
        process = subprocess.Popen([executable, "-q", str(script)], cwd=root,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            bufsize=1, start_new_session=True)
        output: list[str] = []
        try:
            assert process.stdout is not None
            for line in process.stdout:
                output.append(line)
                handle.write(line); handle.flush()
                print(line, end="", flush=True)
            return_code = process.wait()
        except KeyboardInterrupt:
            os.killpg(process.pid, signal.SIGINT)
            try: process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGTERM); process.wait(timeout=15)
            atomic_json(status, {"program":program, "stage":"stopped_by_user",
                                 "return_code":process.returncode, "finished":time.time()})
            checkpoint(root, f"after_{program}_STOPPED")
            raise
    atomic_json(status, {"program": program, "stage": "completed" if return_code == 0 else "failed",
                         "return_code": return_code, "finished": time.time()})
    checkpoint(root, f"after_{program}")
    return return_code, "".join(output)


def markers(output: str, prefix: str, count: int) -> list[str]:
    values = []
    for index in range(1, count+1):
        begin = f"FIN_{prefix}_{index}_BEGIN"
        end = f"FIN_{prefix}_{index}_END"
        if begin not in output or end not in output:
            values.append("<missing marker>")
        else:
            values.append(output.split(begin,1)[1].split(end,1)[0].strip())
    return values


def p486(root: Path) -> None:
    result = subprocess.run([sys.executable, "-u", str(root/"fin_program_486.py")], cwd=root)
    if result.returncode:
        raise RuntimeError(f"P486 failed with code {result.returncode}")
    checkpoint(root, "after_P486")


def p485(root: Path) -> None:
    data = load_input(root)
    names = data["reduced_variables"]
    generators = [singular_text(value) for value in data["reduced_generators"]]
    targets = [singular_text(value) for value in data["p485_targets"]]
    lines = ["option(redSB);", "option(prot);", f"ring R=0,({','.join(names)}),dp;",
             "ideal I=" + ",\n".join(generators) + ";", "ideal G=slimgb(I);",
             'write(\":w FIN_Singular_P485_Basis.txt\",G);']
    for index, target in enumerate(targets, 1):
        lines += [f'print("FIN_REMAINDER_{index}_BEGIN");', f"print(reduce({target},G));",
                  f'print("FIN_REMAINDER_{index}_END");']
    lines += ["quit;"]
    code, output = run_singular(root, "P485", "\n".join(lines)+"\n")
    remainders = markers(output, "REMAINDER", len(targets))
    p486_data = json.loads((root/"FIN_Program_486_Results.json").read_text())
    passed = code == 0 and all(value == "0" for value in remainders) and p486_data["orientation_premises_paid"]
    atomic_json(root/"FIN_Singular_P485_Results.json", {
        "status": "computer_assisted_proof" if passed else "open",
        "all_five_remainders_zero": all(value == "0" for value in remainders),
        "remainders": remainders, "singular_return_code": code,
        "boundary": "Only the declared localized positive branch is tested."
    })
    checkpoint(root, "after_P485_RESULT")


def p487(root: Path) -> None:
    data = load_input(root)
    names = data["reduced_variables"]
    generators = [singular_text(value) for value in data["reduced_generators"]]
    eliminated = "*".join(name for name in names if name != "L")
    lines = ["option(redSB);", "option(prot);", f"ring R=0,({','.join(names)}),dp;",
             "ideal I=" + ",\n".join(generators) + ";",
             f"ideal E=eliminate(I,{eliminated});", "ideal GE=std(E);",
             'write(\":w FIN_Singular_P487_Elimination.txt\",GE);',
             'print("FIN_ELIMINATION_1_BEGIN");', "print(GE);",
             'print("FIN_ELIMINATION_1_END");', "quit;"]
    code, output = run_singular(root, "P487", "\n".join(lines)+"\n")
    elimination = markers(output, "ELIMINATION", 1)[0]
    has_relation = code == 0 and elimination not in ("0", "", "<missing marker>")
    atomic_json(root/"FIN_Singular_P487_Results.json", {
        "status": "exact_elimination_relation_returned" if has_relation else "open",
        "elimination_ideal": elimination, "singular_return_code": code,
        "minimal_polynomial_proved": False,
        "boundary": "Any relation still requires P473 factor isolation, irreducibility, and exact verification."
    })
    checkpoint(root, "after_P487_RESULT")


def p475(root: Path) -> None:
    data = load_input(root)
    names = data["full_variables"]
    generators = [singular_text(value) for value in data["full_generators"]]
    eliminated = "*".join(name for name in names if name != "L")
    lines = ["option(redSB);", "option(prot);", f"ring R=0,({','.join(names)}),dp;",
             "ideal I=" + ",\n".join(generators) + ";",
             f"ideal E=eliminate(I,{eliminated});", "ideal GE=std(E);",
             'write(\":w FIN_Singular_P475_Elimination.txt\",GE);',
             'print("FIN_ELIMINATION_1_BEGIN");', "print(GE);",
             'print("FIN_ELIMINATION_1_END");', "quit;"]
    code, output = run_singular(root, "P475", "\n".join(lines)+"\n")
    atomic_json(root/"FIN_Singular_P475_Results.json", {
        "status": "completed" if code == 0 else "open", "singular_return_code":code,
        "elimination_ideal":markers(output,"ELIMINATION",1)[0],
        "minimal_polynomial_proved":False})
    checkpoint(root, "after_P475_RESULT")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("prepare","p486","p485","p487","p475","recommended"))
    parser.add_argument("--root", default=".")
    args = parser.parse_args()
    root = Path(args.root).resolve(); root.mkdir(parents=True, exist_ok=True)
    os.chdir(root)
    actions = ("p486","prepare","p485","p487") if args.action == "recommended" else (args.action,)
    state_path = root/"FIN_Singular_Campaign_State.json"
    state = {"status":"running", "actions":actions, "started":time.time(), "completed":[]}
    atomic_json(state_path,state)
    try:
        for action in actions:
            state["current"] = action; atomic_json(state_path,state)
            {"prepare":prepare,"p486":p486,"p485":p485,"p487":p487,"p475":p475}[action](root)
            state["completed"].append(action); atomic_json(state_path,state)
        state.update(status="completed", current=None, finished=time.time())
        atomic_json(state_path,state); checkpoint(root,"FINAL")
        return 0
    except KeyboardInterrupt:
        state.update(status="stopped_by_user", current=None, finished=time.time())
        atomic_json(state_path,state); checkpoint(root,"STOPPED"); return 130
    except BaseException as error:
        state.update(status="failed", current=None, finished=time.time(),
                     error_type=type(error).__name__, error=str(error), traceback=traceback.format_exc())
        atomic_json(state_path,state); checkpoint(root,"FAILED"); raise


if __name__ == "__main__":
    raise SystemExit(main())
