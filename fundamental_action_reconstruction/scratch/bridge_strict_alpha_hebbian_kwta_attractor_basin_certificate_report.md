# Scratch strict-alpha Hebbian kWTA attractor-basin certificate probe

Status: conditional d5-teacher Hebbian dynamics certificate with exact tie nonclosure.

- Supports scanned: `792` five-active-node states on `Z_12`.
- Exact first row of centered zero-self Hebbian d5 weights: `['0', '-25/12', '11/12', '-1/12', '-13/12', '23/12', '-25/12', '23/12', '-13/12', '-1/12', '11/12', '-25/12']`.
- Deterministic lexicographic kWTA reaches d5 count: `786/792`.
- Deterministic non-d5 attractor/cycle count: `6`.
- Deterministic max steps to d5 among reached states: `4` with histogram `{'0': 12, '1': 444, '2': 289, '3': 33, '4': 8}`.
- Set-valued all-tie-branches guaranteed basin: `768/792`.
- Set-valued some-branch/existential basin: `792/792`.
- Tie-sensitive blockers: `24/792`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: Hebbian kWTA gives a strong but nonclosed conditional attractor mechanism after a d5 teacher trace is supplied; exact ties still require extra selector data.
- No false pass: no d5-source theorem, no QW-2191 discharge, no ToE closure.
