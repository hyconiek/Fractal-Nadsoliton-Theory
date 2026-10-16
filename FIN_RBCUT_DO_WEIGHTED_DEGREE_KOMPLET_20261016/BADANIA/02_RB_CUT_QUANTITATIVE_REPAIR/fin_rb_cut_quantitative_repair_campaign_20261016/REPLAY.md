# Independent replay

From the extracted new campaign directory:

```
PYTHONDONTWRITEBYTECODE=1 python3 verify_repair.py
PYTHONDONTWRITEBYTECODE=1 python3 verify_package.py
```

`verify_repair.py` uses only Python stdlib: `fractions.Fraction`, exhaustive cut checks on 2/4/8-state artificial networks, and an independent Edmonds–Karp rational max-flow. It does not import previous producer code. Its reported flags check source balance, cut-optimal K, exact primal construction, failed smaller caps, the need for signed corrections, and the false converse to uniform native-relative domination.

`verify_package.py` hashes every package member listed in `MANIFEST.sha256`, checks manifest coverage and bound locally accessible input files in `INPUTS.sha256`. Old inaccessible review materials are not silently replaced. Neither checker proves the all-N actual MT graph claim; the new max-flow theorem is analytic.
