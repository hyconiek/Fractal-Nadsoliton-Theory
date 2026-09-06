# FIN: thirty adaptive falsification studies

Completed programme ST8549--ST8578: one bounded investigation per identifier,
with an English analytical report, numerical results, 33 scientific tests,
figures and a replay manifest. No experiment, physical confirmation,
independent audit or claim of mathematical priority is made.

Initial state-map audit: root AGENTS.md (section index, top-level rules,
current frontier and selected older kernel/duality/refinement rules),
SUMMARY_GROK.md, the current compendium, and ST8547/ST8548 proofs and code.
This is not a fresh line-by-line verification of every historical report.

## Principal result

For a declared finite Markov population on X^N, autonomous singleton rates Q
plus zero nonnegative pair-activity budget

    B(x) = sum_y g(x,y) binom(number_of_changed_coordinates(x,y), 2)

force the generator to be the independent tensor sum of Q. If sup B <= eps,
transition distributions differ from that tensor model by at most
min(1, 2 eps t) in total variation. Initial independence is a separate premise.

The common-transposition counterexample shows that the same FIN singleton
paths, product equilibrium, detailed balance, graph symmetry, exchangeability
and projective consistency do NOT force B=0. Whether an actual FIN composition
law supplies this premise is the primary open source question.

Other results concern the positive double cover of signed legacy, the
measurement-limit W-squared obstruction, nonunique Lindblad/coherence
extensions, exact detector aliases, nonlinear source response and hidden
memory. The report also gives a conditional finite-shot test with exact
rational type-I and type-II error bounds. No physical apparatus is inferred.

## Replay

Requires Python 3, NumPy, SciPy, SymPy and Matplotlib; pdflatex is needed only
to rebuild the report. No additional software is installed by these scripts.
From this directory:

```sh
python3 analysis.py --through 30
python3 -m unittest discover -s . -p 'test_*.py' -v
python3 verify.py
python3 build_artifacts.py
pdflatex -interaction=nonstopmode -halt-on-error report.tex
pdflatex -interaction=nonstopmode -halt-on-error report.tex
python3 build_artifacts.py
sha256sum -c SHA256SUMS.txt
python3 package.py --verify
```

`analysis.py --through 5/10/15/20/25/30` replays the corresponding completed
checkpoint. It overwrites results.json with that checkpoint; return to 30
before building the final report. Replay checks a fixed completed campaign;
it does not invent a new adaptive research sequence.

`build_artifacts.py` can run from a portable extracted package without the
original repository. Unavailable historical sources are then recorded as
unavailable rather than being fabricated. The supplied original provenance
manifest preserves the actual source hashes from this campaign.

Analytic proofs, not floating residuals or status labels, license mathematical
statements. Exact rational tests are explicitly identified. The NumPy/SciPy
values are not interval certificates. Existing archives were preserved. The
preceding ST8547 and ST8548 test suites also passed (15 + 18 tests), but those
older packages are not bundled here.

## Files

- `report.pdf`: complete English report, proofs and thirty-step appendix.
- `report.tex`: source of that report.
- `RESEARCH_LEDGER_EN.md`, `results.json`: hypotheses, outcomes, numerical
  evidence and next-step rationale for all thirty investigations.
- `analysis.py`, `test_analysis.py`: computation and scientific checks.
- `verify.py`, `verification.json`: executed test logs and an exact replay
  comparison with the archived result JSON.
- `build_artifacts.py`, `findings.pdf`, `findings.png`: deterministic figure
  and document-support builder.
- `provenance.json`, `SHA256SUMS.txt`: environment and artifact integrity.
- `package.py`: checks all manifest hashes, builds the ZIP and can replay it
  after extraction into a new temporary directory.

Creator: Krzysztof Żuchowski, ORCID 0009-0002-0909-3613.
Affiliation: Independent Researcher — Fractal Information Theory Project.
Report license: CC BY 4.0. No DOI assigned and no upload performed.
