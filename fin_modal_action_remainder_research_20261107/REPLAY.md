# Replay

From the repository root run with Python and SymPy.

```sh
python fin_modal_action_remainder_research_20261107/check_actions.py
```

The recorded SymPy version is stored in verification.json. The script
uses exact symbolic Poisson brackets and exact rational scaling controls.
It checks the cubic generator, action corrections and radiative work,
then expands the quartic density and its canonical correction to classify
all modal-dependent continuum monomials and their localized coefficients.
The pointwise density algebra commutes with spatial integration because
all Poisson derivatives act only on the four modal coordinates.

The remainder estimates use the analytical inequalities in the proof.
The finite scaling controls are norm histories chosen to test the powers.
INPUTS.json pins predecessors and MANIFEST.json hashes package files.
