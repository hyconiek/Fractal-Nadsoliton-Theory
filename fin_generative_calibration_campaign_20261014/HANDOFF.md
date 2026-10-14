# FIN full handoff — N2 calibration validity and falsification, RG-F01–RG-F05

Date 2026-10-14. **Status: NO_NEW_PHYSICAL_INPUT / CONDITIONAL_TEST_DESIGN**. RG-F06 is architect-only and was NOT executed.

## Results

- F01: reused and checked exact positive difference quotient, nine combinations of three matrices × three kinetics. Independent 144-state replay checks pass; 13 outcomes incl INVALID compose for arbitrarily long finite sequences under a fixed difference-based nondisturbing emission.
- F02: exact model symmetry gate, single h=.01 phase-breaking witness 0.002711545661 in transition-rate defect; invalid-phase dependence and a uniform-valid erasure detector are explicit counterexamples.
- F03: inherited **conditional/hypothetical** gain and detector margin **0.00119114**; disjoint added errors must cost less than this before statistical power. No real calibration bounds and no all-operation dynamic interval cert.
- F04: complete finite-data iid confidence theorem for 13 outcomes. A conservative example with eps=rho/3 and validation alpha=.025 requires **40,278,676 independent records** if that rho is physically validated; not measured or automatically feasible. The continuous-gain composite rejection example distinguishes synthetic ideal FIN/flat and correctly declines erasure data.
- F05: `NO_NEW_PHYSICAL_INPUT`. The mathematical/statistical workflow is a design, NOT physical validation.

## Evidence

Independent verifier: **98/98 PASS**; max independent 144-vs-12 joint-law coordinate discrepancy 1.110e-16. Binary64 reproduction is not a source/physics certificate.

## Source provenance

The original prior RG-E package, as mounted from this conversation, is included under `source_context` unmodified; `source_model.py` is its copied local mathematical reimplementation rather than versioned original `safe_seed.py` bytes. The new uploaded plan is copied byte-for-byte. The RG-E intake supplied by user in conversational text was used as a reviewed claim summary, **not** an independently downloaded file or reverified repo SHA. Repository baseline `08b6d235` is user-reported, not verified against live Git in this run. No predecessor review GAIN_BOX_CERTIFICATE or numerical calibration records are present: its reported sharp bound is inherited explicitly as a conditional premise.

## Replay order

`python run_calibration.py`, `python decision_rule.py`, `python independent_verify.py`, `python write_handoff.py`, `python finalize_package.py`; `python validate_records.py path/to/raw.csv` is syntax-only and does not certify data or iid. See REPLAY.md.

## Disposition

Return RG-F01–F05 artifacts, full probability arrays, proofs, raw-record schema, independent checks, source manifests and no-input certificate for the architect's RG-F06 review. **Do not create an experiment, spend money, move to a new campaign or promote A7 as a physical source.**
