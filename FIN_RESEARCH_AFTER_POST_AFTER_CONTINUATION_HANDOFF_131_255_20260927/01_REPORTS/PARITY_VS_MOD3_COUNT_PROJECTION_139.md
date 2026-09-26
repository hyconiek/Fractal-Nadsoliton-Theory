# PARITY-VS-MOD3-COUNT-PROJECTION-139
## A dynamics-selected three-sector observable retains more information than parity count

Date: 2026-09-26

Status:
exact finite-state calculations for the leave-one-out Gibbs generator,
N=3..6, at

    g=5.145228719489142.

Two equilibrium projections are compared.

### Binary-oriented projection

    Y2 = number of copies on even labels.

Resolved-state count:

    N+1.

### Three-sector projection

    Y3 =
      (number on labels 0 mod 3,
       number on labels 1 mod 3,
       number on labels 2 mod 3).

Resolved-state count:

    (N+1)(N+2)/2.

This projection is motivated by the actual low-barrier partition

    {0,3,6,9},
    {1,4,7,10},
    {2,5,8,11}.

## 1. Strong lumpability

Neither projection is exactly strongly lumpable.

Using the equilibrium-weighted RMS spread of microscopic rates inside each
coarse fiber, the three-sector count projection is systematically better.

    N=3: lump 0.801563->0.743323; memory area 0.960749->0.850885; max semigroup error 0.007857->0.006494
    N=4: lump 0.842473->0.778194; memory area 0.596813->0.589305; max semigroup error 0.016645->0.012350
    N=5: lump 0.861128->0.795720; memory area 0.493084->0.484505; max semigroup error 0.020974->0.016274
    N=6: lump 0.867365->0.804336; memory area 0.429744->0.420467; max semigroup error 0.024929->0.020301

The improvement is modest but uniform over N=3..6.

## 2. Memory

For both projections define the exact Mori-Zwanzig kernel

    K(t)=C^T exp(t QSQ) C.

The three-sector count projection has slightly smaller integrated normalized
memory over t in [0,16] for every tested N.

At N=6:

    parity-count area
      ≈ 0.429743761,

    mod3-count area
      ≈ 0.420467005.

The memory tail also decays slightly faster for mod3 counts.

## 3. Projected-semigroup error

Compare

    B^T exp(tS) B

with the instantaneous Markov closure

    exp[t B^T S B].

At N=6, over t in {0.1,0.5,1,2,4,8}:

    parity-count maximum error
      ≈ 0.024929;

    mod3-count maximum error
      ≈ 0.020301.

So retaining the three residue counts improves the exact-process projection
without fitting any new parameter.

## 4. Interpretation

This result supports the review recommendation that branching number should be
selected by dynamics rather than imposed as binary.

However the three-sector count process is still not exactly Markov.

The next question is therefore not whether memory exists, but whether:
- it decays rapidly;
- its integrated effect can be absorbed into a controlled effective generator.

Reports 141-142 answer that question for the hard metastable-sector projection.
