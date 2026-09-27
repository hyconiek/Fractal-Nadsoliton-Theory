# RECORD-IDENTITY-SOURCE-269
## Minimum content rewrite uniquely selects SWAP among all exact two-record reset dilations

Date: 2026-09-27

Status:
general finite-alphabet theorem;
exhaustively replayed for q=2 and q=3.

Let S be an alphabet of size q>=2.

Consider any bijection

    F:S x S -> S x S

used as a reversible dilation of a full one-record reset.

The reset requirement is:

for every fixed system record x, if the environment record e is uniform on S,
then the FIRST output coordinate is uniform on S.

Equivalently, as e ranges over S, the first output realizes every symbol exactly once.

## 1. Content-preservation theorem

Suppose an elementary transport gate is also required to preserve the unordered multiset of the two record contents:

    boxed:
    {output_1,output_2}
      =
    {x,e}

for every input pair.

Then:

    boxed:
    F(x,e)=(e,x)

for every x,e.

### Proof

Fix x.

For an input (x,e), multiset preservation allows the first output to be only:
- x; or
- e.

But the reset condition requires the first output to run through ALL q symbols as e runs through S.

For every e != x, the symbol e can therefore only be produced by choosing first output=e.

For e=x, the first output is x.

Thus first output=e for all inputs.

Multiset preservation then forces second output=x.

Hence F is SWAP.

## 2. Zero-rewrite variational form

Define the content-rewrite cost of one input/output pair as:

    minimum Hamming edits after optimally matching the two output records
    to the two input records.

This cost is nonnegative.

SWAP has zero cost for every pair.

The theorem above implies:

    boxed:
    SWAP is the UNIQUE zero-cost exact reset dilation.

This does NOT require:
- exchange symmetry as a separate assumption;
- diagonal Z3 covariance;
- involutivity as a separate axiom.

They become consequences/properties of the selected zero-rewrite gate.

## 3. Exhaustive checks

For q=2:

    exact reset bijections:
      16

    zero-rewrite exact reset bijections:
      1.

For q=3:

    exact reset bijections:
      46,656

    zero-rewrite exact reset bijections:
      1.

The survivor is SWAP.

## 4. Relation to information continuity

Abstract Shannon-information preservation is weaker.

Every exact reset bijection preserves global state information, including gates that rewrite record content.

The stronger candidate law is:

    boxed:
    elementary transport may move persistent records but should not rewrite their intrinsic content.

If FIN can derive that notion of persistent record content, SWAP and all additive composition conservation follow automatically.

## 5. Boundary

"Minimize record-content rewrite" is now a sharp source candidate, not yet a theorem of FIN.

The remaining foundational question is:

    what FIN object carries record identity/content strongly enough that local recoding F_± counts as rewriting rather than a gauge-equivalent transport?

That is the correct continuation of the information-continuity lane.
