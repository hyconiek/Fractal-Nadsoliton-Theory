# R13 Partial Host To Canonical Block Overlap Packet

Status: `R13_EXECUTED_PARTIAL_HOST_TO_CANONICAL_BLOCK_OVERLAP_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R12/P19/N22`, one remaining blocker is:

```text
explicit_host_to_submatrix_matching_witness_identifying_the_QW2186_certified_host_operator_with_the_exported_canonical_Psi_block_on_full_declared_transport_support_or_its_declared_control_pullback
```

`R13` does not try to fake that matching witness.

It asks a narrower constructive question:

```text
does the repo already contain a real partial overlap packet between the
QW-2186 host and the exported canonical Psi block, so that the final matching
blocker can be reduced to a smaller finite set?
```

## Inputs reused

1. `QW-2118`
   - frozen `K_total` carrier and cyclic distance profile.
2. `QW-2124`
   - branch-resolved scalar vacuum floor.
3. `QW-2186`
   - certified host `A = K_total + m0^2 I`.
4. `R8`
   - host carrier packet on the declared `12`-slot basis.
5. `R12`
   - explicit canonical `Psi-Psi` block export.

## Result of `R13`

`R13` materializes a real partial overlap:

1. host side:
   - frozen `K_total` off-diagonal kernel on the `12`-slot carrier,
   - broken-branch scalar floor `m0^2 I`,
   - explicit host decomposition `A_host = K_total + m0^2 I`;
2. canonical block side:
   - off-diagonal kernel channel `(K_i_j + K_j_i)/2`,
   - diagonal local channel
     `3*g4_psii*vpsii**2 + 5*g6_psii**4 + 2*gYi*vphi**2 + m2_psii`;
3. light-facing content is explicit on both sides through the kernel channel.

This is a real overlap packet.
It is not yet a full matching witness.

## Honest frontier after `R13`

`R13` does not establish:

- coefficient specialization from the symbolic canonical kernel channel to the
  frozen numeric `K_total` matrix,
- equality between host diagonal floor `m0^2 I` and the canonical local
  diagonal sector,
- full host-to-block matching witness,
- selector-relevant canonicalization inside `QW-2191`.

So the matching blocker is now narrower:

- partial host/block overlap = present,
- kernel specialization witness = absent,
- diagonal equality witness = absent.

## What `R13` does not claim

`R13` does not claim:

- theorem-level PASS,
- full-closure PASS,
- that the host operator already equals the exported canonical block,
- that `QW-2191` is discharged,
- that selector closure is obtained,
- that ToE is closed.

## Recommended next move

The correct next move is now:

1. rerun the host-identification route after this partial overlap packet,
2. accept only:
   - a shorter missing-object list,
   - or the unchanged negative route.
