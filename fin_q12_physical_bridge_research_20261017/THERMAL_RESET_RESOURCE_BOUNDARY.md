# A conditional implementation bridge, not a natural kinetic source

Date: 2026-10-17. Exact comparison constructions; no new FIN physical law
is adopted. Quantum thermal operations here are imported comparison tools,
not a derivation of QM from classical FIN histories.

## 1. Fixed single-ancilla SWAP test

For the cavity r let H_S(r)=diag(epsilon_b(r)) and suppose a q12 ancilla
has one fixed Hamiltonian H_B and thermal state tau_B. A full SWAP sends
the tagged state to tau_B, independently of r. It cannot produce all
the distinct native s(r) without context-dependent preparation.

Moreover

    [SWAP,H_S(r) tensor I + I tensor H_B]=0
       iff H_S(r)-H_B=c(r) I.                       (1)

Proof: conjugating the sum by SWAP exchanges its factors. Their difference
is `D tensor I-I tensor D`, D=H_S-H_B, which vanishes exactly when D is
scalar (compare product-basis matrix entries). The parity perturbation
in the preceding witness changes energy differences between even and odd
labels, not merely a common offset. One fixed matched q12 ancilla cannot
meet (1) for both contexts. This is a restricted implementation obstruction,
not a theorem forbidding larger energy-conserving baths.

## 2. Explicit positive control: a designed finite context bank

At fixed N there are M=binom(N+10,11) possible cavity count vectors.
As an analysis construction, provide one q12 ancilla B_r for every context,
with Hamiltonian H_Br=H_S(r), each independently prepared at inverse beta.
For a tagged update, use the unitary controlled on the background count r
which swaps the tag with B_r. On each r-block the system conditional
Hamiltonian agrees with that ancilla up to a background constant. The
controlled operation commutes with the supplied total energy and resets
the label exactly to s(r).

Fresh independent banks and supplied rate-nu Poisson attempts reproduce
the native tagged update and hence the finite parent when this construction
is supplied for all copies. A straightforward bank for one update has
dimension `12^M`; this is a construction cost, **not a lower bound** or a
practical design. A compact bank or a genuinely autonomous realization
would require a separate theorem.

Everything doing physical work is explicit: energy identification J/beta,
the context-controlled interaction, thermal preparation, fresh independence,
and the arrival clock. Reusing the same finite ancillas does not give fresh
resets or automatically the native Markov path law. This is a designed
realization of an already given target, not an independent origin of HB,
A7 or its clock. No apparatus has been built or licensed.

## 3. Conservation still does not select the activity

Replace the full SWAP by `U_eta=cos(eta) I-i sin(eta) SWAP`, with the same
matched energies. It still commutes with the total energy. For diagonal
classical inputs, tracing out the fresh ancilla gives

    p -> cos²(eta) p+sin²(eta) s(r).

Thus the supplied Poisson stream produces off-diagonal rates
`nu sin²(eta) s_b(r)`. Different pulse choices give different activities
while preserving the same conditional Gibbs law and energy-conserving
collision. The cavity is unchanged in a forward/reverse tagged transition,
so even a fixed chosen function eta(r) preserves the global Gibbs detailed
balance. This remains a diagnostic implementation family, not adopted FIN.

In particular, thermal consistency and a convenient unitary dilation do
not independently choose the complete reset or its attempt clock. Counting
this target-designed bank as an AS/KB provider would be circular.

## 4. Relation to the literature and the next paid gate

Strasberg et al. explicitly supply fresh interaction units and waiting-time
statistics; their Poisson limit gives the generator gamma(J-I), with
energy/entropy/work bookkeeping. That framework can audit this constructed
bank, but it does not source its preparation, controller or clock. Its
conditions are not inferred solely from detailed balance.

The next physical-provider task must supply an evolution/coupling which is
specified independently of the desired reset, account for its resources,
and pass the zero-gap spectator test. A purely mathematical conditional
implementation is a useful bridge to **engineering thermodynamics**. A
fundamental closed-world theory needs an autonomous pre-Markov carrier
and composition, not an unexplained external reservoir below the nadsoliton.
The ontology and all existing physical-source/calibration gates remain intact.
