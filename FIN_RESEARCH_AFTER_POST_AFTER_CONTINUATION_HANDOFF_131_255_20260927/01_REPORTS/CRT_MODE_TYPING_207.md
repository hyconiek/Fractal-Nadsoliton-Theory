# CRT-MODE-TYPING-207
## The hidden H4 modes are mixed Z3×Z4 characters, while retained X7 contains pure base and pure fiber modes

Date: 2026-09-26

Status:
exact representation-theoretic retyping of the existing Fourier modes.

Using

    j=4a+9b mod12,

the Z12 Fourier character

    exp(2 pi i k j/12)

becomes

    exp[
      2 pi i(
        alpha a/3
        +beta b/4
      )
    ]

with

    boxed:
    alpha=k mod3,
    beta=3k mod4.

Therefore every Fourier mode has a definite Z3×Z4 character.

Key modes:

    k=3  -> (0,1)  fiber-only
    k=4  -> (1,0)  base-only
    k=5  -> (2,3)  mixed
    k=6  -> (0,2)  fiber-only

The hidden modes are:

    k=1 -> (1,3) mixed
    k=2 -> (2,2) mixed

plus their conjugates.

## Consequence

The hidden H4 sector contains NO pure-base and NO pure-fiber character.

Every hidden mode transforms nontrivially under both factors.

By contrast the retained X7 sector contains:
- a pure Z3 base pair k=4;
- pure Z4 fiber modes k=3 and k=6;
- one mixed pair k=5.

## Memory interpretation boundary

This makes it mathematically plausible that eliminating H4 can generate
base/fiber cross-memory, because H4 consists of mixed representation sectors.

But representation mixing is not yet a theorem that all observed MZ memory
comes specifically from H4.

The finite-N coarse memory also includes:
- intrabasin degrees;
- transition-state structure;
- occupation fluctuations.

So the next test must decompose the actual memory kernel by CRT character
sector rather than infer its source from labels alone.
