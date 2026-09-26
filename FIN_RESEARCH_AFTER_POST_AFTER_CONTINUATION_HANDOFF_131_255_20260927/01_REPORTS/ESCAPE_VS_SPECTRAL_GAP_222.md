# ESCAPE-VS-SPECTRAL-GAP-222
## Z3 is more metastable by escape/capacity, while Z4 is slower by global quotient relaxation

Date: 2026-09-26

Status:
exact effective-chain identities plus microscopic capacity replay at N=6,8.

The apparent Z3/Z4 conflict comes from comparing two different notions of
"slow".

## 1. Escape / conductance

For one Z3 sector:

    Phi3
      =
      total exit rate
      =
      2(q1+q2+q4+q5).

For one Z4 sector, define

    a=q1+q3+q5
      nearest quotient jump,

    b=2q2+q6
      opposite quotient jump.

Then

    Phi4
      =
      2a+b.

For every N=3..8:

    boxed:
    Phi3 < Phi4.

Results:


    N=3:
      Phi3=0.0880039711107
      Phi4=0.10006781429
      Phi3/Phi4=0.879443

    N=4:
      Phi3=0.0479029220196
      Phi4=0.0545175499938
      Phi3/Phi4=0.878670

    N=5:
      Phi3=0.026816504438
      Phi4=0.0305757123496
      Phi3/Phi4=0.877052

    N=6:
      Phi3=0.0150792747696
      Phi4=0.0172340146722
      Phi3/Phi4=0.874972

    N=7:
      Phi3=0.00842794070859
      Phi4=0.00966145623823
      Phi3/Phi4=0.872326

    N=8:
      Phi3=0.00465994801779
      Phi4=0.00536169520803
      Phi3/Phi4=0.869118

So a trajectory remains longer, on average, inside one mod3 sector.

## 2. Global relaxation gap

Nevertheless:

    gap_Z4 < gap_Z3.

A four-state cyclic quotient can possess a smaller long-wavelength eigenvalue
even though each individual sector leaks faster.

Thus:

    residence/escape time
      !=
    global mixing time.

## 3. Raw microscopic capacity agrees with the escape ranking

At N=6:

    raw cap/pi mod3 ≈ 0.0626081
    raw cap/pi mod4 ≈ 0.0723712.

At N=8:

    raw cap/pi mod3 ≈ 0.0224401
    raw cap/pi mod4 ≈ 0.0266736.

So both the raw microscopic whole-basin capacity and the memory-renormalized
effective conductance rank mod3 sectors as more isolated.

## 4. Interpretation

If the intended object is a metastable "unit" that should preserve identity
before escaping, escape/capacity is the more directly relevant criterion.

If the intended object is the slowest global collective relaxation coordinate,
the Z4 quotient currently wins slightly.

These are different scientific questions and should not be collapsed into one
"best coarse variable".
