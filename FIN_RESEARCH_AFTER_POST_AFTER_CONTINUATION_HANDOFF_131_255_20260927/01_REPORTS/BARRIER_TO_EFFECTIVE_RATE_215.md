# BARRIER-TO-EFFECTIVE-RATE-215
## The q3/q4/q5 hierarchy follows the d3/d4/d5 saddle ordering, but finite-N exponents are not yet asymptotic barrier exponents

Date: 2026-09-26

Status:
finite-N numerical comparison, not an Eyring-Kramers theorem.

At the working gain the mapped localized-state index-one saddles have barriers

    Delta V3 ≈ 0.644851906971
    Delta V4 ≈ 0.662219456770
    Delta V5 ≈ 0.782655029416.

Under CRT:

    d3 = pure Z4/fiber move,
    d4 = pure Z3/base move,
    d5 = mixed move.

The effective rates satisfy throughout N=3..8:

    boxed:
    q3 > q4 > q5.

Ratios:


    N=3:
      q3/q4=1.078459
      exp[N(Delta4-Delta3)]=1.053484

      q4/q5=1.796859
      exp[N(Delta5-Delta4)]=1.435204

    N=4:
      q3/q4=1.083229
      exp[N(Delta4-Delta3)]=1.071940

      q4/q5=1.882142
      exp[N(Delta5-Delta4)]=1.618893

    N=5:
      q3/q4=1.090671
      exp[N(Delta4-Delta3)]=1.090720

      q4/q5=1.979686
      exp[N(Delta5-Delta4)]=1.826091

    N=6:
      q3/q4=1.098995
      exp[N(Delta4-Delta3)]=1.109828

      q4/q5=2.092625
      exp[N(Delta5-Delta4)]=2.059809

    N=7:
      q3/q4=1.108443
      exp[N(Delta4-Delta3)]=1.129272

      q4/q5=2.224818
      exp[N(Delta5-Delta4)]=2.323440

    N=8:
      q3/q4=1.118821
      exp[N(Delta4-Delta3)]=1.149056

      q4/q5=2.376029
      exp[N(Delta5-Delta4)]=2.620813

The agreement is already surprisingly close for some intermediate N, especially
around N=5-6, but the effective finite-N slopes are smaller than the raw saddle
barrier differences.

Descriptive N=3..8 fits give approximately:

    q3 ~ exp(-0.55046 N),
    q4 ~ exp(-0.55790 N),
    q5 ~ exp(-0.61373 N).

So:

    beta4-beta3 ≈ 0.00744

versus the saddle difference

    Delta4-Delta3 ≈ 0.01737,

and

    beta5-beta4 ≈ 0.05584

versus

    Delta5-Delta4 ≈ 0.12044.

Conclusion:

- the ordering is consistent with the saddle landscape;
- N<=8 is not in a regime where one may identify the fitted exponents with the
  asymptotic barriers;
- prefactors, competing paths and finite-N recrossings remain material.

Do not promote the current ratios to an Eyring-Kramers certificate.
