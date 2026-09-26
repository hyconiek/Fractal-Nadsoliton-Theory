# LOCALIZED-12-STATE-MARKOV-GENERATOR-214
## Memory-aware reduction yields a positive D12-circulant effective generator on all 12 localized basins for N=3..8

Date: 2026-09-26

Status:
- D12-circulant form follows exactly from symmetry;
- positivity of all reconstructed shell rates is numerical for N=3..8;
- N=6 is independently validated against exact full-generator eigenmodes.

Keep the 12 localized basins as resolved states, ordered by their dominant
label j in Z12.

For a reflection-symmetric circulant Markov generator, let q_d be the jump
rate to each basin at cyclic distance d=1,...,5, and q_6 the antipodal rate.

Its Fourier eigenvalues are

    lambda_k
      =
      sum_(d=1)^5
      2 q_d[
        cos(2 pi k d/12)-1
      ]
      +
      q_6[
        cos(pi k)-1
      ].

The M0+M1 Mori-Zwanzig rates were inverted to q_d.

Results:


    N=3:
      q1=0.00179034869
      q2=0.00932054908
      q3=0.0227890128
      q4=0.0211310797
      q5=0.0117600081
      q6=0.00874797691
      all positive=True

    N=4:
      q1=0.00126917356
      q2=0.00479252953
      q3=0.0126549862
      q4=0.0116826528
      q5=0.0062071051
      q6=0.0046699612
      all positive=True

    N=5:
      q1=0.00079921456
      q2=0.00254190843
      q3=0.00729500082
      q4=0.00668854178
      q5=0.00337858745
      q6=0.00254628983
      all positive=True

    N=6:
      q1=0.000468550925
      q2=0.00135413049
      q3=0.00425133082
      q4=0.00386837877
      q5=0.0018485772
      q6=0.00138883579
      all positive=True

    N=7:
      q1=0.000258887038
      q2=0.000715830917
      q3=0.00247712206
      q4=0.00223477624
      q5=0.00100447616
      q6=0.000748823886
      all positive=True

    N=8:
      q1=0.000136873225
      q2=0.000373516559
      q3=0.0014327756
      q4=0.00128061241
      q5=0.00053897181
      q6=0.000397420823
      all positive=True

For every tested N:
- all six shell rates are positive;
- q3 is the largest rate;
- q4 is the second largest.

Thus the memory-reduced operator is not merely symmetric:
it is a bona-fide effective Markov generator on the localized-state ring.

## Exact full-generator validation at N=6

For each Fourier sector k=1,...,6, the corresponding microscopic slow
eigenmode has >95% overlap with the localized-basin Fourier subspace.

The M0+M1 relative eigenvalue errors are all below 0.083%.

The exact N=6 shell rates inferred directly from full microscopic eigenvalues
are:


    q1=0.000470645731792

    q2=0.00135372597139

    q3=0.00424680663344

    q4=0.00386463648571

    q5=0.00184808354015

    q6=0.0013893428767


All are positive.

The MZ shell rates differ from the exact spectral rates by at most about 0.45%
for q1 and about 0.11% or less for the dominant q3/q4 rates.

Therefore the 12-state Markov description is dynamically supported, not only
an algebraic fit.
