# SRN-001 — Ideal angular event: no tuned frequency or post-hoc threshold

Frozen feature maps f1(j)=sqrt(2)cos(2pi j/12), f3(j)=sqrt(2)cos(6pi j/12). Frozen raw event |W3|>|W1|; ties are assigned zero, no data-dependent block selection. N is a multiple of 12; k=floor(N^(63/64)); M=floor(N^(1/4)), alpha=k/N, eta=N^(-1/32). Fixed lag/spacing = one model-time unit. For sufficiently large N, 1<=k<=N/2 (this is NOT true for every finite N).

At gamma in {gamma0,gamma+}, the OU reference angular probability equals (2/pi) atan sqrt(1+alpha(1/gamma-1)). The difference grows linearly in alpha at zero, with the exact derivative d_ang=(gamma+^(-1)-gamma0^(-1))/(2pi)>0. Under the available q12 spectral inputs nominal d_ang≈0.0264034996200004. Values are nominal unless certified separately. The gap is asymptotically small: this is a many-record test, not a powerful one-shot discriminator.

The raw event is exactly a common-gain invariant, but the ideal threshold used for distinguishing the two models is a *known-model* analytical midpoint. Invariance of the event does not make the chosen models or their probabilities self-calibrating in nature.
