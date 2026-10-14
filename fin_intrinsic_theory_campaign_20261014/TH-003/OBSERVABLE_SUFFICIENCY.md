# TH-003 — observable sufficiency and the minimal information lift

## 1. Mean closure is not path closure

The deterministic retained mean `mu=X^T p` closes because `A p=X mu`, so

`dot mu=X^T softmax(g X mu)-mu`.

The tilted retained Hamiltonian is instead

`H_mu(p,eta)=H(p,X eta)`

and contains the hidden factor `sum_i p_i exp[-X_i.eta]`. It is generally not a function of mu alone.

A finite exact witness uses the g=0 control and N=120. Let

`h_j=cos(2*pi*2*j/12)` (the hidden k=2 cosine),

`n_+=10+2h`, `n_-=10-2h`.

The vector `2h=(2,1,-1,-2,-1,1,2,1,-1,-2,-1,1)` is integer-valued, so these are genuine count states. Both have `mu=0` because k2 is orthogonal to the retained k3,k4,k5,k6 feature space.

For the retained c5 feature, `x_c5(j)^2=(Lambda5/12)[1+cos(4*pi*j/12)]`. Direct generator evaluation at g=0 therefore gives

`L(mu_c5^2)(p_+)-L(mu_c5^2)(p_-)=Lambda5/(60*120)=Lambda5/7200`.

With the supplied kernel `Lambda5≈2.298606272079097`, this is `0.0003192508711220968`. Thus equal retained means do not imply equal fluctuation generators.

## 2. An exact 11-dimensional linear lift of the occupation simplex

Write raw Fourier moments

`C_k=sum_j p_j cos(2*pi*k*j/12)`, `S_k=sum_j p_j sin(2*pi*k*j/12)`.

The seven retained feature means determine `C_3,S_3,C_4,S_4,C_5,S_5,C_6` after dividing by their known feature normalizations.

Four entries of `S(p)=X^T diag(p) X` determine four raw products:

- `P34c=E[c3 c4]=(C1+C5)/2`,
- `P34s=E[s3 c4]=(-S1-S5)/2`,
- `P35c=E[c3 c5]=(C2+C4)/2`,
- `P35s=E[s3 c5]=(-S2-S4)/2`.

Hence

`C1=2P34c-C5`,
`S1=-2P34s-S5`,
`C2=2P35c-C4`,
`S2=-2P35s-S4`.

All eleven nonconstant real Fourier coordinates are then known. Fourier inversion recovers the full occupation vector:

`p_j=(1/12)[1+2 sum_{k=1}^5 (C_k cos(k theta_j)+S_k sin(k theta_j))+C_6(-1)^j]`.

The numerical linear map consisting of normalization, seven retained means and these four products has rank 12 with smallest singular value about 0.8601.

### Minimality convention

On the continuous interior of the 11-dimensional simplex, no continuous injective coordinate chart can use fewer than 11 real coordinates. Therefore the **seven retained means plus four selected second moments** are dimension-minimal for exact continuous occupation recovery in this linear/Fourier convention. This is not a finite-type coding theorem: for a fixed finite N, an arbitrary real scalar can encode a finite set if one allows nonphysical discontinuous lookup encodings.

The lift recovers p, not the identities of labelled copies. It is therefore an exact predictive-state lift of the exchangeable occupation process, not eleven new material degrees of freedom.
