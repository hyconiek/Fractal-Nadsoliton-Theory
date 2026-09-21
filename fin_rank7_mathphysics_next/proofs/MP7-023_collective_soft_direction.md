# MP7-023 — collective soft direction and scope

Scientific state: **PROVED_SCOPED ANALYTIC CONSEQUENCE**.

On the aligned nonnegative C4 family, use the accepted entrywise-nonnegative covariance
matrix `M4`. Where `M4` is irreducible, Perron--Frobenius supplies a simple leading
eigenvalue with a strictly positive eigenvector. If the entries are strictly positive,
irreducibility is automatic.

The C4 Hessian is `H4=I/g-M4`, so its softest Cartesian direction is the leading
covariance direction. At a simple fold the null vector therefore coincides with that
leading eigendirection (up to sign), because `M4 v=(1/g)v`. MP7-025's positive fold
vector is consistent with this structure.

Target P constrains only the remaining three covariance eigenvalues:
`lambda2(M4)<=67/250`. Thus for `g<250/67`, three Hessian directions have the explicit
lower curvature bound `1/g-67/250>0`; the sign of the leading direction is not fixed by
Target P. A smooth leading spectral projector additionally requires a positive
separation between the first and second covariance eigenvalues.

This is a collective mode in the supplied amplitude coordinates. It does not select a
physical label, orientation, species or time direction, and it does not turn covariance
into a spacetime metric.

## Quantitative isolation at the certified fold

MP7-020 strengthens the global Target-P ceiling to

`lambda2(M4) <= 0.2679999463710577`.

At the certified fold, `H4 v=0` and `||v||=1`, so exactly

`M4 v=(1/g_fold)v`.

The R7P-031 gain box gives

`1/g_fold in [0.2844428483815984, 0.2844428499997532]`.

Therefore the eigenvalue carried by the fold vector is separated from every other
C4 covariance eigenvalue by at least

`0.2844428483815984 - 0.2679999463710577 > 0.0164429020105407`.

All four components of the certified null-vector box are strictly positive.  Hence
this eigenpair is the unique Perron eigenpair, and the leading spectral projector is
smooth in a neighborhood in which this gap persists.  The fold vector is therefore
not merely consistent with a collective Perron mode: it is rigorously identified
with the unique leading covariance direction in the certified fold box.

For orientation only, its cosine with the positive diagonal direction `(1,1,1,1)`
is interval-enclosed in `[0.9926559574,0.9926560369]`; this angle is descriptive and
is not used in the proof.
