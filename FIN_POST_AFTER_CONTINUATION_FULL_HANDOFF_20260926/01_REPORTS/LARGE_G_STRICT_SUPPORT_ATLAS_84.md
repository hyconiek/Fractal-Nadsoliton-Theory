# LARGE-G-STRICT-SUPPORT-ATLAS-84
## Exhaustive D12 quotient of isolated strict zero-temperature supports

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exhaustive enumeration of all 4095 nonempty label subsets for the supplied midpoint strict spectrum;
- every accepted isolated support has large positivity/gap margins relative to the accepted spectral interval widths;
- singular/degenerate supports are separated into report 87.

For a candidate label support S, solve

    A_SS p_S = mu 1,
    1^T p_S = 1.

Accept it as a strict isolated zero-temperature support when

    p_i > 0 for i in S

and

    (A p)_j < mu for every j outside S.

Modulo D12 there are exactly

    76

isolated strict support classes, representing

    1357

labelled supports.

The orbit counts by support size are:

    size 1: 1 D12 class  / 12 labelled supports
    size 2: 6 D12 classes / 66 labelled supports
    size 3: 11 D12 classes / 208 labelled supports
    size 4: 22 D12 classes / 381 labelled supports
    size 5: 19 D12 classes / 408 labelled supports
    size 6: 14 D12 classes / 234 labelled supports
    size 7: 3 D12 classes / 48 labelled supports.

No isolated strict support of size 8 or larger exists.

Every one of the 76 isolated classes is affinely independent:

    affine_rank(X_S)=|S|-1.

Therefore its limiting covariance rank is exactly |S|-1.

The smallest accepted weight over the full atlas is

    0.032017505652

and the smallest outside-field gap is

    0.020214243847.

The worst augmented equal-field matrix has smallest singular value about

    0.22653,

so the classification is numerically far from a singular decision boundary.

The complete machine-readable atlas is supplied in CSV and JSON.
