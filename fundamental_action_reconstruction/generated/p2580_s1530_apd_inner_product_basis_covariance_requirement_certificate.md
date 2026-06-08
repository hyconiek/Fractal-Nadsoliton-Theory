# P2580/S1530 APD inner-product basis covariance requirement certificate

Status: `P2580_APD_INNER_PRODUCT_BASIS_COVARIANCE_REQUIREMENT_CERTIFICATE_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Targets audited: `3`.
- Basis variants audited: `4`.
- Naive Euclidean metric fails basis covariance: `True`.
- Transported metric restores basis covariance: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

P2580 separates a real covariance requirement from a source theorem.  If the metric tensor is transported as `G_new = T^T G_ref T`, the same quadratic form gives the same APD solution in every audited basis.  If Euclidean norm is reset in each basis, the P2578 basis artifact reappears.  Covariance fixes coordinates, not the missing strict origin of the metric.

## Recommended next honest step

Use P2580 only as a covariance requirement, not as a source theorem. A valid APD source must provide a metric tensor/inner product whose components transform covariantly under basis changes; merely resetting Euclidean norm in each coordinate chart reproduces P2578, while choosing a transported metric still leaves the underlying metric unsourced as in P2579. The next honest step is to derive the APD inner product from a strict kinetic or measure term.

## Negative controls

No APD coordinate-covariant inner-product source, Gram-metric source, basis-covariance-law source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`7481961a4751d24d43ad9a67d33037bc8cc99c9dedf42b99edb9ea6e5ab1e31d`
