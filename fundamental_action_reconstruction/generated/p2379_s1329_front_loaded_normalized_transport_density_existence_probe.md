# P2379 S1329: front-loaded normalized transport density existence probe

Status: `OPEN_PROGRESS_FRONT_LOADED_NORMALIZED_TRANSPORT_PROFILE_EXISTS_SOURCE_PROFILE_STILL_OPEN`

## Result

P2379/S1329 refines the P2378 no-go: unit-uniform endpoint transport fails, but a normalized front-loaded transport profile can pass the d5 chamber audit.
For `rho_lambda(s)=1+lambda*(1-2s)`, `int rho_lambda ds=1` and `int rho_lambda*A_s(d) ds=C(d)+lambda*B(d)`.

## Certificate

- Tested front-loading: `lambda=0.8` with density range `[0.19999999999999996, 1.8]`.
- 81x81 lattice max threshold row: `{'eta': 1.8, 'beta_tors': 0.1, 'K1': 0.4699856726450201, 'K5': 0.02413122336363006, 'C1': 0.5978370007556204, 'C5': 2.545243208781978, 'B1': 0.059216335936026265, 'B5': 0.9786115376296642, 'lambda_numerator': 0.6340935880563139, 'lambda_denominator_B5_minus_3B1': 0.8009625298215854, 'lambda_threshold_gt': 0.7916644842269442}`.
- 81x81 lattice min threshold row: `{'eta': 2.0, 'beta_tors': 0.0, 'K1': 0.4699856726450201, 'K5': 0.02413122336363006, 'C1': 0.6931471805599453, 'C5': 3.258096538021482, 'B1': 0.07944154167983575, 'B5': 1.5187442610632012, 'lambda_numerator': 0.20717079822978413, 'lambda_denominator_B5_minus_3B1': 1.280419636023694, 'lambda_threshold_gt': 0.16179914178225746}`.
- Lambda test exceeds lattice max: `True`.
- Max midpoint integral error: `1.2117462055982742e-07`.
- Uniform grid scans select: `[{'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}]`.
- Front-loaded grid scans select: `[{'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}]`.

## Hard limits

- This is a normalized-profile existence/probe result, not a strict variational source theorem deriving `rho_lambda`.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
