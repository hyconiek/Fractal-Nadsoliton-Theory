# P2380 S1330: front-loaded profile rectangle monotonicity certificate

Status: `OPEN_PROGRESS_FRONT_LOADED_PROFILE_RECTANGLE_MONOTONICITY_CERTIFIED_SOURCE_STILL_OPEN`

## Result

P2380/S1330 upgrades the P2379 front-loaded profile from grid/lattice evidence to a rectangle monotonicity certificate for the affine normalized density family.
Interval arithmetic certifies `dT/deta<0` and `dT/dbeta_tors>0` for `T=[3*(K1+C1)-(K5+C5)]/(B5-3*B1)` on the P2376 rectangle.

## Certificate

- Interval boxes: `16384`.
- Worst threshold box: `{'eta_interval': {'lo': 1.8, 'hi': 1.8015625}, 'beta_tors_interval': {'lo': 0.09921875000000001, 'hi': 0.1}, 'threshold_interval': {'lo': 0.7731224845561374, 'hi': 0.8040446708364509}, 'threshold_eta_derivative_interval': {'lo': -3.0826748914481965, 'hi': -2.8244163139794605}, 'threshold_beta_tors_derivative_interval': {'lo': 2.3982665499837355, 'hi': 2.6648216874939807}, 'numerator_interval': {'lo': 0.6291026492071081, 'hi': 0.6362250268671104}, 'denominator_interval': {'lo': 0.791280696139983, 'hi': 0.8137166642725266}}`.
- Worst `dT/deta` upper bound: `-1.3301876312346335`.
- Worst `dT/dbeta_tors` lower bound: `1.1577604521646625`.
- Corner maximum candidate: `{'eta': 1.8, 'beta_tors': 0.1, 'K1': 0.4699856726450201, 'K5': 0.02413122336363006, 'C1': 0.5978370007556204, 'C5': 2.545243208781978, 'B1': 0.05921633593602582, 'B5': 0.9786115376296642, 'lambda_numerator': 0.6340935880563139, 'lambda_denominator_B5_minus_3B1': 0.8009625298215868, 'lambda_threshold_gt': 0.7916644842269429}`.
- Lambda test rectangle-uniform by monotonicity: `True`.
- Sample support replay selections: `[{'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}]`.

## Hard limits

- This certifies the affine profile threshold on the P2376 rectangle; it does not derive the profile from strict dynamics.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
