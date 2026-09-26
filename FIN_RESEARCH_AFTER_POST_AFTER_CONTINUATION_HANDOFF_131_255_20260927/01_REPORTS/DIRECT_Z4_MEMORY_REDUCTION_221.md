# DIRECT-Z4-MEMORY-REDUCTION-221
## Microscopic -> Z4 reduction is as accurate as the previously emphasized Z3 reduction

Date: 2026-09-26

Status:
exact N=6 finite-state Mori-Zwanzig calculation.

Project the microscopic chain directly onto:
- four localized sectors j mod4;
- one residual class.

Masses are approximately:

    0.248124
    0.248124
    0.248124
    0.248124
    residual 0.007504.

The nontrivial Z4 modes are:


    beta=1  (Z12 k=3):
      A=-0.109541167119
      M0=0.0876406910529
      M1=0.0266917468907
      M0/|A|=0.800071
      lambda_MZ=-0.021331111439
      lambda_exact=-0.0213246614497
      relative error=0.030247 %

    beta=2  (Z12 k=6):
      A=-0.131588560029
      M0=0.10446544613
      M1=0.0323240996846
      M0/|A|=0.793879
      lambda_MZ=-0.0262738358109
      lambda_exact=-0.0262621436215
      relative error=0.044521 %

    beta=3  (Z12 k=9):
      A=-0.109541167119
      M0=0.0876406910529
      M1=0.0266917468907
      M0/|A|=0.800071
      lambda_MZ=-0.021331111439
      lambda_exact=-0.0213246614497
      relative error=0.030247 %

Thus direct microscopic->Z4 first-moment closure has errors only about
0.03-0.045%.

This is essentially the same quality as the Z3 closure.

Representation theory explains why:
the k=3,6,9 sectors are symmetry-orthogonal to the other resolved basin
characters, so enlarging the intermediate resolved space does not change their
sector-wise scalar elimination.

Conclusion:

    Z3 is not uniquely selected merely by "best MZ closure".

Both CRT quotients are controlled effective observables.
