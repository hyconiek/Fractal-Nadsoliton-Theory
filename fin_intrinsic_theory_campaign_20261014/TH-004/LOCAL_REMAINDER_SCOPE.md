# TH-004 local remainder scope

MP7-026 already pays the static fold Taylor/reduction remainder on its certified tube: Cartesian radius 0.002, `0<=epsilon<=1e-7`, transverse inverse bound, `T3/T4` tensor bounds, and `|r-r0|<=0.005965415221 epsilon` on the final branch bracket.

For TH-005, the stochastic path is not restricted to the stationary branch bracket. Smoothness of the finite log-sum-exp gives local Taylor remainders of schematic order

`O(epsilon*|xi| + |xi|^3 + |xi|*||y|| + ||y||^2 + epsilon||y||)`

in the soft drift. Under `epsilon=N^-2/3`, `xi=N^-1/3 Z`, these terms vanish after multiplication by `N^(2/3)` **provided the stable coordinate is controlled in the time-integrated/perturbed-test sense**. Paying that stochastic fast-mode estimate is the remaining TH-005 convergence lemma; the static MP7 remainder alone does not prove it.
