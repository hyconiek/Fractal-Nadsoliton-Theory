# Metropolis regularity gate: PASSED in the narrow local theorem sense

The Metropolis acceptance f(d)=min(1,e^d) is globally Lipschitz on the bounded field range and is not C1 at d=0. In the FIN fold, some equal-logit pairs do occur; no unsupported C3 statement is made.

Detailed balance supplies the **exact interior gradient form** b_MET=-K_MET(p,g) grad F_g. The edge logarithmic mean applied to positive, locally Lipschitz intensities is locally Lipschitz. Since grad F_g(p_f) is constant on coordinates and K1=0, multiplication by K(p)-K_f annihilates the constant force and leaves O((|x|+||w||+|delta_g|)(||w||+x²+|delta_g|)). This pays the needed second-order soft drift and fast remainder without differentiating the kink. Finite-N leave-one-out field errors are O(N^-1) by Lipschitz acceptance. No global smoothness is asserted.

This regularity argument would fail for a discontinuous acceptance rule or a cusp whose local modulus decays slower than linearly; those are outside the frozen controls. It also does not license applying a theorem requiring globally C3 microscopic rates by name.
