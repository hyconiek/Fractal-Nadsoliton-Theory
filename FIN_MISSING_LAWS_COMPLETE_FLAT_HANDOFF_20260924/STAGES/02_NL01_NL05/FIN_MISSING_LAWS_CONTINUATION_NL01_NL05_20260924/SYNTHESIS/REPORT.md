# NL-01..NL-05 synthesis

The continuation sharpens the proposed architecture but does not close the
physics source problem.

## What strengthened

1. **Geometry/measure:** on the actually sourced S1 phase geometry, the dual
   measure and face conductance are intrinsic. The 1D operator is conservative,
   self-adjoint and consistent even on strongly irregular dense meshes.
2. **State decomposition:** the correct local architecture is one translation
   coordinate alpha plus a **four-dimensional transverse deformation vector r**,
   not alpha plus one beta scalar.
3. **Memory interpretation:** a reversible extended state remains compatible
   with coarse memory, but compatibility is not sourcehood.

## What failed to source

1. Multipole cutoff L=2 remains a free law choice. Current covariance,
   symmetry, composition and locality do not select it.
2. Spatial refinement plus energy conservation do not select an internal
   kinetic tensor, momentum variable or symplectic bracket.
3. The certified quadratic transverse sector propagates but does not localize.

## Updated architecture

`carrier / mediator -> intrinsic S1 geometry -> (cell mass M, stiffness K) -> alpha + r_4 -> [missing kinetic metric + flow state] -> conditional reversible transport`.

This is narrower than the previous architecture: two of the missing-law
candidates have been promoted to exact/conditional mathematical constructions,
while the two source problems have been isolated more sharply.
