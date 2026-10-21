# What changed in the actual NM gate

Date: 2026-10-21. Original NM files/manifests and statuses are preserved.

The source's [NEXT_OPEN_LEMMA](../fin_nonlinear_mode_memory_20261021/NEXT_OPEN_LEMMA.md)
remains a planning snapshot. The present work changes exactly this part:

- Step1, kernel-level retention: the undamped shape is kept with its exact
  initial-slip primitive. No nonlinear modal normal form is claimed.
- Step2, operator estimate: for given uniformly C2 slow test envelopes,
  the actual projected memory has the uniform O(epsilon^2) bound, including
  original-time derivatives, on0<=t<=S epsilon^-2. The continuum primitive
  is L1 despite the divergent absolute first moment of k_c.
- Actual-envelope regularity/application remains OPEN. Test functions do
  not establish that the real field trajectory is in this slow class.
  In particular a raw carrier envelope may contain small fast shape-slip
  corrections; a C1 bound alone does not pay the theorem's C2 norms.
- Step3, full nonlinear PDE control on epsilon^-2: OPEN. The source
  epsilon^2(1+T)^3 smallness window cannot be reused at that scale.
- Step4, averaged mode law/decay: NOT ACTIVATED. Old Gamma is identified
  only as the already paid susceptibility imaginary part, not a newly
  established damping law for the actual source solution.

The actual threshold counterexample prevents extending the estimate to
arbitrary carriers. Do not remove the shape slip, replace it by noise,
fit new masses, introduce incoming data or backport success into old files.
No physical/source or ontology gate is discharged by this dependency delta.
