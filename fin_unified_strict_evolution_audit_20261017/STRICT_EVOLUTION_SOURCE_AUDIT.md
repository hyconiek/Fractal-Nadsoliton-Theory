# Existing FIN evolution: exact sources, not an assumed fundamental completion

Date: 2026-10-17. This implements Stage A of the expanded objective. The
three named QW programs, their archived JSON/Lean artifacts, all six required
kernel/FAR instruction notes, relevant current state-map and nonlinear-source
audits were read. AGENTS was indexed across its full ledger; the governing
source/action/geometry/nonlinearity sections were inspected verbatim. This
does not claim re-verification of every historical programme.
Original sources and AGENTS remain byte-preserved.

## 1. The actual QW-2163/2165 action

The scripts construct, rather than derive from a kernel, the real polynomial
family with q=(psi_0,...,psi_11,Phi):

    L=1/2 sum_a(q'_a)^2-V(q),
    V= m2_Phi Phi²/2+lambda_Phi Phi⁴/4
       +sum_i[m2_i psi_i²/2+g4_i psi_i⁴/4+g6_i psi_i⁶/6]
       +sum_i y_i Phi² psi_i²+1/2 psi^T K_s psi,
    K_s=(K+K^T)/2, K_ii=0.

All coefficients are independent real SymPy symbols. The2163/2165 field
functions are formal functions without a complex-conjugate variational
calculus;2166 uses real field symbols. The real-valued specialization is
the one audited here, not a derivation that the full legacy fields are real.
The scripts use one real independent variable x, not an exported spacetime manifold, metric or
causal Cauchy surface. The exact Euler--Lagrange equations are

    psi_i''+m2_i psi_i+g4_i psi_i³+g6_i psi_i⁵
       +2y_i Phi² psi_i+sum_j(K_s)_ij psi_j=0,
    Phi''+m2_Phi Phi+lambda_Phi Phi³+2Phi sum_i y_i psi_i²=0.

QW-2163 checks three representative equations; QW-2165 checks all13.
These are real mathematical results about a **supplied template**. Neither
script binds K_ij to the strict radial kernel, derives the nonlinear
coefficients from it, or identifies x with a physical clock. No arbitrary
first-order replacement is made here: the equivalent first-order state
has all13 coordinates and13 momenta.

## 2. What QW-2166 proves and what it does not

Its potential Hessian has entries

    H_ii=m2_i+3g4_i v_i²+5g6_i v_i⁴+2y_i v_Phi²,
    H_ij=(K_s)_ij for distinct psi indices,
    H_i,Phi=4y_i v_i v_Phi,
    H_Phi,Phi=m2_Phi+3lambda_Phi v_Phi²+2sum_i y_i v_i².

The script derives eta''+H(v)eta=0 for its quadratic fluctuation template.
It does not check grad V(v)=0 before calling v a vacuum. A constant v with
nonzero gradient is not a reference solution of the full equations; its
tadpole cannot be discarded. Symmetry of H alone is not positivity,
stability, continuum propagation or a physical mass identification.

The three reports explicitly keep their full-action/full-continuum theorem
flags **false**. Their Lean files prove conjunctions from supplied Prop
hypotheses. That logical wrapper is not a Lean formalization of polynomial
variation, well-posedness or source provenance. We independently reconstruct
the actual symbolic action, all13 equations and Hessian entries; no historical
writer is run against the source repository.

## 3. The legacy field action is a different typed object

The legacy TeX, lines289 onward and1814 onward, supplies12 **complex**
scalar Psi fields, a real Phi, a spacetime derivative and polynomial
potentials. It is not the same state type as13 real one-coordinate fields.
Before a real invariant slice, it has25 real configuration fields, not13.

For real symmetric K, realification Psi_i=u_i+i v_i gives a real action
in25 fields, with radial polynomials in u_i²+v_i². The slice v_i=dot v_i=0
is invariant and restricts to the QW potential after explicit coefficient
and sign mapping. A complex Hermitian mixing with nonzero imaginary part
need not preserve that slice. Deleting the imaginary fields without these
conditions would lose dynamics and possible charges.

There are also literal convention differences: the legacy psi mass term
is negative, while QW m2_i is an unrestricted signed coefficient. One TeX
section has sextic delta/6, another delta/8. The respective maps are
m2_i=-m_i² and g6_i=delta or3delta/4, not one raw identity. Generation-shared
couplings versus independently indexed y_i require a stated map. These are
source/convention obligations, not licence to select the most useful version
after seeing a particle-like solution.

## 4. Kernel split and current FAR status

K_legacy=alpha_geo cos(omega_L d+phi_L)/(1+beta_tors d) remains an intermediate
legacy object. The strict working tuple is(0.18575,.16250,1,9/5), not the
legacy tuple. The required K1/K2/F2/F3/S2 notes and later state-map guardrails
forbid silent identification and legacy role transfer. Later state-map rules
supersede older generic bridge-first recommendations; no closed selector,
unit or lower-boundary lane is rerun here.

The current STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md contains real
conditional scalar/fermion/gauge/gravity formulas. Its own promotion table
calls moment coefficients effective/scaffold inputs and does not export
role-bearing L_total. P2685/P2697 retain nonproxy reverse-closure and freshness
limits. P3075 exports an internal dimensionless Dirichlet gradient source,
not spacetime Hamiltonian dynamics. P507/P517 already show that the quadratic
core and information language do not select the attractive fourth-order jet.
These results are reused, not rediscovered or upgraded by the present audit.

## 5. What is missing for one strict T_t

Once the QW coefficients and interpretation are supplied, it defines a
well-posed local finite-dimensional family, as proved in the composition
document. What is not exported is the **source map** from strict nadsoliton
data to this full action/state/domain, its coefficient selection, admissible
Cauchy data and inter-object/field composition. The FAR map
lambda_eff=lambda(1+c1²) retains the free base lambda even after every strict
moment is fixed; it does not source that freedom.

The q12/A7 heat bath is another declared parent. Neither a typed dynamical
projection from this full field state to occupation probabilities nor a
derivation of its kinetic clock/noise from the full action is supplied by
QW-2163/2165/2166. The coincidence of13 real fields and a twelve-state-plus-
intermediate comparison is not such a map. The fundamental closed-world
claim therefore remains underdetermined on these actual inputs.
