# MICROSCOPIC-VS-METASTABLE-FIBER-SCALE-125
## A microscopic maximum-entropy fiber is exponentially faster than the metastable geometry

Date: 2026-09-26

Status:
- exponential-scale theorem within the conditional finite-N heat-bath model;
- clarifies the meaning of the competing fiber-rate source candidates.

## 1. Metastable coarse scale

For the twelve localized minima, the slow global communication is controlled
at large g by the d=4 saddle class.

Write its effective rate as

    k_4
      =
      A_4(N,g)
      exp[-N Delta V_4(g)],

where A_4 is subexponential on the N scale.

Reports 110-112 give the coarse spectral gap

    delta_meta
      ~ 3 k_4.

At leading large-g order,

    Delta V_4(g)
      =
      alpha_4 g-log 2+o(1),

with

    alpha_4=0.281135945088269.

## 2. If the fiber uses the microscopic refresh clock

The declared heat-bath microscopic event rate is an O(1) dimensionless clock
rho_micro.

If a binary refined fiber uses the same maximum-entropy refresh process,

    rho_f=rho_micro=O(1),

and report 124 gives

    2mu=rho_micro.

Therefore

    (fiber relaxation)/(coarse metastable gap)
      =
      rho_micro/delta_meta

      ~
      [rho_micro/(3A_4)]
      exp[N Delta V_4(g)].

Hence:

    boxed:
    the fiber is exponentially faster than the metastable coarse geometry.

In the metastable limit it equilibrates essentially instantaneously compared
with transitions between localized minima.

## 3. Spectral matching means something different

The conditional candidate

    2mu=delta_meta

requires

    rho_f=delta_meta.

That means the fiber itself must update on the metastable slow scale.

It cannot simultaneously be an ordinary microscopic degree of freedom updated
at the O(1) heat-bath rate unless its own updates are suppressed by a new
barrier/mechanism.

So there are two physically different interpretations:

### microscopic fiber
    rho_f=O(1);
    fiber is fast;
    c=mu/delta_meta is exponentially large;

### emergent/metastable fiber
    rho_f~delta_meta;
    fiber is slow;
    c=O(1).

## 4. Consequence

The spectral-matching rule is not merely a choice of numerical constant.

It implicitly chooses the ONTOLOGICAL LEVEL of the fiber.

If children are microscopic, spectral matching is generically wrong on the
metastable N scale.

If children are emergent metastable states, a new incidence/barrier law is
needed to explain why their rate tracks delta_meta.
