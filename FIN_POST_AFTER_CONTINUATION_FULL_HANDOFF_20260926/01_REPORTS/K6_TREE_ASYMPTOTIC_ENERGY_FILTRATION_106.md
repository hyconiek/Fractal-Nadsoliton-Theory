# K6-TREE-ASYMPTOTIC-ENERGY-FILTRATION-106
## The support 2–6 branch tree is ordered simultaneously by energy slope and Morse index

Date: 2026-09-26

Status:
- direct consequence of the large-g energy expansion for the five support
  endpoints in reports 101-104;
- this ordering is asserted only for this branch tree, not for all 76 support
  classes.

For any isolated large-g support,

    V_g
      =
      -(mu/2) g
      +D(p*||uniform)
      +O(1/g).

For the five support strata attached to the k6/k4 tree, the leading data are:

    support 2:
      V_g ≈ -0.280888322492 g + 1.791759469228
      index 1

    support 3:
      V_g ≈ -0.229777844729 g + 1.398273263269
      index 2

    support 4:
      V_g ≈ -0.143415269409 g + 1.098612288668
      index 3

    support 5:
      V_g ≈ -0.118073179261 g + 0.878388483861
      index 4

    support 6:
      V_g ≈ -0.097590918381 g + 0.693147180560
      index 5.

Thus within this tree:

    boxed:
    mu_2 > mu_3 > mu_4 > mu_5 > mu_6,

so for sufficiently large g

    V_2 < V_3 < V_4 < V_5 < V_6.

At the same time

    index = 1,2,3,4,5.

Therefore lower-support states are both:
- energetically deeper at leading O(g);
- lower in Morse index.

Higher-support branches sit progressively higher in the stationary landscape
and possess progressively more unstable covariance directions.

This gives the k6 branch tree a genuine asymptotic filtration:

    deep / simple
        support 2, index 1
          <
        support 3, index 2
          <
        support 4, index 3
          <
        support 5, index 4
          <
        support 6, index 5
    high / complex.

The stable one-label minima lie below this entire saddle tree.

Important boundary:
the full 76-class support atlas does not obey a single total order solely by
support size; ranges from different support geometries can overlap.  The
monotone result here is specific to this connected symmetry-organized branch
tree.
