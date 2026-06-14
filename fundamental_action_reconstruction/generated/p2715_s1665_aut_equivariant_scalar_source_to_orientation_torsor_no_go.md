# P2715/S1665 Aut-equivariant scalar-source to orientation-torsor no-go

Status: `P2715_AUT_TRIVIAL_SCALAR_SOURCE_TO_ORIENTATION_TORSOR_NO_GO`

## Finite equivariance result
- one_point_equivariant_maps=0
- two_point_trivial_branch_equivariant_maps=0

## Scalar source matrix
- `entropy_or_uv_scale_scalar`: breaks_orientation_torsor=False. Aut-trivial scalar data can map equivariantly to the orientation torsor only through an Aut-fixed torsor point; P2714/P2715 find none.
- `alpha_geo_amplitude_normalization_scalar`: breaks_orientation_torsor=False. Aut-trivial scalar data can map equivariantly to the orientation torsor only through an Aut-fixed torsor point; P2714/P2715 find none.
- `positive_beta_or_z_beta_damping_scalar`: breaks_orientation_torsor=False. Aut-trivial scalar data can map equivariantly to the orientation torsor only through an Aut-fixed torsor point; P2714/P2715 find none.
- `scalar_lagrangian_density_or_metric_residual`: breaks_orientation_torsor=False. Aut-trivial scalar data can map equivariantly to the orientation torsor only through an Aut-fixed torsor point; P2714/P2715 find none.

## Decision
P2715 tests whether ordinary strict scalar data can break the P2714 orientation torsor.  The finite equivariance calculation finds zero Aut-equivariant maps from a trivial scalar source domain to the orientation torsor, because orientation-reversing units 7 and 11 require an Aut-fixed torsor point and none exists.  Entropy/UV scale, alpha_geo amplitude, positive beta damping, and scalar Lagrangian data therefore remain orientation-blind unless a new strict pseudoscalar/chiral source is supplied.

## Next honest step
Do not seek the missing sign in Aut-trivial scalar quantities.  The next admissible move must supply a genuinely strict inversion-odd pseudoscalar/chiral source with a nonzero signed value coupled to the orientation torsor, or pivot to a different typed object outside the closed lanes.  Without that, preserve the P2697-P2715 no-new-live-frontier certificate.
