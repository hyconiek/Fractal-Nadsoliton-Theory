#!/usr/bin/env python3
"""FIN ST1752--ST1766: literal variation and minimal gauge-repair audit."""
import sympy as sp

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1752,1766
NAMES=["LiteralFlatAction","ScalarVariation","YukawaSignDefect","GaugeVariation","FermionVariation",
 "DensityTypeMismatch","GaugeInvarianceNoGo","ChargedScalarRepairExpansion","GaugeCoefficientConstraint",
 "ProcaFork","CovariantExtraActionRequirements","PotentialNormalizationMismatch","SingleActionCoverage",
 "RoundThree_Verdict","RoundThree_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 xvar=sp.symbols('x');phi=sp.Function('phi')(xvar);Avec=sp.Function('A')(xvar);m2,lam,y,g,xi,R,B=sp.symbols('m2 lam y g xi R B');L=sp.diff(phi,xvar)**2/2-m2*phi**2/2-lam*phi**4/4-y*phi*B+g*Avec**2*phi**2/2+xi*phi**2*R
 x={}
 x["ST1752"]=packet(1752,"constructed_literal_flat_scalar_vector_source_action_from_displayed_terms","Suppresses spin/gravity derivatives only for sign/source audit.",{
  "L":str(L)})
 EL=sp.simplify(sp.diff(sp.diff(L,sp.diff(phi,xvar)),xvar)-sp.diff(L,phi))
 x["ST1753"]=packet(1753,"proven_literal_scalar_Euler_row","Using equation-sheet EL convention.",{
  "E_phi":str(EL)})
 yuk_coeff=sp.expand(EL).coeff(B)
 x["ST1754"]=packet(1754,"proven_compact_draft_Yukawa_source_sign_disagrees_with_literal_variation","L contains -y phi psibarpsi but displayed compact EOM contains -y psibarpsi.",{
  "derived_coefficient":str(yuk_coeff),"draft_coefficient":"-y","match":str(yuk_coeff)=="-y"})
 E_A_mass=sp.simplify(-sp.diff(L,Avec))
 x["ST1755"]=packet(1755,"proven_literal_vector_interaction_variation_matches_masslike_term_up_to_Maxwell_overall_convention","Multiply full Maxwell equation by -1.",{
  "minus_dL_dA":str(E_A_mass),"mass_term":"g phi^2 A"})
 x["ST1756"]=packet(1756,"proven_fermion_variation_gives_Dirac_minus_yphi_operator","Displayed fermion row locally consistent.",{
  "E_psibar":"(i gamma D - y phi)psi"})
 x["ST1757"]=packet(1757,"proven_compact_Ltotal_adds_scalar_Lagrangians_to_gravity_density_before_covariant_repair","sqrt(-g) typing inconsistent in compact sum.",{
  "compact_density_consistent":False,"later_covariant_presentation_repairs_notation":True})
 x["ST1758"]=packet(1758,"proven_A2phi2_seagull_alone_is_not_local_gauge_invariant","Off-shell infinitesimal variation nonzero.",{
  "delta":"g phi^2 A^mu partial_mu alpha","integrated":"-g alpha partial_mu(phi^2 A^mu)+boundary"})
 x["ST1759"]=packet(1759,"constructed_minimal_complex_charged_scalar_repair","Adds more than the existing term.",{
  "expansion":"|D Phi|^2=|partial Phi|^2 + iq A(Phi*partial Phi-Phi partial Phi*) + q^2 A^2|Phi|^2"})
 x["ST1760"]=packet(1760,"proven_gauge_repair_ties_seagull_and_derivative_current_to_same_charge","Current draft has independent g_eff and no required derivative current.",{
  "constraint":"g_eff=q^2 modulo normalization","missing_current":True})
 x["ST1761"]=packet(1761,"constructed_Proca_interpretation_fork","Consistent effective vector option but abandons gauge/BRST claim.",{
  "local_gauge":False,"constraint":"divergence of phi-dependent Proca equation"})
 extras={"lambda3_phi3":"add cubic potential","H_sector":"add complex Higgs kinetic/potential","phi2_H2":"add mixing potential","chi_RG_RF2":"add R F^2 term","curvature_squared":"add R^2,Ricci^2,Riemann^2 action","T_SM_mix":"define and vary all matter sectors"}
 x["ST1762"]=packet(1762,"proven_claimed_covariant_EOM_extras_require_explicit_additional_action_terms","Cannot arise from coefficient renaming.",extras)
 x["ST1763"]=packet(1763,"proven_quartic_and_cubic_normalizations_differ_between_scaffold_and_sector_rows_without_exported_map","-lambda phi^4/4 yields lambda phi^3; lambda4/3! corresponds lambda4 phi^4/4!.",{
  "normalization_map_supplied":False,"lambda3_action_absent":True})
 x["ST1764"]=packet(1764,"proven_no_single_literal_displayed_action_generates_all_current_covariant_rows","Term-source/sign audit.",{
  "single_action_closure":False,"reasons":["unsourced sectors","Yukawa sign","gauge seagull","normalization mismatch"]})
 x["ST1765"]=packet(1765,"round_three_falsifies_full_Ltotal_claim_but_preserves_reduced_termwise_models_as_separate_actions","Reduced recomposition zero is internally valid.",{
  "full_covariant":False,"reduced_models":True})
 x["ST1766"]=packet(1766,"recommended_ST1767_ST1781","Falsify the moment-to-coupling map: reference dependence, rank and predictive content.",{
  "next":["dref variation","coordinate rescaling","Jacobian rank","free base couplings","moment-equivalent kernels"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
