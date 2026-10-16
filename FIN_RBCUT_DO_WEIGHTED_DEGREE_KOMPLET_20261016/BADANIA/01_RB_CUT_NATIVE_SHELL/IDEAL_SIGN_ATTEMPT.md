# RB-CUT actual sign attempt — exact stopped line

The proof attempt uses the complete genuine finite source `s=A_U c0` and separates all allowed and removed edges, rather than replacing the true MT domain by a toy cylinder.

1. **Exact deleted-edge geometry:** all forbidden interior phase-crossing edges are confined to the moving Fourier-radius shell of `NATIVE_SHELL_THEOREM.md`. This uses the parity/Fourier lattice relation, not a Gaussian approximation.
2. **Exact finite downward defect control:** for every H-ideal `L`, `sum_L s=cut_up(H,L)+cut_up(F,L)-cut_down(F,L)`. Entire `cut_down(F,L)` is bounded above by `N Z^W pi_D(Shell)`. The positive terms are not known to dominate it pointwise on every irregular ideal.
3. **MT asymptotics:** stationary mass of the shell tends to zero; hence every normalized ideal violation is at most a vanishing amount. This is NOT a proof of positivity of the exact source for any finite N.
4. **Attempted repair:** deleting forbidden native edges yields nonnegative reversible conductances with uniform rate cap 1. But the leftover signed source `A_F c_F` is shell-supported and has total variation at most `2N Z^W pi_D(Shell)`. An exact flow repair would have to solve `A_H delta_c=A_F c_F`, with `c_H^0+delta_c>=0`, on the unchanged actual graph, while respecting one uniform cap. The mass estimate alone does not guarantee existence of such a signed correction or nonnegativity.
5. **Exact validity obstruction:** no injection pairing every forbidden downward-crossing native edge with admissible H-upward traffic was constructed. No actual negative H-ideal on the true `D_N` was constructed either. A synthetic rational witness from RB-002 cannot be promoted to q12.

**Terminal disposition**: `OPEN_ACTUAL_H_IDEAL_SIGN_WITH_BOUNDARY_SHELL_LOCALIZATION`.

The next *exactly formulated* atom is to control the sign of every true H_N-lower ideal with the known frozen Fourier-shell geometry, by a legal native-weight pairing / flow reconstruction / exact negative ideal certificate. Neither a small N LP nor a surrogate threshold resolves the requested all-large-N problem. The independent uniform rate cap awaits a positive sign result.
