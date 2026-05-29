# P2351 strict Track-B quartic support global convexity/Jacobian identity witness

Status: global convexity certificate for the explicit quartic support fixture; no universal boundary theorem.

- `b_GB = 13152087*log(2)/(320000000*pi**2)`.
- Support function: `h(u)=1 + eps*u1**4 on S^3, eps=1/10` with `y = u1**2 in [0,1]`.
- Radii: `r_perp = -(3*y**2 - 10)/10` and `r_gradient = -(15*y**2 - 12*y - 10)/10`.
- Lower-bound certificates: `r_perp-7/10 = -3*(y - 1)*(y + 1)/10`; `r_gradient-7/10 = -3*(y - 1)*(5*y + 1)/10`.
- Determinant Jacobian: `-(3*y**2 - 10)**2*(15*y**2 - 12*y - 10)/1000`; lower bound `343/1000`.
- Pointwise density-Jacobian residual: `0`.
- Integrated replay: `Integral sigma3 dA = 2*pi**2`; boundary functional `32*pi**2`; residual `0`.
- Normalized Track-B pairing: `13152087*log(2)/10000000`; residual `0`.
- No general support-function theorem, no arbitrary-boundary theorem, no general Chern-Gauss-Bonnet theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.
