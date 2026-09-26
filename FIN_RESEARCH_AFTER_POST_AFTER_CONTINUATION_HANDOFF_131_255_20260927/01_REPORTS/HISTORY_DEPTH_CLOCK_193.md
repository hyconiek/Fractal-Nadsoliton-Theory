# HISTORY-DEPTH-CLOCK-193
## The reversible natural extension supplies an exact additive transformation-depth clock, but not a physical duration scale

Date: 2026-09-26

Status:
exact theorem for the natural extension constructed in report 188.

Let sigma be the invertible two-sided history shift.

For two points on one shift orbit,

    omega_b = sigma^k omega_a,

define the oriented history depth

    Delta n(omega_a,omega_b)=k.

## 1. Exact additivity

If

    omega_b=sigma^k omega_a,
    omega_c=sigma^m omega_b,

then

    omega_c=sigma^(k+m) omega_a

and therefore

    boxed:
    Delta n(a,c)
      =
      Delta n(a,b)+Delta n(b,c).

So transformation depth is exactly composable.

It is independent of which state labels happen to be observed along the path.

## 2. Reversal

Because sigma is invertible,

    Delta n(b,a)
      =
      -Delta n(a,b).

Thus the global reversible history has an oriented coordinate but no preferred
thermodynamic direction.

Choosing increasing n as "future" is an orientation convention unless a
boundary/preparation condition breaks the reversal symmetry.

## 3. Reparametrization

Any positive constant a gives another additive duration

    tau=a Delta n.

All path/order statements are unchanged.

Therefore the natural extension determines:

    order + integer depth

but leaves one global positive scale unfixed.

This is precisely the clock gauge already identified in PHY-009 and reports
34/40.

## 4. Result

The new history construction PASSES:
- order;
- additivity;
- composability.

It FAILS to manufacture:
- seconds;
- an intrinsic arrow;
- a preferred positive conversion factor from update depth to physical time.

So it upgrades FIN from "no time object" to a rigorous DIMENSIONLESS history
clock, not to a calibrated physical clock.
