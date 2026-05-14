# P1570 S520 Strict Closed-Form Jacobian/Hessian And Stability-Bound Packet (No Legacy Bridge)

Status: `P1570_EXECUTED_CLOSED_FORM_EXPORT_AND_LOCAL_STABILITY_BOUND`
As of: `2026-05-14`

## Cel

Po `P1569` wykonujemy formalny eksport `E1569A` i start `T1569B`:

1. jawny analityczny Jacobian mapy 4-obserwowalnej,
2. jawny analityczny Hessian skalaru `xi_phase_curv`,
3. lokalny bound stabilności inwersji oparty o `||J^{-1}||` i residual obserwowalny.

## Mapa strict-only

`F(p) = (lambda_sm_eff, kappa_gr_eff, epsilon_mix_eff, xi_phase_curv)`
na `p=(omega,phi,beta,eta)`.

Eksport nie używa legacy bridge.

## Kryterium PASS/FAIL

- `PASS_E1569A`: analityczny Jacobian zgadza się z różnicami skończonymi
  w tolerancji numerycznej.
- `PASS_T1569B_CANDIDATE`: lokalny bound `||delta p|| <= ||J^{-1}|| * ||delta F||`
  pozostaje skończony i dodatni w punkcie strict.

## Wynik

`PASS_E1569A_AND_PASS_T1569B_CANDIDATE`.

## Brakujące kolejne obiekty do final strict-core closure

1. `T1570C`: theorem o jednolitym boundzie stabilności na rozszerzonej domenie.
2. `T1570D`: patching lokalnych kart inwersji do semiglobalnego atlasu strict.
3. `W1570E`: witness kompatybilności boundu z pełnym torem `L_SM + L_GR`.

## Rekomendacja następnego uczciwego kroku

Uruchomić `P1571`: semiglobalny test stabilności (`T1570C`) na siatce domeny
admissible, z raportem najgorszego punktu uwarunkowania.

## Omówienie dla laika

To krok „inżynierski”: nie tylko wiemy, że lokalnie da się odwrócić mapę, ale
też mamy wzór, który mówi jak błędy pomiaru przekładają się na błąd parametrów.
To kluczowe, by model był użyteczny fizycznie, a nie tylko formalnie.
