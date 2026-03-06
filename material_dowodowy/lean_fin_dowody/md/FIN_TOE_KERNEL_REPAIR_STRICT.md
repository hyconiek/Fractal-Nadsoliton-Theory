# FIN-ToE Strict Kernel Repair (Nadsoliton Characteristics First)

Data: 2026-03-02
Status: roboczy, bez zmian w TeX

## 1) Zasada przewodnia
Najpierw sa charakterystyki nadsolitona, potem dopiero postac kernela.
Nie odwrotnie.

## 2) Minimalny, rygorystyczny ledger statusow
- [D] alpha_geo = 4 ln 2 (2.7725887222): silna kotwica geometryczno-informacyjna.
- [A] omega = pi/4: obecnie ansatz symetrii oktawowej, nie zamkniete wyprowadzenie mikrodynamiczne.
- [A] phi = pi/6: obecnie ansatz symetrii heksagonalnej, nie zamkniete wyprowadzenie mikrodynamiczne.
- [C] beta_tors ~ 0.01: stabilna kalibracja robocza, ale brak derywacji first-principles.
- [M] 1/(1+beta d): mechanistycznie mozliwe przez path-sum/topologie, ale obecnie tuning-sensitive.

## 3) Co trzeba uznac za niespojne (do odciecia narracyjnego)
1. Jednoczesne uzywanie phi=0 i phi=pi/6 jako "kanonicznych".
2. Jednoczesne uzywanie beta_tors=0.01 i beta_tors=0.05 jako "zamrozonych" baz.
3. Traktowanie wzorcow wezlow [2,5,8,11] i [2,8,14] jako rownowaznych.
4. Twierdzenie, ze dla (omega,phi)=(pi/4,pi/6) wezly sa co 3 oktawy.

## 4) Poprawiona, spojna wersja kernela (vR1)
Przyjmujemy (tymczasowo) gałąz zgodna z charakterystykami nadsolitona:

K(d) = alpha_geo * cos(omega d + phi) / (1 + beta_tors d)

z:
- alpha_geo = 4 ln 2,
- omega = pi/4,
- phi = pi/6,
- beta_tors wstepnie 0.01 (do niezaleznej derywacji).

Konsekwencja matematyczna (obowiazkowa):
- zera cos: d_n = (pi/2 - phi + n pi)/omega = 4/3 + 4n.
- To NIE daje wezlow [2,5,8,11].
- Jesli ktos chce wezly [2,5,8,11], musi porzucic co najmniej jedno z: omega=pi/4, phi=pi/6.

## 5) Dwie rozlaczne galezie (bez mieszania)
- Galaz A (charakterystyki-first): utrzymuje (pi/4, pi/6), poprawia narracje wezlow.
- Galaz B (wezly-first): dopasowuje (omega,phi) do wezlow, ale wtedy odchodzi od obecnej osi charakterystyk.

Do czasu badan 1734-1737 obie galezie nie moga byc laczone w jednym ciagu "dowodu".

## 6) Wyniki badan wykonanych teraz
- QW-1729: mapa charakterystyk->kernel, score=0.522, not closed.
- QW-1730: chrono-audit, high-risk inconsistency (risk_points=9).
- QW-1731: node compatibility, strongly inconsistent.
- QW-1732: path->hyperbolic, plausible but tuned.
- QW-1733: closure gate, hard FAIL, readiness open.

## 7) Warunek domkniecia ToE dla kernela
Kernel uznajemy za domkniety dopiero gdy jednoczesnie:
1. (omega,phi,beta_tors) sa wyprowadzone, nie tylko dobrane.
2. Wzorzec wezlow jest jeden, statystycznie preferowany i zgodny z priorem charakterystyk.
3. Ten sam zestaw parametrow przechodzi wspolnie flavor + GW + OOS stabilnosc.

## 8) Powiazane artefakty
- report_qw1729_nadsoliton_kernel_characteristics_map.json
- report_qw1730_nadsoliton_kernel_chrono_audit.json
- report_qw1731_nadsoliton_kernel_node_compatibility.json
- report_qw1732_fractal_path_to_hyperbolic_test.json
- report_qw1733_nadsoliton_kernel_closure_gate.json
