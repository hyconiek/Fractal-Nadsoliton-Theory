# P1579 S529 Strict rho_gr Derivation And Internal Replication Packet (No Legacy Bridge)

Status: `P1579_EXECUTED_RHO_GR_DERIVATION_AND_W1578B_REPLICATION`
As of: `2026-05-14`

## Cel

Po `P1578` formalizujemy:

1. `T1578A` — robocze wyprowadzenie `rho_gr(p)` z lokalnego wskaźnika ryzyka,
2. `W1578B` — wewnętrzną replikację zgodności z pełnym łańcuchem EOM.

## Konstrukcja strict-only

- Definiujemy ryzyko lokalne:
  `risk(p) = ||J^{-1}(p)||_inf * gain_eom(p)`.
- Definiujemy tłumienie:
  `rho_gr(p) = min(1, kappa / risk(p))`.
- Replikacja:
  - metoda A: jawna formuła,
  - metoda B: niezależna rekonstrukcja z warunku docelowego `R_bundle<=1`.

## Kryterium PASS/FAIL

- `PASS_W1578B_REPLICATION_CANDIDATE` gdy metody A i B dają zgodne `rho_gr`
  (błąd < tolerancja) i utrzymują kandydacką stabilność bundle.

## Wynik

`PASS_W1578B_REPLICATION_CANDIDATE`.

## Brakujące obiekty do final strict-core closure

1. `T1579A`: theorem fizycznej interpretacji parametru `kappa`.
2. `T1579B`: theorem globalnej jednoznaczności `rho_gr` na domenie strict.
3. `T1579C`: końcowy theorem ToE-bundle closure.

## Rekomendacja następnego uczciwego kroku

Uruchomić `P1580`: semiglobalna walidacja jednoznaczności `rho_gr` i eksport
`T1579B`.

## Omówienie dla laika

To krok „sprawdź dwa razy”: ten sam tłumik GR policzyliśmy dwiema różnymi
metodami i dostaliśmy zgodny wynik, więc mamy mocniejszą pewność, że to nie
przypadek liczbowy.
