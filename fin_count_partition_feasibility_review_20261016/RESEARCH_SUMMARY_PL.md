# FP po audycie: konieczny limit aktywności i bardziej konkretny kierunek

Data:2026-10-16. [Pełny intake](INTAKE_20261016.md).

Kampania wykonała poprawną próbę, ale nie skonstruowała szybkości i nie
wykazała niewykonalności zamrożonego podziału. FP-003 pozostaje nieaktywne;
uniwersalne IBA nadal jest otwarte. Nie wolno zamieniać tego wyniku w PASS
przez wykorzystanie samej dodatniości prądów granicznych.

Nowy przyjęty rezultat jest warunkowy: jeżeli takie zamknięte fazy istnieją,
to graniczny profil ich średniej kwadratowej wariacji jest jednoznaczny.
Wymusza to C>=C_star, liczbowo około .64702198747. Niezależnie wyprowadziłem
ścisłe ograniczenie wymierne 23/36<C_star<47/72. Nie certyfikuję wszystkich
cyfr wydrukowanych przez kwadraturę. Ten skończony warunek konieczny nie
rozwiązuje problemu istnienia przepływów; dopuszczalne C nie ma ustalonego
z góry tak małego limitu.

Sprawdzono48/48 wpisów manifestu,7/10 pierwotnych odniesień i cztery
dodatkowe składniki archiwów. Brak dawnych opakowań BA/KR i AGENTS
pozostaje jawnym ograniczeniem pochodzenia. Replay daje identyczny wynik.
Producent ma27 kontroli mieszających rachunek dokładny z diagnostyką
numeryczną; własny weryfikator ma23 ograniczone kontrole. To nie liczby
niezależnych dowodów ani symulacje krytycznego q12.

Żeby nie zlecać ponownie tego samego nierozstrzygniętego lematu, dodałem
[nową przygotowaną drogę](BOUNDARY_LEAKAGE_AND_IDEAL_CUT_SEED.md):

- Sam gładki podział promienia Fouriera nie blokuje ruchu parzystości
  w żadnym kierunku poza rzeczywistym brzegiem domeny MT.
- Dokładne zamknięcie części jest silniejsze niż potrzeba. Gdyby udało się
  skonstruować szybkości z przeciekiem tylko przy wykładniczo rzadkim brzegu,
  oczekiwana liczba zmian części w czasie sqrt(N) nadal zanikałaby, a test
  IBA zachowałby dodatnie dolne ograniczenie.
- Przed ograniczeniem szybkości należy zbadać dokładne cięcia typu
  dolnego ideału. Dodatnie przekroje samej parzystości nie wystarczają.

**Takich nowych szybkości również jeszcze nie skonstruowano.** To jawnie
zmieniony pomocniczy test w tej samej klasie, nie ukończone stare FP i nie
nowe prawo FIN. Następny agent otrzymuje [łańcuch RB](../FIN_RARE_BOUNDARY_FLOW_NEXT_CHAIN_20261016.md):
graf i źródło -> wszystkie istotne cięcia -> jednostajny limit szybkości
-> konsekwencja czasowa dopiero po konstrukcji. Zachowano stare źródła;
aktualizacja AGENTS obejmuje tylko poprawne zakresy i otwarte bramki.

Nie wyprowadzono niezależnego źródła aktywności/A7, jednostek, przestrzeni,
zachowanej energii, aparatury, QM/grawitacji ani ToE. Brzeg MT pozostaje
obiektem analizy, nie fizyczną ścianą lub sterownikiem.
