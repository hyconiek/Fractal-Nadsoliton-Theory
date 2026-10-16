# FIN — wynik następnego etapu po RB-CUT

**Werdykt:** `SHARP_NATIVE_DOMINATED_CUT_CONSTRUCTION_PROVED; ACTUAL_Q12_BOUND_OPEN`.

Zamiast kolejnego oszacowania masy cienkiej warstwy, uzyskano dokładną **formułę najmniejszego współczynnika zwiększenia dozwolonych natywnych przewodności**, który pozwala odtworzyć cały pierwotny dryf parzystości.

Dla każdego skończonego grafu liczności `H_N` niech `s_x` będzie pierwotnym źródłem wynikającym z natywnego heat-bath, a `c^0_e` — jego rzeczywistą przewodnością. Wtedy

`K_N^* = max(0, sup_S (sum_(x in S) s_x)/ (sum_(e wychodzące z S w H) c^0_e))`,

przy czym dodatnia liczba dzielona przez zero oznacza nieskończoność. Wynika to bezpośrednio z twierdzenia o maksymalnym przepływie.

Jeżeli dla wszystkich dostatecznie dużych `N` ten wskaźnik pozostawałby ograniczony przez jedną liczbę `K`, istnieje konstrukcja nieujemnych, odwracalnych przewodności, które zachowują dokładny dryf we wszystkich stanach i spełniają ograniczenie całkowitego tempa `KN`.

**Ale takiego ograniczenia dla prawdziwego q12/A7 nie udało się udowodnić.** Otwarty pozostaje nawet znak źródeł na wszystkich nieregularnych dolnych ideałach grafu MT.

Nowe kontrole wymierne pokazują, że:

- usuniętych przejść nie można po prostu zastąpić monotonicznymi ścieżkami o tych samych końcach — skok parzystości ma zawsze jeden krok;
- można jednak kompensować ich łączny efekt;
- czasami trzeba zarówno zwiększać, jak i zmniejszać niektóre dozwolone przewodności;
- nawet rozbieżny współczynnik `K_N^*` **nie dowodzi**, że nie da się spełnić słabszego warunku całkowitego tempa `CN`.

Dowody, testy i ograniczenia są w dołączonym handoffie. Nie zmieniono prawa FIN, poprzednich pakietów ani `AGENTS.md`. Fizyczne źródło kinetyki, QM/GR/ToE pozostają nieudowodnione.
