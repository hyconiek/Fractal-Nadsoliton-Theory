# Nowe badania po SC: odrzucenie zbyt silnego lematu i jawny obiekt transportowy

Rozpoczęto: 2026-10-16; końcowa kontrola: 2026-10-17. Dowody analityczne warunkowe względem przyjętego modelu
A7/MT, nie formalizacja ani niezależna recenzja naukowa.

## Wynik

Nie zakończyłem pracy na ponowieniu OPEN. Powstały cztery obiekty dowodowe:

1. **Rzeczywisty kontrprzykład do mocniejszego lematu osłaniania.** Dla fazy
   A z usuniętym jednym rzadkim stanem parzystości iloraz proponowanego
   mocniejszego oszacowania rośnie do nieskończoności. Właściwy signed
   residual jest jednak dokładnie niedodatni. To nie kontrprzykład do SC
   lub IBA: pokazuje błąd kierunku dowodu, który usuwał potrzebne znoszenie
   przepływów. Dowód działa na rzeczywistej rodzinie liczności FIN, nie
   tylko na grafie zabawkowym; korzysta z podciągu N podzielnych przez 12000.
2. **Poprawna kontrola defektów całej fazy.** Przekroje różniące się od
   zbilansowanej fazy o ograniczoną liczbę orbit we wnętrzu mają jednolicie
   ograniczony wynik przekroju. Ogólnie wystarcza ograniczona skierowana
   głębokość defektów. Nie zakładamy, że każdy przekrój ją ma.
3. **Lokalne prawo na pojedynczym poziomie parzystości.** Przeliczenie
   rzeczywistej dziesięciowymiarowej siatki daje kowolumen 6N^-5, a z krokiem
   miękkim pełne 12N^-23/4. Pozwala to udowodnić ograniczenie przekrojów
   faza-przecięta-prefiksem w zwartym oknie krytycznym, z granicą g_P(z)/288.
   Nie rozciąga to wyniku na dalekie ogony rosnące z N.
4. **Jawne nieujemne marginesy prądu.** Z dokładnego źródła mikroskopowego
   i miar poszczególnych włókien zbudowałem prądy wejścia/wyjścia. Poprawka
   na prawdziwej warstwie wyjścia zachowuje dokładnie masę i punktowy bilans
   dryfu. W zwartym oknie marginesy są nieujemne i ograniczone przez 2NW
   oraz 3NW. To gotowy obiekt do szukania sprzężeń po dozwolonych krawędziach.

## Co pozostaje do udowodnienia

Nieujemne marginesy nie gwarantują sprzężenia. Następny dokładny warunek
to ważone nierówności Halla pomiędzy **całymi** sąsiednimi włóknami, przy
podporze H obejmującej także dozwolone przejścia między fazami przy E_N.
Jeśli sprzężenia istnieją, dają we wnętrzu badanego pasma dokładny dryf
i całkowitą szybkość co najwyżej 5N. Pozostaje też sklejenie z dalekimi
poziomami i brzegiem całej dziedziny. Nie skonstruowano pełnego generatora.

Wybór marginesów jest jawną metodą dowodową. Jego ewentualna porażka nie
oznacza niemożliwości wszystkich przepływów ani dowodu uniwersalnego IBA.
Tym bardziej nie zmienia praw podstawowego FIN. Zachowano oryginalne
progi, A7, g_c, zegar, źródło i dziedzinę; nie dopisano kąpieli, pamięci,
szumu ani kontrolera.

## Weryfikacja i przekazanie

[Dowód kontrprzykładu](ACTUAL_COUNTERCUT_TO_STRONG_SHIELDING.md),
[defekty fazy](BALANCED_PHASE_DEFECT_THEOREM.md),
[lokalne poziomy](COMPACT_PHASE_PREFIX_LOCAL_LIMIT.md),
[konstrukcja marginesów i sklejanie](EXACT_COMPACT_CURRENT_MARGINALS_AND_GLUING.md).

[Weryfikator](verify_research.py) ma 20 ograniczonych kontroli: wymierne
tożsamości przepływów, siatkę i poprawkę marginesów, dodatnie/ujemne testy
Halla oraz oddzielne nominalne kontrole A7 i modelu q2. Nie są to 20 dowodów
asymptotycznych ani numeryczne certyfikaty całego grafu. Małe przykłady nie
udają rzeczywistego obszaru metastabilnego. Wyniki są w RESULTS.json.

[Zlecenie kolejnego łańcucha](../FIN_BALANCED_PHASE_CURRENT_NEXT_CHAIN_20261016.md)
jest gotowe dla agenta: sprzężenie jawnych marginesów -> dokładny przepływ
na paśmie -> jawna bramka sklejenia z resztą D_N -> decyzja. Nie zleca
udowadniania obalonego lematu o samej masie osłoniętych stanów.

Postęp dotyczy odporności statystycznej teorii efektywnej. Nie wyprowadzono
naturalnego źródła aktywności/A7, jednostek, przestrzeni, zachowanej energii,
aparatury, fizycznego zegara, QM/GR ani ToE.
