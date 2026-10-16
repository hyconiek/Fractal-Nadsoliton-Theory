# FIN SC-001–SC-004 — podsumowanie po polsku

**Wynik:** `OPEN_SIGNED_INTERFACE_LEMMA`. Nie udało się ani udowodnić jednolitej nierówności dla wszystkich rzeczywistych nieregularnych przekrojów orbitowych, ani skonstruować rosnącego kontrprzykładu q12. To nie jest kolejna rzekomo zamknięta bramka przepływów.

**Nowy, ścisły wynik lokalny:** w zamrożonym FIN q12/A7 każdy koniec usuniętej krawędzi między fazami leży w cienkiej powłoce Fourierowskiej, ale poza właściwą warstwą odejścia MT. Cztery zmiany całkowitych sum `(C,S)` przy ruchu w każdą stronę parzystości to `(±1,±1)`; każda jest realizowana przez dziewięć różnych zamian etykiet. Dla dostatecznie dużego N dziewięć ruchów w kierunku wewnętrznym (faza A) albo zewnętrznym (faza B) na pewno pozostaje w tej samej fazie. Są to zatem **rzeczywiste legalne przejścia H**, nie nowe krawędzie.

Otrzymano dokładne oszacowanie dla dowolnego zbioru orbit V:

`Fout(V)-Uin(V) <= N B(V)+N W(Z_F(V))-Uin(V)`.

Tu `Z_F` to stany wysyłające zabroniony strumień poza V, których wszystkie dozwolone ruchy wzrostu parzystości prowadzą do wnętrza V. Dla przekrojów bez takich stanów nierówność kompensacji SC zachodzi już ze stałą `L=1`.

**Otwarta przeszkoda:** dozwolone lokalnie ruchy mogą być *wewnętrzne* względem nieregularnego przekroju. Nie wiadomo, czy ich masę można kontrolować strumieniem natywnym wchodzącym do przekroju i liczbą dozwolonych krawędzi wychodzących, jednolicie względem N. Bez tego SC-002 i SC-003 pozostają nieuruchomione.

**Pochodzenie danych:** starsze pliki kampanii i archiwalna (nieaktualna) wersja AGENTS były dostępne i zahaszowane. Najnowsze pliki przeglądu SC C01–C03 oraz innych przeglądów nie były dostępne jako bajty; nie deklarujemy ich niezależnego odczytu.

Pakiet zawiera dowód matematyczny, uściślony brakujący lemat, 148 ograniczonych kontroli (większość całkowitych/wymiernych, część przykładowo zmiennoprzecinkowych), wszystkie statusy i manfest. Testy nie są 148 niezależnymi dowodami ani certyfikacją całej rodziny N. Nie zmieniono FIN, poprzednich kampanii ani AGENTS; status fizycznego źródła aktywności i teorii fundamentalnej pozostaje otwarty.
