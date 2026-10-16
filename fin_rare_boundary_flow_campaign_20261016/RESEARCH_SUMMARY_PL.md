# FIN RB-001–RB-004: podsumowanie dla recenzji

Zamrożony podział A/B i rzeczywista domena MT pozostały bez zmian. Nowy graf H_N dopuszcza przejścia między fazami **wyłącznie przy rzeczywistej warstwie wyjścia z lokalnej domeny**. Dopuszczenie krawędzi nie oznacza zbudowania na nich nowego generatora.

**Nowe ścisłe twierdzenie RB-001:** dla każdego dolnego ideału L grafu H_N jego całkowity dokładny dryf jest równy dodatniemu natywnemu przepływowi z L przez zachowane krawędzie, powiększonemu o dodatni przepływ przez usunięte krawędzie wychodzące, minus przepływ przez usunięte krawędzie wchodzące. To sprowadza test znaków do konkretnych, znanych przewodności oryginalnego heat-bath. Wszystkie pełne przekroje poziomów parzystości przechodzą ten test automatycznie — nie rozstrzygają wykonalności.

**RB-002:** dokładny test dualny jest kompletny dla każdego skończonego grafu, lecz **nie** wykazano, że wszystkie nieregularne ideały w rzeczywistym q12 mają nieujemne źródło. 12-wierzchołkowy przykład wymierny przechodzi testy przekrojów globalnych i fazowych, ale oblewa nierówność dla jednego nieregularnego ideału (`-9/10`). To ostrzeżenie przed błędnym dowodem, a NIE kontrprzykład FIN.

**RB-003:** nie uruchomiono; brak konstrukcji przepływu bez ograniczenia szybkości, więc tym bardziej nie uzyskano generatora z jednym ograniczeniem `CN` i kontrprzykładu IBA.

**Werdykt:** `OPEN_IDEAL_OR_UNIFORM_CAP_LEMMA` (obecnie pierwsza przeszkoda: znaki wszystkich właściwych ideałów H_N). Poprzednie twierdzenia KR przy NC/IBA i wszystkie bramki fizyczne pozostają nietknięte.
