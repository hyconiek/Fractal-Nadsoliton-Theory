# Jawna kanoniczna mapa modów z zachowaniem całego kontinuum

Data: 2026-10-21. Warunkowy dowód analityczny, audyt własny.

Wyprowadziłem mapę współrzędnych z rzeczywistego Hamiltonianu zadanej
rodziny RP. Usuwa nierezonansowe sprzężenia modów związanych w rzędzie
trzecim oraz ich zależne od fazy wyrazy czwartego rzędu. Nie usuwa
kontinuum, jego pędów ani sprzężenia odpowiedzialnego za promieniującą
odpowiedź. Nie jest to zastąpienie PDE czterema zmiennymi.

Generator sześcienny ma jawną formułę; generator czwartego rzędu daje
odwracalny układ19-wymiarowy. Dowód obejmuje cały dotychczasowy przedział
parametru. Przekształcono również prawdziwe dane początkowe. Przy środkowym
a0=-1/2 swobodny składnik modu kształtu w drugim rzędzie znosi się, ale
wymuszona odpowiedź pozostaje. Nie wolno mylić tych dwóch stwierdzeń.

[Pełny dowód](PROOF_OR_OPEN_LEMMA.md) i
[weryfikator](verify_research.py) zawierają dokładny rachunek wymierny,
w tym tożsamości symboliczne w parametrze. Kontrole nie są dowodem
istnienia powolnej amplitudy rzeczywistego PDE na czasie epsilon^-2.
Pozostaje oszacowanie mieszanych sprzężeń z kontinuum i reszt przestrzennych;
[następny konkretny brak](STATE_MAP_AND_REMAINING_GATE.md) jest jawny.

Nie dodano prawa ewolucji, szumu, kąpieli, przygotowania ani zegara.
Stare źródła i AGENTS zachowano. F01 oraz fizyczne bramki QM/GR/ToE
pozostają otwarte: mapa założonej dynamiki nie jest jej naturalnym źródłem.
