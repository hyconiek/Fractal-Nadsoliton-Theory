# FIN — KI-001–KI-005: wyniki i ograniczenia

Przeprowadzono lokalną kampanię analityczną dla trzech **zadanych** kinetyk odwracalnych, zachowujących tę samą miarę Gibbsa FIN: heat-bath, Metropolisa i Barkera. Założenia źródłowego A7, folda i zegara modelowego pozostają niezależnymi przesłankami.

1. **To samo ukształtowanie folda.** W pełnym stanie zajętości występuje jeden kierunek miękki i dziesięć stabilnych. Odpowiednio znormalizowane współczynniki potencjału spełniają dokładnie a_F=a i b_F=b.
2. **Różne czasy i natężenia fluktuacji.** Mobilności obliczone z mikroprzejść: HB 3,5156447; Metropolis 0,4607712; Barker 0,8694091. Przy każdej z tych kinetyk w kanale miękkim zachodzi D=2m_K, a nie arbitralnie wybrane natężenie szumu.
3. **Wspólne prawo krytyczne dopiero po prawidłowym przeskalowaniu.** Lokalny, kompaktowo zatrzymany limit został wyprowadzony dla specjalnej klasy startów uzyskanych przez zaokrąglenie rzeczywistych liczności. Nowy argument czwartego momentu kontroluje różnicę między kinetetyczną współrzędną miękką a tym samym odczytem termodynamicznym.
4. **Wynik negatywny w sprawie unikalności FIN.** Kanoniczne prawdopodobieństwo późniejszego wyjścia górą przedziału ±1 przy lambda=0 jest takie samo dla tych trzech kinetyk, około 0,41685565169838. Samo nie identyfikuje A7 ani reguły kinetycznej.

To jest wyprowadzenie konsekwencji **modelu matematycznego**, nie eksperymentalny test FIN. Przedstawiony dowód jest analityczny i warunkowy na przyjętym wcześniej dokładnym foldzie/certyfikacie; nie został automatycznie zweryfikowany w asystencie dowodowym ani niezależnie zrecenzowany. Liczby z programów są tylko diagnostyką `binary64`.

Kampania kończy się handoffem dla architekta. Nie włączono nowych praw fizycznych, testów laboratoryjnych ani szerokich symulacji.
