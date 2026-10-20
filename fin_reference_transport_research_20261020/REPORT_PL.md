# Nowy wynik: spójne prawo przy zmianie opisu stanu

Data: 2026-10-20. Badanie tej samej warunkowej hipotezy
FIN-HYP-HIST-ME-001, bez nowego postulatu i bez przyjęcia jej jako fizyki.

## Co zostało wyprowadzone

Nie wykonywano kolejnej recenzji własnego wyniku ani szukania większego N,
które „naprawiłoby” HB. Pytanie było nowe: **kiedy różne opisy tej samej
historii dają dokładnie to samo wybrane prawo?**

Powstały cztery obiekty:

1. **Pełny test równoważności miar odniesienia.** Dwie dodatnie miary na
   tym samym grafie dają to samo prawo krokowe wtedy i tylko wtedy, gdy
   różnią się czynnikami przypisanymi wierzchołkom. Przy dodatnich pętlach
   wystarczy sprawdzić dokładne ilorazy `Kxy^2/(Kxx Kyy)` na krawędziach.
   To test prawa krokowego, nie automatycznie niezmiennik zegara fizycznego.
2. **Wyprowadzona miara dla opisu zredukowanego.** Należy przenosić
   rzeczywiste masy mikroskopowych włókien i użyć `V^T K V`. Warunki
   niezmienniczości i zgodność miary stanów płacą wtedy całe prawo historii,
   a nie tylko średnie szybkości. Sama statyczna kompresja nie wystarcza.
3. **Pełne kryterium zamknięcia w tej klasie.** Po kanonicznym zbilansowaniu
   zadanej miary odniesienia można sprawdzić, czy zredukowana historia jest
   Markowowska. Zamknięcie dla jednej dodatniej zgodnej miary stanów oznacza
   w tej konkretnej rodzinie zamknięcie dla wszystkich takich miar.
   Nie trzeba dobierać ich pod pożądany wynik. Normalizacja bazowa jest
   narzędziem dowodu, nie nowym fizycznym przygotowaniem układu.
4. **Dokładne składanie redukcji.** Po dwóch etapach trzeba zachować masy
   z pierwszego; zastąpienie ich ponownie jedynkami może zmienić prawo,
   a nie tylko zapis macierzy. Podano jawny wymierny kontrprzykład.

## Konkretny obiekt dla FIN

Z liczby rzeczywiście dozwolonych zmian jednej oznaczonej kopii otrzymano
dla dowolnego skończonego N miarę odniesienia procesu liczności:

    Kcount(n,n-e_a+e_b)=sqrt(n_a(n_b+1)),
    Kcount(n,n)=1.

To czynniki normalizacji orbit, **nie** szybkości HB, nawias martyngału,
operatory kreacji kwantowej ani fizyczne cząstki. Wykorzystuje się prawdziwe
wagi Gibbsa z czynnikiem wielomianowym N!/prod n_a!, a nie ich zamiennik.
Dowód dla wszystkich skończonych N jest kombinatoryczny; testy nie udają
enumeracji dużych krat liczności.

Dla q12/N2 sprawdzono dwie drogi:

    144 stany oznaczone -> 12 różnic -> 7 klas inwersji,
    144 stany oznaczone -> 78 typów liczności -> te same 7 klas.

Z poprawnymi wagami obie dają to samo prawo. W klasach dwuelementowych
pojawia się masa własna 3: zawiera też prawdziwe zmiany niewidoczne w tym
opisie, nie tylko puste odświeżenia. Nie wolno ich zgubić.

## Wyniki przeciwne zbyt szerokim wnioskom

Wyprowadzono dokładne przykłady, w których kompresja tworzy pamięć,
spłaszczanie mas zmienia klasę odniesienia, a zbyt proste kryterium
normalizacji diagonalnej zawodzi przy wewnętrznych przejściach.
Pełne kryterium zbilansowane zachowuje te rozróżnienia.

Niezgodność pierwotnej hipotezy z HB **nie znika po usunięciu informacji
o orientacji różnicy**. W wymiernym przykładzie q12 pierwsze przejście
między wskazanymi klasami ma prawdopodobieństwo 1/42 zamiast 19/764,
z różnicą -17/16044. To przykład algebraiczny, nie liczby rzeczywistego A7.
Właściwy wniosek A7 jest jakościowy i warunkowy względem odziedziczonej
dodatniej, parzystej i niejednorodnej miary p.

## Weryfikacja i znaczenie dla F01

Przechodzi **55 ograniczonych kontroli** w rachunku wymiernym i Q(sqrt(5)).
Są to kontrole skończone: grafy 4, 50 i 144 stanów, 78 typów q12/N2,
normalizacje i kontrprzykłady. Nie są to 55 dowodów ogólnych ani niezależna
recenzja. Zakres szukania wcześniejszych odpowiedników jest zapisany osobno;
nie przypisuje się tej pracy globalnego pierwszeństwa matematycznego.

Nowy postęp polega na **wyprowadzeniu brakującej reguły przenoszenia
i składania opisu**, gdy dane mikroskopowe są już zadane. Nadal nie
wyprowadzono, dlaczego natura miałaby wybrać właśnie te dane lub zasadę
entropijną. Strict F01/F03/F06 pozostają otwarte; AGENTS, źródła i poprzednie
plany zachowano bez zmian. Cel „badaj teorię” pozostaje aktywny.

Pełny dowód: [PROOF_OR_OPEN_LEMMA.md](PROOF_OR_OPEN_LEMMA.md).
