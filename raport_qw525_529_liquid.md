# Raport Badań QW-525 do QW-529: Liquid Crystal Dynamics

**Data:** 2025-12-04
**Cel:** Weryfikacja hipotezy, że Jądro $K(d)$ jest Ciekłym Kryształem (Nadciekłym), a cząstki są defektami topologicznymi.

---

## Wyniki Symulacji

### **QW-525: Stabilność Wiru (Vortex Stability)**
*   **Cel:** Czy wir (defekt topologiczny) jest stabilny w polu Jądra?
*   **Wynik:** 🔴 **PORAŻKA (Niestabilność)**.
    *   Początkowy wir ($m=1$) uległ inwersji lub rozpadowi ($m_{final} = -1.0$).
    *   **Wniosek:** W obecnej fazie ("Szkło") proste wiry są niestabilne. Struktura medium wymusza ich rekonfigurację.

### **QW-526: Nadciekłość (Superfluidity)**
*   **Cel:** Czy przepływ informacji jest bezoporowy?
*   **Wynik:** 🟢 **SUKCES**.
    *   Zachowanie przepływu (Overlap) jest bardzo wysokie.
    *   **Wniosek:** Medium zachowuje się jak **Nadciecz**. Informacja propaguje się bez strat (rozpraszania), mimo "sztywnej" struktury Jądra. Potwierdza to wniosek z `gemini_sum.md`.

### **QW-527: Uporządkowanie Nematyczne (Order)**
*   **Cel:** Czy spiny/fazy porządkują się w domeny?
*   **Wynik:** 🟡 **SZKŁO (Glassy)**.
    *   Funkcja korelacji oscyluje i zanika. Brak długozasięgowego porządku (ferro/nematic).
    *   **Wniosek:** Oscylacje Jądra ($K(d)$ zmienia znak) powodują **frustrację**. System nie może się zdecydować na jeden kierunek, tworząc stan "Szkła Spinowego" (Spin Glass).

### **QW-528: Topnienie (Melting)**
*   **Wynik:** Parametr porządku pozostaje niski ($\sim 0.1$) niezależnie od temperatury.
*   **Wniosek:** System jest "zamrożony" w stanie nieuporządkowanym. Brak wyraźnego przejścia fazowego w badanym zakresie.

### **QW-529: Stałe Elastyczności (Frank Energy)**
*   **Cel:** Czy medium stawia opór deformacjom (skręceniu)?
*   **Wynik:** 🟢 **SUKCES (Sztywność)**.
    *   Stała elastyczności skręcania $K_{twist} \approx 5.37 > 0$.
    *   **Wniosek:** Medium jest **Elastyczne**. Reaguje na deformacje jak ciało stałe lub ciekły kryształ, a nie jak zwykła ciecz (gdzie $K=0$).

---

## Podsumowanie Generalne: "Nadciekłe Szkło"
Badania QW-525+ ujawniły niezwykłą naturę Jądra przy zamrożonych parametrach:
1.  Jest **Nadciekłe** (przepuszcza informację bez oporu).
2.  Jest **Elastyczne** (ma sztywność skręcania).
3.  Jest **Szkliste/Sfrustrowane** (brak uporządkowania dalekiego zasięgu, niestabilność wirów).

**Werdykt:**
To nie jest zwykły Ciekły Kryształ (uporządkowany), ale **Nadciekłe Szkło (Superfluid Glass)**.
*   **Dobra wiadomość:** Wyjaśnia propagację fal (światła) i brak oporu próżni.
*   **Zła wiadomość:** Frustracja (szklistość) uniemożliwia formowanie się stabilnych, prostych cząstek (wirów). Cząstki w takim medium muszą być bardziej złożonymi strukturami (np. splotami, które "oszukują" frustrację), a nie prostymi wirami.

**Rekomendacja:** Aby uzyskać stabilne cząstki, musimy znaleźć sposób na "lokalne roztopienie" szkła lub znalezienie konfiguracji topologicznych odpornych na frustrację (np. Skyrmiony w szkle spinowym).
