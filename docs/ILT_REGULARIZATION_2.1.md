# DOSY Lab 2.1 · Automatyczny dobór λ i stabilność ILT

## Jak używać

W oknie **Analizuj ILT**, **Global ILT** albo **Mapa DOSY 2D → Adaptive Local ILT** otwórz **Auto λ / stabilność**. Program przygotuje skan w tle. Możesz zmienić granice, liczbę prób od 7 do 41 i kryterium. Domyślnie skan obejmuje trzy dekady poniżej i powyżej bieżącego λ, w granicach od 10⁻¹² do 100. W AILT środek domyślnego zakresu to wspólne ręczne λ. Zaznacz **Potwierdzam wiarygodne σ(I) z pomiaru** tylko przy kompletnych, niezależnie oszacowanych błędach standardowych całek. Opcja **Sprawdź siatkę D** sprawdza szerszy zakres i inną liczbę komórek.

Wykresy pokazują krzywą L lub warunkową, reszty, normę kary, pozycje maksimów i podgląd rozkładu. Kliknięcie punktu pokazuje zapisaną próbę, liczby i masy przedziałów. Biała kreska oznacza bieżące λ, zielona sugestię, a żółta podgląd. Na samej krzywej zaznaczona jest najbliższa obliczona próba. Okno umożliwia przesuwanie, prostokątny zoom, pełny ekran i eksport PNG. **Zastosuj rekomendację** zmienia λ. **Przywróć poprzednie λ** przywraca wartość sprzed ostatniego zastosowania. Ręczna regulacja pozostaje dostępna.

W AILT skan dotyczy poszczególnych zgodnych grup. Zastosowanie zapisuje tylko wyniki z dostępną sugestią. Grupy bez sugestii i niezależnie liczone ROI korzystają ze wspólnego ręcznego λ. Wybór **Ręczna λ** wraca do jednego λ dla grup. Grupowanie, tolerancja D, balans amplitud oraz polityka mapy są osobnymi ustawieniami. Zmiana danych lub warunków grupowania unieważnia wcześniejsze rekomendacje grupowe.

## Funkcje celu i kryteria

W indywidualnym ILT każdy zanik jest dzielony przez największą bezwzględną aktywną całkę. Wektor wag z aktywnych punktów ma normę 1. Przy włączonym ważeniu proporcje wag wynoszą 1/σ. Dla komórek siatki o szerokości 1/P w znormalizowanej współrzędnej D albo log D solver NNLS minimalizuje ważoną kwadratową resztę plus λ razy kwadrat normy kary. Przy rzędzie zerowym kara wynosi √P razy wektor mas komórek. Dla drugiej pochodnej jest to √P razy (0,05 P)² razy druga różnica mas komórek.

Global ILT i grupowy AILT używają wspólnej siatki, skali oddzielnej dla każdego ROI oraz kary rzędu zerowego. Funkcja celu to suma kwadratów ważonych reszt, λ P razy kwadrat normy Frobeniusa F i koszt wspólnej rzadkości μ √M razy suma norm wierszy F. M to liczba ROI. AILT może dodatkowo zmniejszyć wagę kolumny na podstawie balansu amplitud. Solver grupowy kończy po osiągnięciu zadanej względnej normy kroku proksymalnego albo limitu iteracji. Próby bez zbieżności nie są kandydatami. Wartości λ z indywidualnej kary krzywizny i z grupowej kary zerowego rzędu nie są wprost porównywalne.

**Krzywa L.** Krzywiznę krzywej log–log wyznaczam różnicami na osi log λ. Odrzucam krańce skanu, płaskie odcinki, brak wymiany między resztą i karą oraz konkurencyjne maksima krzywizny o wysokości co najmniej 75% najwyższego maksimum. Dwa konkurencyjne załamania dają zakres zamiast wartości. Odrzucam też rozwiązanie z ważoną resztą co najmniej 75% reszty najlepszego modelu stałego. Są to reguły przesiewowe, a zwrócone załamanie jest punktem startowym do oceny. Dla Global ILT i AILT μ jest stałe. Na osi poziomej krzywej warunkowej używam pierwiastka z sumy kwadratu ważonej reszty i kosztu μ. Na osi pionowej pozostaje norma kary λ. To przekrój zadania z dwiema karami przy ustalonym μ. [Hansen i O’Leary, 1993](https://doi.org/10.1137/0914086)

**Zasada rozbieżności.** Dla błędów standardowych σ w oryginalnych jednostkach, znormalizowanych wag W i skali całek s oczekiwany kwadrat normy szumu to suma składników (W σ / s)² dla aktywnych obserwacji. Kwantyle 2,5% i 97,5% to przybliżenie rozkładem chi-kwadrat dopasowanym do dwóch momentów sumy kwadratów przy niezależnych błędach gaussowskich. To przedział odniesienia dla szumu, a nie przedział ufności dla λ lub D. Do porównania służy sama ważona reszta danych, bez kosztu μ. Wybieram tylko próbę w paśmie szumu przy jednoznacznym przejściu przez oczekiwaną normę. Spośród dwóch sąsiednich prób preferuję mniejszy udział masy maksimów niestabilnych, a następnie bliższą normę. Brak przejścia, brak kompletu σ lub ostrzeżenie QC o systematycznym wpływie kształtu, fazy lub bazy ograniczają wnioskowanie. Brak dokładnego przejścia przy zaufanych σ nie uruchamia automatycznego zastępstwa krzywą L. Osobna wskazówka z krzywej pozostaje widoczna. [Hansen, zasada rozbieżności w dyskretnych problemach odwrotnych](https://www2.compute.dtu.dk/~pcha/AIRtoolsII/Tutorial/DublinDay1.pdf)

**Stabilność.** Łączę maksima między sąsiednimi próbami wyłącznie przy wzajemnie jednoznacznym dopasowaniu w odległości log D najwyżej 0,12 lub 1,5 maksymalnego kroku siatki log D, zależnie od tego, która granica jest większa. Złączenia i rozszczepienia są nierozstrzygnięte. Stabilna ścieżka występuje w co najmniej 80% prób i ma rozpiętość log D najwyżej 0,12. Umiarkowanie wrażliwa występuje w co najmniej połowie prób przy rozpiętości do 0,30. Pozostałe trwałe ścieżki są silnie zależne od λ, a pojedyncze lub brzegowe pozostają nierozstrzygnięte. Udziały liczę przez sumowanie mas komórek, a wartości całek składowych podaję w oryginalnej skali dla b = 0. Udział zbiorczy jest sumą mas po normalizacji każdej kolumny. Podgląd grupowy rysuje znormalizowaną normę L₂ każdego wiersza, zgodnie z definicją rozkładu zbiorczego w Global ILT. Przy silnie tłumionych składowych I(0) jest ekstrapolacją. Dwie kontrole siatki stosują to samo λ i nie są ponowną optymalizacją.

**Walidacja predykcyjna.** Benchmark ocenia predykcję na oddzielnej, gęstej siatce analitycznej z wartościami bez szumu. Walidacja krzyżowa na pomiarze pozostaje poza automatycznym doborem. Podział wymagałby grupowania powtórzonych b, co najmniej sześciu różnych b w zbiorze uczącym, odpowiedniej reprezentacji małych i dużych b oraz ponownego wyznaczania normalizacji i wag wyłącznie w zbiorze uczącym. Bez tego błąd predykcji byłby mylący dla typowego krótkiego DOSY.

## Wyniki

W benchmarku syntetycznym porównano 25 wartości od 10⁻¹⁰ do 1 z analitycznym zanikiem bez szumu. Najmniejszy błąd pełnego skanu jest wynikiem odniesienia dostępnym tylko dlatego, że znana jest prawda symulacyjna. Maksima porównano w log D z tolerancją 0,20, bez wielokrotnego przypisywania. Wiersze poniżej pokazują osobne realizacje.

| Przypadek | Sugerowane λ | RMS predykcji / najlepsze RMS w skanie | Fałszywe / pominięte maksima |
| --- | ---: | ---: | ---: |
| Monoeksponenta | zakres 1,0 × 10⁻¹⁰–8,25 × 10⁻⁸ | — / 0,00080 | — |
| Dwie rozdzielone | zakres 1,0 × 10⁻¹⁰–1,47 × 10⁻⁶ | — / 0,00102 | — |
| Dwie bliskie | zakres 1,0 × 10⁻¹⁰–2,15 × 10⁻⁷ | — / 0,00051 | — |
| Słaba druga składowa | 1,47 × 10⁻⁶ | 0,00287 / 0,00146 | 0 / 1 |
| Niski S/N | 0,147 | 0,02980 / 0,01923 | 1 / 2 |
| Rozkład ciągły | 1,21 × 10⁻³ | 0,00370 / 0,00231 | nie dotyczy |
| Słabe pokrycie b | zakres 1,0 × 10⁻¹⁰–8,25 × 10⁻³ | — / 0,01467 | — |
| Szum heteroscedastyczny | 5,62 × 10⁻⁷ | 0,00223 / 0,00104 | 0 / 0 |
| Global ILT, nakładanie i różne amplitudy ROI | 3,83 × 10⁻⁶ | 0,00387 / 0,00260 | 0 / 0 |

W 20 niezależnych realizacjach szumu na przypadek otrzymałem sugestię w 14 przypadkach mono, 13 rozdzielonych, 19 bliskich i 4 przy niskim S/N. Dla rozdzielonych składowych żadna z 13 rekomendacji nie dodała fałszywego maksimum ani nie pominęła składowej w użytej tolerancji. Dla bliskich składowych w każdej z 19 rekomendacji brakowało jednego maksimum. Przy niskim S/N cztery rekomendacje dały łącznie cztery fałszywe i osiem pominiętych maksimów. Podwojenie deklarowanego σ w heteroscedastycznym przykładzie przesunęło wybór z 5,62 × 10⁻⁷ do 3,16 × 10⁻³, a zmniejszenie σ o połowę usunęło jednoznaczną sugestię.

Na zapisanym projekcie strychniny skan AILT objął jedną grupę 18 ROI, 13 prób od 10⁻⁹ do 10⁻³, μ = 10⁻⁴. Krzywa warunkowa sugeruje **λ = 10⁻⁵** jako punkt startowy, z sąsiednim zakresem **3,16 × 10⁻⁶–3,16 × 10⁻⁵**. Brak niezależnych σ oznacza, że zasada rozbieżności nie ocenia tego wyniku. Przy λ = 10⁻⁶ ważona norma reszty wynosi 0,02009, a przy sugestii 0,02818. Maksimum zbiorcze jest w obu próbach w tej samej komórce, 2,148 × 10⁻⁹ m²/s. Wynik skanu zależy od danych i ustawień, w tym μ i siatki. [Wyniki maszynowe](../benchmarks/ilt_regularization_results.json)

Przeszło 107 testów, w tym 69 testów obliczeń i 38 testów GUI. Przy tym samym ręcznym λ benchmark map 2.1 nadal zwraca identyczne wejścia, mono D i numeryczne rozwiązania ILT jak przed zmianą. Sprawdziłem w Qt podgląd, zastosowanie i cofnięcie λ, eksport oraz przełączenie AILT na wspólne ręczne λ. Global ILT zachował μ podczas zastosowania. Zbieżność mówi o rozwiązaniu funkcji celu. Reszty mówią o zgodności z pomiarem. Stabilność opisuje reakcję na badane parametry. Żadna z tych cech sama nie potwierdza fizycznego istnienia osobnej populacji.

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
& '.venv/Scripts/python.exe' benchmarks/ilt_regularization.py 'projects/Strychnina_CDCl3_bpp_ste_led_diffusion-1-1.proj'
```
