# Instrukcja użytkownika DOSY Lab

[Powrót do README](../README.md)

1. [Otwarcie danych](#1-otwarcie-danych)
2. [Nawigacja i processing](#2-nawigacja-i-processing)
3. [Korekcja fazy](#3-korekcja-fazy)
4. [Linia bazowa](#4-linia-bazowa)
5. [Całki](#5-całki)
6. [Wartości b](#6-wartości-b)
7. [Dopasowanie i kontrola jakości](#7-dopasowanie-i-kontrola-jakości)
8. [Statystyka D](#8-statystyka-d)
9. [Projekty i eksport](#9-projekty-i-eksport)
10. [Najczęstsze problemy](#10-najczęstsze-problemy)
11. [Mapa DOSY 2D](#11-mapa-dosy-2d)

## 1. Otwarcie danych

1. Uruchom aplikację zgodnie z [README](../README.md#pierwsze-uruchomienie-na-windows).
2. Kliknij **Otwórz DOSY** i wybierz surowy plik JEOL `.jdf`. Program przeliczy widma z zespolonych FID i wyświetli warstwę z najniższym dostępnym **b**. Gdy **b** jest nieznane, warstwę referencyjną wybiera na podstawie gradientu, jeśli ten jest dostępny.
3. Otwórz **Info**, aby sprawdzić nazwę pliku, parametry akwizycji i aktualnego processingu, pochodzenie **b** oraz metadane zapisane w JDF.
4. Pole **Warstwa** wybiera pojedyncze widmo. **Pokaż referencję** wraca do warstwy referencyjnej. **Nałóż warstwy** przyjmuje numery od 1, oddzielone przecinkami, np. `1,8,16`. Pierwsza wskazana warstwa wyznacza wartość całki pokazywaną przy sygnale i w tabeli. **Pokaż Im** dołącza do podglądu część urojoną.

Po zmianie warstwy podglądu wszystkie ROI nadal obejmują te same zakresy ppm w każdej warstwie. W tabeli widzisz całkę dla aktualnego podglądu, a model zaniku korzysta z całej serii aktywnych punktów.

## 2. Nawigacja i processing

| Gest lub przycisk | Działanie |
| --- | --- |
| Rolka nad widmem | Zmienia skalę intensywności w pionie. |
| **Zoom** i przeciągnięcie lewym przyciskiem | Powiększa zaznaczony zakres ppm. Zaznaczenie jest widoczne w trakcie przeciągania. |
| Przeciągnięcie prawym przyciskiem nad widmem głównym | Powiększa zakres ppm także podczas pracy z narzędziami całkowania. |
| **Przesuwaj** i przeciągnięcie | Przesuwa widok w pionie i poziomie. |
| Dwuklik lub **Cały zakres** | Przywraca pełny zakres widma. |
| **Pełny ekran** | Otwiera powiększony widok. **Esc** wraca do głównego okna. |
| **Zapisz PNG** | Zapisuje bieżący widok 1D w jasnej kolorystyce z granatowym widmem, zachowując przybliżenie. |
| **Widok** | Ustawia kolor i grubość linii widma oraz widoczność siatki. |

Gdy zatrzymasz kursor na przycisku głównego okna, pojawi się krótka podpowiedź.

**Processing** pokazuje ustawienia przetwarzania. Kolejność operacji to apodyzacja FID, dopełnienie zerami, FFT, wspólna korekcja fazy, opcjonalne odjęcie bazy, a następnie całkowanie rzeczywistej części widma. Zmiany przeliczają wszystkie warstwy na podstawie oryginalnych FID.

| Opcja | Znaczenie |
| --- | --- |
| **LB [Hz]** | Wykładnicza apodyzacja FID przez czynnik `exp(−π · LB · t)`. |
| **Gauss: FWHM [Hz]** | Gaussowska apodyzacja FID przez czynnik `exp(−(π · GB · t)² / (4 ln 2))`. Oba okna można łączyć. Zero wyłącza daną apodyzację. |
| **Zero filling** | Dopełnienie zerami do 1×, 2×, 4× lub 8× liczby punktów. Zagęszcza siatkę FFT. Rozdzielczość fizyczną wyznacza akwizycja. |
| **PH0, PH1, Pivot** | Wspólna faza opisana w następnym rozdziale. |
| **Referencja: przesunięcie ppm** | Dodaje wspólne przesunięcie osi chemicznej. Po zmianie sprawdź położenie zapisanych granic ROI. |
| **Opóźnienie filtra [próbki]** | Stosuje oddzielną rampę fazową do korekcji opóźnienia filtra cyfrowego. Zero wyłącza tę korekcję. |
| **Przytnij pasmo filtra według JDF** | Ogranicza pasmo do wartości zapisanej w JDF. Program uwzględnia liczbę próbek i czas akwizycji, gdy plik zgłasza przycięcie, a metadane wskazują na pełny FID. |
| **Warstwa referencyjna** | Ustala punkt startowy podglądu fazowania i warstwę wskazywaną przyciskiem **Pokaż referencję**. |
| **Zapisz profil / Wczytaj profil** | Przenosi ustawienia processingu jako profil JSON. Po wczytaniu na innym widmie sprawdź zwłaszcza fazę, referencję i zakresy bazy. |

## 3. Korekcja fazy

![Interaktywna korekcja fazy z oznaczeniem pivotu i podglądem przed korekcją](images/korekcja-fazy.png)

*Podgląd wspólnej fazy w syntetycznej serii gradientowej.*

1. Kliknij **Korekcja fazy**. Okno pokazuje początkowo warstwę referencyjną.
2. Kliknij widmo w pobliżu wybranego piku, żeby ustawić **pivot**, albo wpisz jego pozycję w ppm. To punkt, w którym korekcja ręczna ma wartość PH0.
3. Przytrzymaj **lewy przycisk** i przesuwaj mysz góra–dół, aby zmieniać **PH0**. Przytrzymaj **prawy przycisk**, aby tak samo zmieniać **PH1**. **Shift** włącza dokładniejsze strojenie. Czułość obu gestów można ustawić w oknie.
4. Skorzystaj z **Zoom zakresu** i **Przesuń widmo**, jeśli potrzebujesz zbliżenia lub linia bazowa wyszła poza widok. Wyłącz wybrany tryb, aby wrócić do gestów fazowania. Rolka zmienia skalę intensywności. Opcje **Pokaż Im** oraz **Pokaż fazę sprzed edycji** pomagają ocenić efekt.
5. Sprawdź w polu **Podgląd** jedną z warstw o większym **b**. **Zastosuj wspólną fazę do wszystkich warstw** zatwierdza jeden zestaw parametrów dla całej serii. **Anuluj** porzuca edycję w tym oknie.

PH1 opisuje różnicę fazy w stopniach na szerokość całego akwizycyjnego pasma SW, w kierunku rosnących ppm. Ręczna faza w punkcie `ppm` to `PH0 + PH1 · (ppm − pivot) / SW_ppm`. Zoom nie zmienia tej zależności. Wyświetlane ppm zwykle maleją od lewej do prawej. Wartość liczbowa PH1 zapisana w innym programie może wymagać przeliczenia, jeśli używa on innego pivotu, kierunku lub znaku.

Jeśli metadane JDF są odpowiednie, DOSY Lab osobno szacuje opóźnienie filtra na podstawie zapisanych etapów decymacji. **Info** pozwala sprawdzić, czy model został zastosowany. Okno fazowania pozwala też ręcznie dodać rampę `360° × opóźnienie w próbkach` do PH1. Przy korzystaniu z tej funkcji sprawdź, czy niezależna korekcja filtra nie obejmuje już tej samej rampy.

## 4. Linia bazowa

**Zaznacz bazę** pozwala wskazać obszary, w których sygnał powinien odpowiadać lokalnemu tłu. Wybierz kilka przedziałów bez pików, najlepiej po obu stronach sygnałów całkowanych. Niebieskie pasma oznaczają dane użyte do dopasowania wielomianu, a żółta linia pokazuje podgląd bazy.

W **Processing** domyślny rząd to 1, czyli linia prosta. Możesz ustawić rząd od 1 do 5. Domyślnie zaznaczony **Osobny wielomian dla każdej warstwy** pokazuje podgląd baz dopasowanych na tych samych zakresach ppm. Przed zatwierdzeniem możesz przełączyć podgląd na bazę referencyjną. **Zatwierdź odejmowanie bazy** zawsze odejmuje od każdej warstwy jej własną bazę przed całkowaniem. Dotyczy to także wczytanych projektów z aktywną korekcją. Upewnij się, że niebieskie pasma nie zachodzą na ROI. Poza skrajnymi wskazanymi obszarami wielomian jest ekstrapolowany, więc warto obejrzeć tam jego podgląd.

Aby skasować jeden przedział, zaznacz go w liście panelu **Processing** i kliknij **Usuń obszar**. **Cofnij korekcję** wyłącza odejmowanie i zachowuje maskę. **Wyczyść bazę** nad widmem usuwa wszystkie obszary i wyłącza odejmowanie.

## 5. Całki

ROI to zakres ppm, na którym program całkuje dany sygnał w każdej warstwie. Całki mają znak i powstają z części rzeczywistej po wspólnej korekcji fazy oraz ewentualnym odjęciu bazy. Jednostką jest intensywność w jednostkach umownych razy ppm.

- **Integracja myszą** aktywuje tworzenie zakresów. Kliknij jedną granicę, przytrzymaj przycisk, przeciągnij i puść na drugiej granicy. Kolorowe pasmo pokazuje zakres na żywo.
- **Dodaj całkę ręcznie** przyjmuje dwie granice w ppm. Dla liczb z polskim przecinkiem wpisz np. `3,60; 3,68`.
- **Modyfikuj całkę** pozwala przeciągnąć środek ROI, by przesunąć cały zakres, lub krawędź, by zmienić tylko granicę.
- **Podziel całkę** wymaga przeciągnięcia jednej krawędzi do środka ROI. Po puszczeniu w miejscu podziału powstają dwa sąsiadujące zakresy.
- **Usuń całkę** usuwa pozycję wybraną w tabeli. Każda całka ma własne ID. Dwuklik w kolumnie **Grupa / związek** pozwala wpisać nową etykietę lub wybrać wcześniej użytą z rozwijanej listy. Wspólna etykieta oznacza świadomą decyzję o zestawieniu tych sygnałów w statystyce.

Identyfikator, np. **S1**, znajduje się nad sygnałem, a wartość całki poniżej. W kolumnie **Całka** można nadać wybranemu sygnałowi wartość odniesienia, np. **1,00**. Program przeliczy współczynnik skali i proporcjonalnie odświeży pozostałe całki na widmie oraz w tabeli. Zaokrąglenie do dwóch miejsc dotyczy prezentacji. Dopasowania zachowują oryginalne całki.

W syntetycznym przykładzie możesz utworzyć dwie grupy i sprawdzić, jak działa porównanie.

| Grupa | Zakresy ppm | Zadane D [m²/s] |
| --- | --- | --- |
| A | 3,60–3,68 · 3,74–3,82 · 5,26–5,33 | około 3,8 × 10⁻¹⁰ |
| B | 4,66–4,74 | około 8,0 × 10⁻¹⁰ |

## 6. Wartości b

**Lista b** przypisuje każdej warstwie wartość **b** w **s/m²**. Kolumna **g**, jeśli jest dostępna, pokazuje gradient w **T/m**. Dwie warstwy z takim samym **b** nadal stanowią dwa odrębne punkty pomiarowe.

1. Sprawdź, czy program rozpoznał gotową oś **b** lub zweryfikował model fabryczny. Dla protonowych BPP-STE-LED (w tym obsługiwanego wariantu z presaturacją) i BPP-DSTE-LED rozpoznaje przebieg kodowania w osadzonym JXP i porównuje zapisany `st_coef` z modelowymi wartościami, kiedy współczynnik jest obecny. Nazwa JXP jest widoczna w **Info**. Obsługiwane nazwy kształtów w tym modelu to `SQUARE`, `SINE`, `SQUARE_SINE`, `TRAPEZOID` i `S_RECTANGLE`. Każdy ma osobne współczynniki.
2. Dla sekwencji spoza modelu podaj własne **b**, wyznaczone dla faktycznej sekwencji. Możesz użyć **Import CSV** z jedną kolumną albo pliku z kolumną `b`, `b_s_m2` lub `b_s_per_m2`. Poniżej minimalny przykład formatu. Liczby służą wyłącznie do pokazania zapisu.
3. Wpisz dokładnie jedną skończoną, nieujemną wartość na każdą warstwę, w kolejności akwizycji. Zachowaj powtórzenia. Kliknij **Zastosuj jako ręczną listę** i opisz źródło wartości w polu uzasadnienia.

~~~csv
b_s_m2
5.70e7
2.99e8
5.41e8
~~~

**b** obliczone z nominalnych amplitud i timingów ma charakter nominalny. Aplikacja sprawdza spójność zapisu i modelu. Rzeczywistą kalibrację pola gradientowego ustala się niezależnie.

## 7. Dopasowanie i kontrola jakości

Dla każdego ROI DOSY Lab dopasowuje model **I(b) = I₀ exp(−bD)** do całek ze znakiem. Wartości **D** i SE(D) raportuje w m²/s. Tabela główna wyświetla obie wielkości w skali 10⁻¹⁰ m²/s. **R²** pomaga ocenić dopasowanie, a **QC** sygnalizuje m.in. problem z liczbą punktów lub niestabilnością wyniku. Wymagane są co najmniej trzy różne aktywne wartości **b**.

- **Zaniki i reszty** pokazują dane, krzywą oraz różnicę pomiar − model. **Podgląd I / pierwsza całka** normalizuje prezentację. Dopasowanie nadal używa całek w ich pierwotnej skali.
- **Całki / punkty fitu** pokazuje poszczególne warstwy i ich udział w dopasowaniu. Punkty możesz wyłączyć ręcznie, pozostają one w projekcie i eksporcie. Ważenie fitu według niepewności całek stosuj, gdy masz wiarygodne wartości σ(I) dla aktywnych warstw.
- **Widma w aktywnym ROI** pomaga porównać kształt sygnału w wybranych warstwach. Jeśli zakres obejmuje nakładające się składniki, dobre **R²** samo w sobie nie potwierdza, że **D** opisuje jedną substancję.
- Na wykresie zaniku lub reszt przeciągnij prostokąt, by powiększyć wycinek. Dwuklik przywraca pełny zakres. Dedykowane przyciski otwierają wykresy na pełnym ekranie.

## 8. Statystyka D

![Okno Statystyka D z tabelą sygnałów, wykresem zbiorczym i diagnostyką](images/statystyka-d.png)

*Dane syntetyczne. Kliknięcie punktu D pokazuje oryginalne punkty zaniku, dopasowanie i reszty dla odpowiadającego mu ROI.*

1. Nadaj sygnałom etykiety w kolumnie **Grupa / związek** w głównym oknie. Otwórz **Statystyka D** i wybierz grupę w polu **Związek**. Tabela pokazuje sygnały w eksperymencie, a podsumowanie korzysta z sygnałów z wybranej grupy.
2. Użyj kolumny **Uwzględnij**, żeby włączać i wyłączać ROI w statystyce. Zmiana dotyczy zestawienia. Dopasowanie pojedynczego ROI i jego punkty pozostają w projekcie.
3. Kliknij punkt **D** lub odpowiadający mu wiersz. Panel **Zanik wybranego ROI** pokazuje zapisany model, **D**, SE(D), **R²**, korekcję bazy oraz rzeczywiste punkty I(b). Włącz **Pokaż reszty**, jeśli chcesz ocenić strukturę błędu. Po wskazaniu punktu zaniku dostępne wartości obejmują **g**, **b**, całkę, przewidywanie, resztę i S/N.
4. Użyj **Zachowaj kopię do porównania** lub **Wczytaj zapis statystyki**, aby zestawić osobne eksperymenty. W zakładce porównania każdy pomiar ma własne podsumowanie.
5. **Eksport CSV + zapis statystyki** zapisuje tabele i migawkę JSON. Przyciski **PNG** i **SVG** eksportują wykres wybrany w polu **Eksportowany wykres** z białym tłem, czarnymi opisami i granatowymi danymi. **Wyłącz siatkę** działa osobno dla wykresu D oraz pary wykresów zaniku i reszt.

Wykres **D** pokazuje wyłącznie zaznaczone sygnały wybranego związku. Punkty mają słupki ±SE(D), a wykres zawiera średnią i pasmo ±1 SD. **Linie ±2 SD** można włączyć dodatkowo. Gdy liczba uwzględnionych sygnałów na to pozwala, podsumowanie zawiera średnią, medianę, SD próby, CV, opisowe SEM, minimum, maksimum i **N**. Dla N < 2 brakuje podstaw do wyznaczenia SD próby, CV i SEM.

- **SE(D)** opisuje niepewność jednego dopasowania wynikającą z jego modelu i danych.
- **SD** opisuje rozrzut wyników **D** pomiędzy wybranymi sygnałami. **SEM** to opisowe SD/√N. Sygnały z jednego eksperymentu mogą mieć skorelowane błędy, np. wspólne odchylenie kalibracji gradientu, więc to SEM ma ograniczoną interpretację jako niepewność średniej.
- Opcjonalna **Średnia ważona wg SE(D)** używa wag 1/SE². Włącza się ją po potwierdzeniu poprawności i porównywalności SE. Potrzebne jest dodatnie i skończone SE dla każdego uwzględnionego ROI. Domyślnym wynikiem jest średnia arytmetyczna.

## 9. Projekty i eksport

**Zapisz projekt** zachowuje ustawienia processingu, wspólną fazę, ręczne **b**, granice ROI, grupy, maski punktów oraz kalibrację całek. Projekt JSON przechowuje odwołanie do źródłowego JDF i jego SHA-256. **Wczytaj projekt…** odtwarza analizę z tym źródłem. Przy przenoszeniu projektu zachowaj oryginalny JDF i sprawdź ścieżkę zapisaną w JSON. Profil processingu i migawka statystyki pełnią odrębne funkcje.

**Eksport Statystyki D** w głównym oknie zapisuje pakiet wyników analizy. Poniżej znajduje się lista najważniejszych plików.

| Plik | Zawartość |
| --- | --- |
| `roi_results.csv` | Granice ROI, grupy, **D**, SE, **R²** i QC. |
| `integrals_and_fits.csv` | Całki według warstwy, **g**, **b**, przewidywania, reszty i udział w fitach. |
| `compound_statistics.csv` | Podsumowania według jawnie zdefiniowanych grup. |
| `DOSY_results.xlsx` | Arkusze ROI, punktów, statystyki, processingu i eksperymentu. |
| `experiment.json`, `processing.json`, `fit_results.json` | Metadane, opcje i szczegóły dopasowań. |
| `project.json`, `statistics_snapshot.json` | Projekt i migawka wyników do porównań. |
| `complex_data.npz` | FID, ppm, widma zespolone, baza oraz wartości gradientów i **b**. |
| `embedded_sequence.txt` | Tekst sekwencji, jeśli występuje w JDF. |
| `Raport.md` i pliki PNG/SVG | Podsumowanie i wykresy. |

Pakiet eksportu może zawierać pełne dane pomiarowe i treść osadzonej sekwencji. Przed udostępnieniem wybierz konkretne pliki odpowiednie dla odbiorcy.

## 10. Najczęstsze problemy

| Objaw | Co sprawdzić |
| --- | --- |
| Widmo otwiera się bez wartości D | Otwórz **Listę b** i sprawdź źródło wartości, a w tabeli punktów liczbę różnych aktywnych **b**. Gdy **b** jest nieznane, wykres nadal pokazuje całki względem gradientu albo numeru warstwy. |
| Całka jest ujemna lub silnie zależy od granic | Sprawdź fazę absorpcyjną, bazę i możliwe nakładanie się pików. Porównaj kilka warstw tym samym processingiem. |
| Linia bazowa faluje już przed korekcją wielomianową | Zacznij od metadanych i korekcji filtra cyfrowego w **Info** oraz **Processing**. Następnie sprawdź wspólne PH0 i PH1 na kilku pikach. |
| Dopasowanie ma dobre R², ale **D** różni się między sygnałami | Sprawdź reszty, widma w ROI, kalibrację **b**, granice całek i przypisanie sygnałów do związków. |
| Statystyka pokazuje N = 0 | Sprawdź grupę **Związek**, pole **Uwzględnij** i QC sygnałów. |
| Projekt nie znajduje JDF | Sprawdź ścieżkę do oryginalnego pliku w projekcie i jego dostępność. |
| Skrypt startowy zamyka się lub brakuje bibliotek | Spróbuj [instalacji ręcznej](../README.md#pierwsze-uruchomienie-na-windows) i uruchom `app.py` z tego samego środowiska. |

## 11. Mapa DOSY 2D

![Mapa DOSY 2D z połączonym widmem 1D, listą całek i prowadnicami](images/mapa-dosy-2d.png)

*Dane syntetyczne. Cztery całki grupy A tworzą pasmo w pobliżu 3,9 × 10⁻¹⁰ m²/s. Prowadnice wskazują sygnał grupy B przy 4,70 ppm.*

1. Wyznacz całki i poczekaj na zakończenie dopasowań. Nad widmem, obok **Wyczyść bazę**, kliknij **Mapa DOSY 2D**.
2. Z prawej strony wybierz, które całki mają się znaleźć na mapie. Dostępne są całki z poprawnie zakończonym dopasowaniem i dodatnim, skończonym **D**. Wybór dotyczy mapy i nie zmienia projektu ani statystyki.
3. Przeciągnij prostokąt w trybie **Zoom**. Zmiana zakresu ppm na mapie od razu zmienia zakres widma 1D nad nią. **Przesuwaj** włącza przesuwanie, a **Cały zakres** lub dwuklik przywraca widok wszystkich wyświetlanych sygnałów.
4. Włącz **Punkt / prowadnice** i kliknij mapę. Wartość **D** pojawi się przy przecięciu linii i zostanie zapisana w PNG. Pionowa prowadnica pojawia się też na widmie 1D. Odczyt dotyczy miejsca kliknięcia, bez przyciągania do najbliższego piku.
5. Rolka nad **1D** zmienia wysokość widma. Rolka nad **2D** zmienia **Podłogę**, czyli próg widoczności w procentach najsilniejszego sygnału. Aby uwidocznić słabsze sygnały, obniż podłogę i przesuń suwak **Wzmocnienie 2D**. Suwak rozjaśnia punkty powyżej progu, a najsilniejsze mogą osiągnąć maksymalną jasność.
6. W menu **Widok** możesz osobno ustawić kolor mapy 2D i widma 1D. Przycisk **Siatka** włącza lub wyłącza siatkę obu wykresów.
7. **Zapisz widok PNG** zachowuje aktualny zakres, intensywność, siatkę i prowadnice. Plik ma białe tło, czarne osie i opisy oraz granatowe widmo i sygnały. Eksport renderuje tę wersję w tle, z zachowaniem kolorystyki roboczej w aplikacji.

Mapa umieszcza dodatnią intensywność sfazowanego widma referencyjnego w granicach każdej całki przy współczynniku **D** z jej zapisanego dopasowania. Domyślny kolor to pomarańczowy. Ujemna intensywność widma, np. fragment dyspersyjnej linii, nadal jest widoczna w widmie 1D. Wzmocnienie, podłoga i kolor mapy wpływają wyłącznie na podgląd. Jeśli wybrane całki się nakładają, ich wspólne punkty są pomijane. Wyłączenie jednej z nich przywraca wspólny fragment w drugiej.

**Szerokość pasma D** steruje graficznym rozmyciem w pionie jako FWHM w procentach największego widocznego **D**. Położenie pasma wynika z fitu. Każda całka ma jedno przypisane **D**, więc mapa nie rozdziela związków, których sygnały leżą w tym samym ROI.
