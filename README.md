Zaimplementuj i przetestuj pipeline filogenetyczny prowadzący do obliczenia zbioru drzew genów i dalej drzewa genomów, począwszy od genomów z podziałem na geny. Można zastosować własne pomysły – do uzgodnienia.

---

## Część I. Proteomy

Wybrać zbiór genomów. Np.:

- min. kilkudziesięciu genomów wirusów najlepiej pokrewnych (np. koronawirusów),
- ~30 bakterii lub archea o rozmiarze liczonym w genach do 5000,
- min. 10 organizmów wyższych; tu można nieco ograniczyć liczbę genów,
- zestaw genomów z różnych królestw.

Pobrać z NCBI całe genomy jako listy sekwencji białkowych (proteomy).  
Warto dobrze zastanowić się jakie genomy wybrać. W projekcie należy wykazać się znajomością własności wybranych genomów.

---

## Część II. Klastrowanie

Dla podanych proteomów wyliczyć plik BLAST z porównaniami sekwencji genowych metodą *każdy z każdym* i poklastrować metodą MCL, MMSeqs lub alternatywną.

---

## Część III. Rodziny genów (przypadek 1-1)

Wybór/wygenerowanie klastrów (rodzin genów) o jednoznacznych nazwach genomów i zastąpienie nazw genów przez nazwy genomów.  
Część klastrów jest pomijana, np. o małej liczności lub z tylko jednym genomem.

---

## Część IV. Multiuliniowienia

Oblicz multiuliniowenie dla każdej rodziny.

---

## Część V. Drzewa rodzin

Obliczenie drzew rodzin z multiuliniowień metodą NJ (lub inną np. ML, MP).

Użycie drzew z programu uliniawiającego jest możliwe, ale niezalecane  
(obniżenie punktacji bazowej o 13 pkt).

---

## Część VI. Drzewo genomów

Oblicz drzewo genomów metodą konsensusową i superdrzewową.

Dla drzew konsensowych potrzebne są klastry ortologiczne, tzn. takie które mają po dokładnie jednej sekwencji z każdego proteomu.  
Jeśli takich klastrów nie ma albo jest ich mało, należy je wygenerować przez przetworzenie klastrów nieortologicznych (zaproponować metodę).

---

## Część VII. Analiza

Analiza biologiczna.  
Porównanie z istniejącymi poglądami nt. drzewa dla tego zbioru (literatura), porównać do taksonomii NCBI.  
Napisać dokument z podsumowaniem użytych metod i wnioskami biologicznymi.

---

## Część V.a (+4 punkty) – opcjonalne

Wyeliminuj słabo wspierane drzewa z użyciem bootstrappingu.  
Sprawdź, czy to daje lepsze wyniki niż metoda bez eliminacji.

---

## Część V.b. Rodziny paralogów (+4 punkty) – opcjonalne

Zastosuj klastry bez usuwania sekwencji, czyli dopuszczamy paralogi.  
Oblicz superdrzewo i porównaj z wynikami dla klastrów ortologicznych.

---

## Raport

Raport powinien mieć formę krótkiego artykułu z podziałem na:

- **Wstęp**
- **Metody**
- **Wyniki**
- **Wnioski**

z bibliografią na końcu.

**Wstęp** powinien zawierać informację o wybranych genomach i znanych hipotezach dot. ich relacji.  
Np. tutaj powinny być odniesienia do powiązanego artykułu z krótkim podsumowaniem, jak autorzy zrekonstruowali drzewo.

**Metody** powinny opisywać zastosowany pipeline z podsumowaniem narzędzi i innych skryptów, opis zastosowanych niestandardowych rozwiązań (np. filtracji, ukorzeniania itp.) oraz opis zasobów komputerowych (jaki komputer, pamięć, czas działania).

Sekcja **Wyniki** powinna zawierać końcową analizę porównawczą:  
przedstawienie drzew wynikowych, porównanie do tych z literatury, taksonomii NCBI oraz timetree.org.

**Wnioski** powinny podsumowywać podejście, zawierać element krytycznej analizy oraz sugerować co zrobić lepiej.

**Uwaga:** plik prezentacji nie jest równoważny raportowi.

---

## Dodatkowe informacje

Wyślij na Moodle wyniki całego pipeline’u, który powinien być w dużym stopniu zautomatyzowany i najlepiej zrównoleglony; można używać narzędzi do pipelinowania np. Snakemake.

Wymagane:

- skrypty (bash, python, R, Makefile, Snakefile, etc.),
- raport z wynikami (format: odt, doc, pdf),
- README – opis jak używać skrypty,
- drzewa genów w jednym pliku (a nie w tysiącach pojedynczych) + plik/pliki z wynikowym drzewem gatunków/consensusu/etc. w formacie Newick,
- jeśli więcej wariantów obliczeń (np. V.a, V.b) – dodać więcej plików drzew,
- paczka na Moodle powinna mieć rozmiar do kilku megabajtów, dlatego **nie załączać**:
  - plików fasta,
  - proteomów,
  - multiuliniowień,
  - pythonowych site-packages itp.,
  - dostępnych programów i repozytoriów,
- w przypadku prezentacji na ostatnich zajęciach – dołączyć w dniu prezentacji (może być po).

---

## Cennik

- maksimum **40 pkt**, w tym:
  - 27 pkt bazowe,
  - 2 × 4 pkt bonusów (opisanych powyżej),
  - 5 pkt – prezentacja na ostatnim wykładzie/labie,
- wariant z prezentacją 15–25 min na ostatnich zajęciach (wymagana dla terminu 0):
  - wysyłka wyników na Moodle maks. na dzień przed prezentacją,
  - plik prezentacji można dosłać w dniu prezentacji,
- wariant bez prezentacji na zajęciach:
  - wysyłka w sesji I lub II,
  - wymagana prezentacja osobista do końca sesji I (termin I) lub do końca sesji II (termin II) po uprzednim przesłaniu projektów i umówieniu się,
- za brak porównania do literatury: **–3 pkt**,
- za brak porównania do timetree.org lub taksonomii NCBI: **–2 pkt**,
- w przypadku większej liczby zainteresowanych w sesji wyznaczony zostanie termin dodatkowy na prezentacje, ale bez bonusu 5 pkt.
