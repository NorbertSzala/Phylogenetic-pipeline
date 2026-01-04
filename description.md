
# Half-raport

W tym pliku zawieram swoje luźne przemyślenia na temat projektu, tak by pomogło mi w stworzenu docelowego raportu.

---

wskazówki:

- zrób mały zestaw genomów testowych, np w obrębie Aspergillus
- Użyj nextflow
- zrób 5a

---

TODO: dodaj moduł pobierający sekwencje z organizmów
TODO: BLAST/MMseqs2 na tych sekwencjach

---

## Wybór genomów

1.1 Zastanawiam się między wyborem genomów grzybowych a wyborem genomów z różnych królestw (to może prowadzić do trudności w klastrowaniu, w wybranych organizmach może po prostu brakować wspólnych genów)
    Za to grzyby są mi całkiem bliskie. Na pewno klasycznie:


|              Gatunek              | Phylum            | Dikarya/nondikarya | Lifestyle              |
| :-------------------------------: | ----------------- | ------------------ | ---------------------- |
|     Saccharomyces cerevisiae      | Ascomycota        | Dikarya            | nectar_tap_saprotroph  |
|   Schizosaccharomyces japonicus   | Ascomycota        | Dikarya            | nectar_tap_saprotroph  |
|     Schizosaccharomyces pombe     | Ascomycota        | Dikarya            | nectar_tap_saprotroph  |
|         Candida albicans          | Ascomycota        | Dikarya            | nectar_tap_saprotroph  |
|         Serpula lacrymans         | Basidiomycota     | Dikarya            | wood_saprotroph        |
|        Coprinopsis cinerea        | Basidiomycota     | Dikarya            | soil_saprotroph        |
|       Aspergillus nidulans        | Ascomycota        | Dikarya            | unspecified_saprotroph |
|         Aspergillus niger         | Ascomycota        | Dikarya            | unspecified_saprotroph |
|        Gigaspora margarita        | Glomeromycota     | non-dikarya        | arbuscular_mycorrhizal |
|          Gigaspora rosea          | Glomeromycota     | non-dikarya        | arbuscular_mycorrhizal |
|       Zymoseptoria tritici        | Ascomycota        | Dikarya            | plant_pathogen         |
|     Pneumocystes carinii B80      | Ascomycota        | Dikarya            | animal_parasite        |
|       Neolecta irregularis        | Ascomycota        | Dikarya            | plant_pathogen         |
|        Malassezia globosa         | Basidiomycota     | Dikarya            | soil_saprotroph        |
|         Puccinia graminis         | Basidiomycota     | Dikarya            | plant_pathogen         |
|          Ramaria rubella          | Basidiomycota     | Dikarya            | ectomycorrhizal        |
| Melampsora larici-populina 98AG31 | Basidiomycota     | Dikarya            | plant_pathogen         |
|     Choanephora cucurbitarum      | Mucoromycota      | non-dikarya        | litter_saprotroph      |
|          Absidia glauca           | Mucoromycota      | non-dikarya        | soil_saprotroph        |
|         Mucor circinatus          | Mucoromycota      | non-dikarya        | soil_saprotroph        |
|          Mucor ambiguus           | Mucoromycota      | non-dikarya        | soil_saprotroph        |
|       Umbelopsis isabellina       | Mucoromycota      | non-dikarya        | soil_saprotroph        |
|        Umbelopsis vinacea         | Mucoromycota      | non-dikarya        | soil_saprotroph        |
|          Podila humilis           | Mortierellomycota | non-dikarya        | soil_saprotroph        |
|      Mortierella polycephala      | Mortierellomycota | non-dikarya        | soil_saprotroph        |
|      Mortierella antarctica       | Mortierellomycota | non-dikarya        | soil_saprotroph        |
|      Mortierella hygrophila       | Mortierellomycota | non-dikarya        | soil_saprotroph        |
|       Modicella reniformis        | Mortierellomycota | non-dikarya        | soil_saprotroph        |
|      Synchytrium microbalum       | Chytridiomycota   | non-dikarya        | plant_pathogen         |
|     Synchytrium endobioticum      | Chytridiomycota   | non-dikarya        | plant_pathogen         |
|        Nosema apis BRL 01         | Microsporidia     | non-dikarya        | animal_parasite        |
|          Nosema ceranae           | Microsporidia     | non-dikarya        | animal_parasite        |
|         Nosema granulosis         | Microsporidia     | non-dikarya        | animal_parasite        |

Do tego chciałbym dołączyć:  Busco N50, L50, genome_size, sum_protein_length, mean_prot_length


Wybrałem łącznie genomów, z czego to *dikarya* a reszta to *non-dikarya*. Ten podział w taksonomii grzybów jest bardzo znaczący. Do *Dikarya* wchodzą *Basidiomycota* oraz *Ascomycota*, a skolei *non-dikarya* to cała reszta.



## Blast i klastrownia

Zastanawiam się nad wyborem MMseqs2 zamiast Blast - jedno narzędzie do alignemntu (znacznie szybsze niż blast) oraz do klastrowania

## Filtracja klastrów - odrzucenie tych pojedynczych

Własny skrypt pythonowy
Przykladowy output klastrowania wraz z interpretacja:

| cluster_id                | specie      | gene_id                   |
| :------------------------ | :---------- | :------------------------ |
| Saccharomyces__cerevisiae | NP_009557.1 | Saccharomyces__cerevisiae | Saccharomyces__cerevisiae|NP_009557.1
| Saccharomyces__cerevisiae | NP_009557.1 | Aspergillus__niger        | Aspergillus__niger|XP_001394994.1
| Saccharomyces__cerevisiae | NP_009557.1 | Malassezia__globosa       | Malassezia__globosa|XP_001731973.1

- Mamy klaster `Saccharomyces__cerevisiae | NP_009557.1` (*cluster_id*), którego reprezentantem jest sekwencja genu **NP_009557.1** pochodząca z *Saccharomyces cerevisiae* (`specie`). Do tego klastra należą również inne geny (`gene_id`) z podanymi informacjami o organiźmie oraz o identyfikatorze genu.
- Razem te geny tworzą homologi,- mają wspólnego przodka i zostały rozdzielone na zasadzie specjacji
- Zostały razem zgrupowane na podstawie podobieństwa sekwencji
- 

## MSA

Zapewne MAFFT

## Drzewa genomów

Fastqtree lub iq-tree (co przystępniejsze) + bootstrap 5A.
Zastanów się nad 5B

## Konsensus

Do rozważenia

## Porównanie z NCBI taxonomy

Robinson folds method, jaccard similarity i coś jeszcze numerycznego
