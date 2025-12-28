(project_gp) norbert@norbert:~/Studia_VisualCode/Genomika_Porównawcza/projekt/results/clusters/mmseqs2$ head -n 2 *

==> all_proteomes.faa <==
>Saccharomyces__cerevisiae|NP_001018029.1 Coq21p [Saccharomyces cerevisiae S288C]
MRNELYQLWCVASAARGVAKSSFVRANSAMCEYVRTSNVLSRWTRDRQWEAAKALSQRVKKEYAAN

File with merged sequences from all organisms. Input for mmseqs2



---------------------------------------------------------------
==> clusters_all_seqs.fasta <==
>Saccharomyces__cerevisiae|NP_001018029.1
>Saccharomyces__cerevisiae|NP_001018029.1 Coq21p [Saccharomyces cerevisiae S288C]

unused.


---------------------------------------------------------------
==> clusters_cluster.tsv <==
Saccharomyces__cerevisiae|NP_001018029.1	Saccharomyces__cerevisiae|NP_001018029.1
Saccharomyces__cerevisiae|NP_001018031.2	Saccharomyces__cerevisiae|NP_001018031.2

The most important file. Key to 3rd step - gene families
format: <member_sequence_id>    <representative_sequence_id>
Sequence <representative_sequence_id> belongs to <member_sequence_id>
In this example, those clusters contain only one sequnce each (singletons)


---------------------------------------------------------------
==> clusters_rep_seq.fasta <==
>Saccharomyces__cerevisiae|NP_011256.1 Vel1p [Saccharomyces cerevisiae S288C] 
MSFLSIFTFFSVLISVATTVRFDLTNVTCKGLHGPHCGTYVMEVVGQNGTFLGQSTFVGADVLTESAGDAWARYLGQETRFLPKLTTIASNETKNFSPLIFTTNINTCNPQSIGDAMVPFANTVTGEIEYNSWADTADNASFITGLANQLFNSTDYGVQVASCYPNFASVILSTPAVNIFGKDDTLPDYCTAIQLKAVCPPEAGFD

Representative sequences for clusters
one representative sequence for one cluster