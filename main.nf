#!/usr/bin/env nextflow

// nextflow run main.nf -resume

nextflow.enable.dsl=2

workflow {
    // 1. input taxonomy file
    // taxonomy_ch = channel.fromPath(params.taxonomy)

    /*
     2. Select genome assemblies (JNOS summaries)
            input : species
            output: tuple(species, summaries.tsv)
     */
    species_ch = Channel
        .fromPath(params.taxonomy)
        .splitText()
        .map { it.trim().replaceAll(/\r/, '') }
        .filter { it }
    
    n_genomes = species_ch.count().first()





    summaries_ch = select_genomes(species_ch)
    
    /*
      3. Select best assemblies
            input : tuple(species, summaries.tsv)
            output: tuple(species, assemblies.tsv)
     */
    assemblies_ch = select_best_assemblies(summaries_ch)

    /*
     4. Download proteomes
            input : tuple(species, assemblies.tsv)
            output: tuple(species, proteome.faa, qc.tsv)
     */
    download_genomes(assemblies_ch)

    proteomes_ch = download_genomes.out.faa   // (species, faa)
    qc_ch        = download_genomes.out.qc    // (species, qc)

    proteomes_with_qc_ch =
        proteomes_ch
            .map { species, faa -> tuple(species, faa) }
            .join(
                qc_ch.map { species, qc -> tuple(species, qc) }
            )

    /*
    5. Create gene to spciecies mapping
        input:      tuple(species, proteome.faa, qc.tsv)
        output:     gene_to_species.tsv
    */
    gene_maps_ch = create_mapping(proteomes_with_qc_ch)
    //merge all outputs to single file
    gene_maps_file_ch = gene_maps_ch
        .collectFile(
            name: "gene_to_species.tsv",
            storeDir: "${params.results}/mapping",
            newLine:true
        )

    /*
    6. run mmseqs2 to make 'BLAST' and create gene families
        collect all proteomes into single file 

        input:tuple (species, proteome.faa, qc.tsv)
        output: tuple(clusters_cluster.tsv, clusters_rep_seq.fasta)
    */
    merged_faa_ch = proteomes_ch
        .map { species, proteome -> proteome }
        .collectFile(name: 'all_proteomes.faa', sort: true)


    clusters_ch = run_mmseqs(merged_faa_ch)



    /*
    7. Filter orthologous genes (1 to 1)

        input:      clusters from mmseqs2
        output:     orthologs_strict.tsv
    * STRICT orthologs (true 1-to-1)
    */
    strict = orthology_pipeline(
        clusters_ch.clusters,
        merged_faa_ch,
        "strict",
        n_genomes
    )



    /*
    * STRICT → consensus (IQ-TREE)
    * 1) strip gene IDs -> taxa == species
    * 2) count taxa
    * 3) filter only complete trees
    * 4) concat into one .tre
    */
    strict_stripped_trees = strict.trees \
        | map { cid, tree -> tree } \
        | strip_gene_ids_from_trees \
        | map { cid, tree -> tree }

    strict_complete_trees = strict_stripped_trees \
        | map { tree -> tuple(tree.simpleName, tree) } \
        | count_taxa_in_tree \
        | map { cid, tree, ntaxa_file -> tuple(cid, tree, ntaxa_file.text.trim().toInteger()) } \
        | filter { cid, tree, ntaxa -> ntaxa > 5 } \
        | map { cid, tree, ntaxa -> tree }

    consensus_input = strict_complete_trees.collect() | concat_trees
    consensus_tree(consensus_input)


    /*
    * 18. Build supertree (ASTRAL) from STRICT gene trees
    */
    astral_input_ch = strict.trees
        .map { cid, tree -> tree }
        .collectFile(
            name: 'all_gene_trees.strict.tre',
            newLine: true
        )

    astral_clean_ch = prepare_astral_trees(astral_input_ch)
    astral_tree_ch  = astral_tree(astral_clean_ch)



}

// subworkflow
workflow orthology_pipeline {

    take:
        clusters_ch
        merged_faa_ch
        label
        n_genomes

    main:

        /*
         * Filter orthologs with given missing-data tolerance
         */

        filt = filter_orthologs_strict(clusters_ch, n_genomes)


        orthologs_ch = filt.orthologs
        stats_ch     = filt.stats



        /*
        * 8. Extract cluster IDs
            input:      orthologs_strict.tsv
            output:     Channel< tuple(cluster_id, safe_cluster_id) >
        */
        cluster_ids_ch = orthologs_ch
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def cid = row.cluster_id
                def safe = cid.replace('|','_').replace('.','_')
                tuple(cid, "${safe}.${label}")
            }
            .distinct { it[0] }
            .take(params.max_clusters)


        /* 9. Create per-cluster FASTA files
        input:
            tuple(cluster_id, safe_cluster_id)
            orthologs_strict.tsv
            all_proteomes.faa
        output:
            tuple(cluster_id, cluster_sequences.faa)
        */
        fasta_ch = cluster_ids_ch.combine(orthologs_ch).combine(merged_faa_ch) | make_clusters_fasta_strict


        /* 10. Make alignemnt (MAFFT) on clusters)
            input:      fassta_ch - one FASTA file with sequences for a single cluster
            output:     alignments
        */
        aln_ch = mafft_align_strict(fasta_ch)


        /*11. Count ML trees on given alignments. 
            input:      tuple val(cluster_id), path(aln)
            ouptut:     tuple val(cluster_id, treefiles)    
        */
        trees_ch = gene_tree_ml_strict(aln_ch)


    emit:
        trees = trees_ch
        alignments = aln_ch
}


process select_genomes {
    publishDir "results/selected_genomes"

    tag "${species}"


    input:
    val species

    output:
        tuple val(species), path("*.summary.json")

    script:
    """
    bash ${projectDir}/scripts/select_genomes.sh \
        "${species}" \
        "${species.replaceAll(' ', '_')}.summary.json"
    """
}



process select_best_assemblies {
    publishDir "results/selected_assemblies"
    tag "${species}"

    input:
        tuple val(species), path(summary_json)

    output:
        tuple val(species), path("*.assemblies.tsv")

    script:
    """
    python3 ${projectDir}/scripts/select_best_assemblies.py \
        --input ${summary_json} \
        --output ${species.replaceAll(' ', '_')}.assemblies.tsv
    """
}

process concat_trees {

    input:
        path trees

    output:
        path "trees_for_consensus.tre"

    script:
    """
    cat ${trees} > trees_for_consensus.tre
    """
}


process count_taxa_in_tree {

    input:
        tuple val(cluster_id), path(tree)

    output:
        tuple val(cluster_id), path(tree), path("ntaxa.txt")

    script:
    """
    ntaxa=\$(grep -o ',' ${tree} | wc -l)
    ntaxa=\$((ntaxa + 1))
    echo \$ntaxa > ntaxa.txt
    """
}




process download_genomes {

    input:
        tuple val(species), path(assemblies_tsv)

    output:
        tuple val(species), path("qc_summary.tsv"), emit: qc
        tuple val(species), path("*.faa"), emit: faa

    script:
    """
    python3 ${projectDir}/scripts/download_genomes.py \
        --input ${assemblies_tsv} \
        --outdir .
    """
}


process create_mapping {
    publishDir "results/created_mapping"
    tag "${species}"

    input:
    tuple val(species), path(proteome), path(qc_tsv)
    
    output:
    tuple val(species), path("*.gene_to_species.tsv")


    script:
    """
    python3 ${projectDir}/scripts/create_mapping_gene_to_specie.py \
        --input ${proteome} \
        --output "${species.replaceAll(' ', '_')}.gene_to_species.tsv"
    """
}


process run_mmseqs {
    publishDir "${params.results}/clusters/mmseqs2", mode: 'copy'

    input:
    path merged_faa

    output:
    path "clusters_cluster.tsv", emit: clusters
    path "clusters_rep_seq.fasta", emit: repseq

    script:
    """
    mkdir -p tmp

    mmseqs easy-cluster \
        ${merged_faa} \
        mmseqs_out \
        tmp \
        --min-seq-id ${params.min_seq_id_mmseqs} \
        -c ${params.c_mmseqs} \
        --cov-mode ${params.cov_mode_mmseqs}
    
    mv mmseqs_out_cluster.tsv clusters_cluster.tsv
    mv mmseqs_out_rep_seq.fasta clusters_rep_seq.fasta
    """
}



process make_clusters_fasta_strict {

    tag "${safe}.strict"
    publishDir "${params.results}/clusters/fasta/strict", mode: 'copy'
    cpus 4

    input:
        tuple val(cluster_id), val(safe), path(orthologs_tsv), path(all_proteomes)

    output:
        tuple val(cluster_id), path("${safe}.faa")

    script:
    """
    python3 ${projectDir}/scripts/create_fasta_from_clusters.py \
        --cluster_id '${cluster_id}' \
        --orthologs ${orthologs_tsv} \
        --proteomes ${all_proteomes} \
        --pick ${params.pick_representative} \
        --output ${safe}.faa

    """
}



process mafft_align_strict {
    publishDir "${params.results}/alignments/strict", mode: 'copy'
    cpus 4

    input:
        tuple val(cluster_id), path(fasta)

    output:
        tuple val(cluster_id), path("${fasta.simpleName}.aln.faa")

    script:
    """
    mafft --thread ${task.cpus} ${fasta} > ${fasta.simpleName}.aln.faa
    """
}


process gene_tree_ml_strict {

    tag "${cluster_id}"

    publishDir "${params.results}/gene_trees/strict", mode: 'copy'

    cpus 4
    memory '4 GB'
    time '2h'

    errorStrategy 'ignore'
    maxRetries 0

    input:
        tuple val(cluster_id), path(aln)

    output:
        tuple val(cluster_id), path("${aln.simpleName}.treefile")

    script:
    def bootstrap_flag = params.bootstrap
        ? "-b ${params.bootstrap_reps}"
        : ""


    """
    # --- minimalna walidacja alignmentu ---
    nseq=\$(grep -c "^>" ${aln})

    if [ "\$nseq" -lt 3 ]; then
        echo "[SKIP] ${aln} has only \$nseq sequences" >&2
        exit 0
    fi

    # --- budowa drzewa ---
    iqtree2 \
        -s ${aln} \
        -m MFP \
        -nt ${task.cpus} \
        ${bootstrap_flag} \
        -pre tree \
        -quiet

    # --- sanity check ---
    if [ ! -f tree.treefile ]; then
        echo "[SKIP] IQ-TREE failed for ${aln}" >&2
        exit 0
    fi

    mv tree.treefile ${aln.simpleName}.treefile
    """
}



process filter_orthologs_strict {
    publishDir "${params.results}/stats/filtering/strict", mode:'copy'

    input:
        path clusters
        val n_genomes

    output:
        path "orthologs.tsv", emit: orthologs
        path "filter.stats.tsv", emit: stats

    script:
    """
    python3 ${projectDir}/scripts/filter_cluster.py \
        --input ${clusters} \
        --n_genomes ${n_genomes} \
        --output orthologs.tsv
    """
}




process tree_stats {

    publishDir "${params.results}/stats/tree_stats", mode:'copy'

    input:
        path gene_trees

    output:
        path "tree_stats.tsv"

    script:
    """
    awk '
    {
      gsub(/[();]/,"");
      n=split(\$0,a,",");
      print NR "\t" n
    }' ${gene_trees} > tree_stats.tsv
    """
}


process strip_gene_ids_from_trees {

    input:
        path treefile
    
    output:
        tuple val(treefile.simpleName), path("${treefile.simpleName}.stripped.tre")
    
    script:
    """
    sed -E 's/\\|[^,:)]*//g' ${treefile} > ${treefile.simpleName}.stripped.tre
    """
}

process consensus_tree {

    /*
     * Build a consensus species tree from multiple gene trees.
     *
     * This step implements a majority-rule consensus (50%) over
     * individual gene trees inferred from orthologous clusters.
     *
     * Input :
     *   - gene_tree_list : text file with one gene tree path per line
     *
     * Output:
     *   - consensus.treefile : consensus species tree
     *
     * Notes:
     *   - Works best when gene trees contain the same set of taxa
     *     (e.g. strict 1-to-1 orthologs).
     */

    publishDir "${params.results}/species_tree/consensus", mode:'copy'


    input:
        path treefile

    output:
        path "consensus.treefile"

    script:
    """
    iqtree2 -con ${treefile} -minsup 0.5 -pre consensus -quiet
    mv consensus.contree consensus.treefile
    """
}

process prepare_astral_trees {

    input:
        path gene_trees

    output:
        path "astral_input.tre"

    script:
    """
    # 1) remove empty lines
    # 2) strip gene IDs (|NP_xxx)
    grep -v '^\\s*\$' ${gene_trees} \
      | sed -E 's/\\|[^,:)]*//g' \
      > astral_input.tre
    """
}

process astral_tree {

    /*
     * Build a species tree using a supertree / summary method (ASTRAL).
     *
     * ASTRAL infers a species tree under the multispecies coalescent
     * model using a collection of gene trees, allowing for:
     *   - missing taxa
     *   - gene tree discordance
     *
     * Input :
     *   - all_gene_trees : single file containing all gene trees (Newick)
     *
     * Output:
     *   - astral.treefile : species tree inferred by ASTRAL
     *
     * Notes:
     *   - Particularly useful when strict 1-to-1 orthologs are scarce.
     *   - Complements the consensus-based species tree.
     */
    publishDir "${params.results}/species_tree/astral", mode:'copy'

    input:
        path all_gene_trees

    output:
        path "astral.treefile"

    script:
    """
    java -jar ${projectDir}/Astral/astral.5.7.8.jar \
            -i ${all_gene_trees} \
            -o astral.treefile
    """
}