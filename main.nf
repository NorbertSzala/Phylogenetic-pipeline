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
    
    species_ch.view { "SPECIES='${it}'" }



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
    proteomes_with_qc_ch.view { "DEBUG: $it" }

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
        output:     orthologs_relaxed.tsv
    * STRICT orthologs (true 1-to-1)
    */
    strict = orthology_pipeline(
        clusters_ch.clusters,
        merged_faa_ch,
        0,
        "strict"
    )

    /*
    * RELAXED orthologs (N–k, pseudo-single-copy)
    */
    relaxed = orthology_pipeline(
        clusters_ch.clusters,
        merged_faa_ch,
        params.max_missing,
        "relaxed")


    /*
    * STRICT → consensus (IQ-TREE)
    */
    strict_species_tree_ch = strict.trees
        .map { cid, tree -> tree }
        | strip_gene_ids_from_trees
        | collectFile(name: "gene_trees.strict.list", newLine: true)
        | consensus_tree



    /*
    * 18. Build supertree (ASTRAL) from RELAXED gene trees
    */
    astral_input_ch = relaxed.trees
        .map { cid, tree -> tree }
        .collectFile(
            name: 'all_gene_trees.relaxed.tre',
            newLine: true
        )

    astral_clean_ch = prepare_astral_trees(astral_input_ch)
    astral_tree_ch  = astral_tree(astral_clean_ch)


    // Count statistics after run
    relaxed_tree_list_ch = relaxed.trees
    .map { cid, tree -> tree }
    .collectFile(name: "gene_trees.relaxed.tre", newLine:true)

    tree_stats(relaxed_tree_list_ch)

}

// subworkflow
workflow orthology_pipeline {

    take:
        clusters_ch
        merged_faa_ch
        max_missing
        label

    main:

        /*
         * Filter orthologs with given missing-data tolerance
         */
        filt = label == "strict"
            ? filter_orthologs_strict(clusters_ch)
            : filter_orthologs_relaxed(clusters_ch)

        orthologs_ch = filt.orthologs
        stats_ch     = filt.stats



        /*
        * 8. Extract cluster IDs
            input:      orthologs_relaxed.tsv
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
            orthologs_relaxed.tsv
            all_proteomes.faa
        output:
            tuple(cluster_id, cluster_sequences.faa)
        */
        fasta_ch = label == "strict"
            ? cluster_ids_ch.combine(orthologs_ch).combine(merged_faa_ch) | make_clusters_fasta_strict
            : cluster_ids_ch.combine(orthologs_ch).combine(merged_faa_ch) | make_clusters_fasta_relaxed


        /* 10. Make alignemnt (MAFFT) on clusters)
            input:      fassta_ch - one FASTA file with sequences for a single cluster
            output:     alignments
        */
        aln_ch = label == "strict"
            ? mafft_align_strict(fasta_ch)
            : mafft_align_relaxed(fasta_ch)


        /*11. Count ML trees on given alignments. 
            input:      tuple val(cluster_id), path(aln)
            ouptut:     tuple val(cluster_id, treefiles)    
        */
        trees_ch = label == "strict"
            ? gene_tree_ml_strict(aln_ch)
            : gene_tree_ml_relaxed(aln_ch)


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
    tuple val(species), path("*gene_to_species.tsv")


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

process make_clusters_fasta_relaxed {

    tag "${safe}.relaxed"
    publishDir "${params.results}/clusters/fasta/relaxed", mode: 'copy'
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

process mafft_align_relaxed {
    publishDir "${params.results}/alignments/relaxed", mode: 'copy'
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


// remember to make two executions - with and without bootstrap
// process gene_tree_ml {
//     // Adjust path depending on params.bootstrap. Returns folder path
//     // ? means or. If params.bootstrap is True, use first value, else second one
//     publishDir {
//         params.bootstrap ?
//         "${params.results}/gene_trees/bootstrap" :
//         "${params.results}/gene_trees/no_bootstrap"
//     }, mode: 'copy'

//     cpus 4
//     memory '4 GB'
//     time '2h'

//     input: //one alignemnt = one tree
//     path aln
    
//     output: //remove last extension and add your own A.aln.faa -> A.aln.treefile
//     path "${aln.simpleName}.treefile"

//     script:
//     def bootstrap_flag = params.bootstrap ? "-B ${params.bootstrap_reps}" : ""

//     """
//     iqtree2 \
//         -s ${aln} \
//         -m MFP \
//         -nt ${task.cpus} \
//         ${bootstrap_flag} \
//         -pre ${aln.simpleName} \
//         -quiet
//     """
// }



process gene_tree_ml_strict {

    tag "${cluster_id}"

    publishDir "${params.results}/gene_trees/strict", mode: 'copy'

    cpus 4
    memory '4 GB'
    time '2h'

    /*
     * KLUCZOWE:
     * - nie zabijaj pipeline'u na pojedynczym błędzie
     * - brak retry (błąd jest deterministyczny)
     */
    errorStrategy 'ignore'
    maxRetries 0

    input:
        tuple val(cluster_id), path(aln)

    output:
        /*
         * emit tylko wtedy, gdy treefile faktycznie powstał
         * (Nextflow sam odfiltruje brakujące outputy)
         */
        tuple val(cluster_id), path("${aln.simpleName}.treefile")

    script:
    def bootstrap_flag = params.bootstrap
        ? "-B ${params.bootstrap_reps}"
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


process gene_tree_ml_relaxed {

    tag "${cluster_id}"

    publishDir "${params.results}/gene_trees/relaxed", mode: 'copy'
    cpus 4
    memory '4 GB'
    time '2h'

    /*
     * KLUCZOWE:
     * - nie zabijaj pipeline'u na pojedynczym błędzie
     * - brak retry (błąd jest deterministyczny)
     */
    errorStrategy 'ignore'
    maxRetries 0

    input:
        tuple val(cluster_id), path(aln)

    output:
        /*
         * emit tylko wtedy, gdy treefile faktycznie powstał
         * (Nextflow sam odfiltruje brakujące outputy)
         */
        tuple val(cluster_id), path("${aln.simpleName}.treefile")

    script:
    def bootstrap_flag = params.bootstrap
        ? "-B ${params.bootstrap_reps}"
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

    output:
        path "orthologs.tsv", emit: orthologs
        path "filter.stats.tsv", emit: stats

    script:
    """
    python3 ${projectDir}/scripts/filter_cluster.py \
        --input ${clusters} \
        --n_genomes ${params.n_genomes} \
        --max_missing 0 \
        --output orthologs.tsv
    """
}



process filter_orthologs_relaxed {
    publishDir "${params.results}/stats/filtering/relaxed", mode:'copy'

    input:
        path clusters

    output:
        path "orthologs.tsv", emit: orthologs
        path "filter.stats.tsv", emit: stats

    script:
    """
    python3 ${projectDir}/scripts/filter_cluster.py \
        --input ${clusters} \
        --n_genomes ${params.n_genomes} \
        --max_missing ${params.max_missing} \
        --output orthologs.tsv \
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
        path "stripped.treefile"
    
    script:
    """
    sed -E 's/\\|[^,:)]*//g' ${treefile} > stripped.treefile
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
        path gene_tree_list

    output:
        path "consensus.treefile"

    script:
    """
    iqtree2 \
        -t ${gene_tree_list} \
        -con \
        -pre consensus \
        -quiet

    # IQ-TREE writes consensus as *.contree
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
     *   - relaxed orthology (after paralog filtering)
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

// process  {
//     input:
//     path
    
//     output:
//     path

//     script:
//     """
    
//     """
// }


// TODO: Konsensus drzew genowych
// TODO: Astral - supertree
// TODO: opracować metodę gdy będzie mało klastrów 1-1 (poczekać)