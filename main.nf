#!/usr/bin/env nextflow

// nextflow run main.nf -resume

nextflow.enable.dsl=2

params.taxonomy = "data/short_taxonomy.csv"
params.scripts = "scripts"
params.results = 'results'
params.bootstrap = false
params.bootstrap_reps = 1000
params.min_support = 70
params.max_bad_frac = 0.3
params.min_seq_id_mmseqs = 0.3
params.c_mmseqs = 0.8
params.cov_mode_mmseqs = 1
params.max_clusters = 50


workflow {
    // 1. input taxonomy file
    // taxonomy_ch = channel.fromPath(params.taxonomy)

    // 2. select genome assemblies (JNOS summaries)
    // CSV file -> species stream. Parse species one by one
    species_ch = channel
        .fromPath(params.taxonomy)
        .splitText()
        .map { it.trim() }
        .filter { it }

    summaries_ch = select_genomes(species_ch)
    
    // 3. select best assemblies
    assemblies_ch = select_best_assemblies(summaries_ch)

    //4. download proteomes
    proteomes_ch = download_genomes(assemblies_ch)
    
    //5. rCreate gene to species mapping
    gene_maps_ch = create_mapping(proteomes_ch)

    //merge all outputs to single file
    gene_maps_file_ch = gene_maps_ch
        .collectFile(
            name: "gene_to_species.tsv",
            storeDir: "results/mapping",
            newLine:true
        )

    //6. run mmseqs2 to make 'BLAST' and create gene families
    // collect all proteomes into single file 
    merged_faa_ch = proteomes_ch.collectFile(
        name: 'all_proteomes.faa'
    )
    clusters_ch = run_mmseqs(merged_faa_ch)

    //7. filter orthologous genes (1-to-1). returns orthologs1to1.tsv
    orthologs_ch = filter_orthologs(clusters_ch.clusters)

    //8 Extract cluster ids  from othologs (tsv file)
    cluster_ids_ch = orthologs_ch
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def cid = row['cluster_id']
            def safe = cid.replace('|','_').replace('.','_')
            tuple(cid, safe)
        }
        .distinct { it[0] }
        .take(params.max_clusters)



    //9. Make alignemnt (MAFFT) on clusters)
    // orthologs_ch and merged_faa_ch are singleton channels
    fasta_inputs_ch = cluster_ids_ch

    
    fasta_ch = make_clusters_fasta(fasta_inputs_ch)
    aln_ch = mafft_align(fasta_ch)

    //10. Count ML trees on given alignments. One alignment per one script execution
    trees_ch = gene_tree_ml(aln_ch)

    //11. Filter trees only if they were made using bootstrap method. Otherwise return trees_ch (previous output)
    trees_stream_ch = trees_ch.flatten()

    filtered_trees_ch = params.bootstrap \
        ? filter_gene_trees(trees_stream_ch) \
        : trees_stream_ch

    // 12. Trim fasta headers to unambigous names and concatenate alignments
    // tuple (cluster_id, path) -> path)
    trimmed_aln_ch = trim_alignment_headers(aln_ch)

    trimmed_paths_ch = trimmed_aln_ch
        .map { cluster_id, aln -> aln }
        .collect()

    species_order_file_ch = channel.fromPath(params.taxonomy)

    supermatrix_ch = concatenate_alignments(
        trimmed_paths_ch,
        species_order_file_ch
    )


    // 13. Make species_tree  
    species_tree(supermatrix_ch)

    // 14. Consensus tree. Collect all gene trees into one channel
    tree_paths_ch = trees_ch.map { cid, tree -> tree }
    
    filtered_tree_paths_ch = params.bootstrap \
        ? filter_gene_trees(tree_paths_ch) \
        : tree_paths_ch

    consensus_tree_ch = gene_tree_consensus(filtered_tree_paths_ch.collect())

}


process select_genomes {

    tag { species }
    // Copy output to that path
    publishDir "data/proteomes/candidates", mode: 'copy'

    input:
    val species

    output:
    path "${species.replaceAll(' ', '_')}.summary.json"


    script:
    def safe = species.replaceAll(' ', '_')

    """
    ls ./
    bash ${projectDir}/scripts/select_genomes.sh \
        "${species}" \
        "${safe}.summary.json"
    """
}


process select_best_assemblies {
    tag { summary.simpleName }

    publishDir "data/proteomes/assemblies", mode: 'copy'

    input:
    path summary

    output:
    file "${summary.simpleName}.assemblies.tsv"

    script:
    """
    python3  ${projectDir}/scripts/select_best_assemblies.py \
        --input ${summary} \
        --output "${summary.simpleName}.assemblies.tsv"
    """
}


process download_genomes {
    tag { best_assembly_tsv.simpleName }
    
    publishDir "data/proteomes/sequences", mode: 'copy'

    input:
    path best_assembly_tsv
    
    output:
    path "*.faa"

    script:
    """
    python3 ${projectDir}/scripts/download_genomes.py \
        --input ${best_assembly_tsv} \
        --output_zipped zipped \
        --output_sequences .
    """
}

process create_mapping {
    tag { proteome.simpleName }

    publishDir "results/mapping/gene_to_species", mode: 'copy'
    
    input:
    path proteome
    
    output:
    path "${proteome.simpleName}.gene_to_species.tsv"


    script:
    """
    python3 ${projectDir}/scripts/create_mapping_gene_to_specie.py \
        --input ${proteome} \
        --output "${proteome.simpleName}.gene_to_species.tsv"
    """
}


process run_mmseqs {

    publishDir "results/clusters/mmseqs2", mode: 'copy'

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
        clusters \
        tmp \
        --min-seq-id ${params.min_seq_id_mmseqs} \
        -c ${params.c_mmseqs} \
        --cov-mode ${params.cov_mode_mmseqs}
    """
}


process filter_orthologs {
    publishDir "results/clusters", mode: 'copy'
    
    input:
    path clusters

    output:
    file "orthologs1to1.tsv"

    script:
    """
    python3 ${projectDir}/scripts/filter_cluster.py \
        --input ${clusters} \
        --output orthologs1to1.tsv
    """
}

process make_clusters_fasta {

    tag { safe }
    publishDir "results/clusters/fasta", mode: 'copy'
    cpus 4

    input:
    tuple val(cluster_id), val(safe)

    output:
    tuple val(cluster_id), path("${safe}.faa")

    script:
    """
    python3 ${projectDir}/scripts/create_fasta_from_clusters.py \
        --cluster_id '${cluster_id}' \
        --orthologs ${workflow.projectDir}/results/clusters/orthologs1to1.tsv \
        --proteomes ${workflow.projectDir}/results/clusters/mmseqs2/all_proteomes.faa \
        --output ${safe}.faa
    """
}



process mafft_align {
    publishDir "results/alignments", mode: 'copy'
    cpus 4

    input:
    tuple val(cluster_id), path(fasta)

    output:
    tuple val(cluster_id), path("${fasta.simpleName}.aln.faa"), emit: aln

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
//         "results/gene_trees/bootstrap" :
//         "results/gene_trees/no_bootstrap"
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


process gene_tree_ml {
    publishDir {
        params.bootstrap ?
        "results/gene_trees/bootstrap" :
        "results/gene_trees/no_bootstrap"
    }, mode: 'copy'

    cpus 4
    memory '4 GB'
    time '2h'

    input:
    tuple val(cluster_id), path(aln)

    output:
    tuple val(cluster_id), path("${aln.simpleName}.treefile"), emit: tree

    script:
    def bootstrap_flag = params.bootstrap ? "-B ${params.bootstrap_reps}" : ""

    """
    iqtree2 \
        -s ${aln} \
        -m MFP \
        -nt ${task.cpus} \
        ${bootstrap_flag} \
        -pre tree \
        -quiet

    # IQ-TREE always writes: tree.treefile
    mv tree.treefile ${aln.simpleName}.treefile
    """
}


process filter_gene_trees  {
    publishDir "results/filtered_trees", mode:"copy"

    input:
    path tree
    
    output:
    path "${tree.simpleName}.filtered.treefile"

    script:
    """
    python3 ${projectDir}/scripts/reject_weak_trees.py \
        --input ${tree} \
        --output ${tree.simpleName}.filtered.treefile \
        --min_support ${params.min_support} \
        --max_bad_frac ${params.max_bad_frac}

    """
}


process trim_alignment_headers {

    input:
    tuple val(cluster_id), path(aln)

    output:
    tuple val(cluster_id), path("${aln.simpleName}.trimmed.faa")

    script:
    """
    python3 ${projectDir}/scripts/trim_alignment_headers.py \
        --input ${aln} \
        --output ${aln.simpleName}.trimmed.faa
    """
}


process concatenate_alignments {

    publishDir "results/species_tree", mode: 'copy'

    input:
    path trimmed_alignments
    path species_order

    output:
    path "supermatrix.faa"

    script:
    """
    python3 ${projectDir}/scripts/concatenate_alignments.py \
        --inputs ${trimmed_alignments} \
        --species-order ${species_order} \
        --output supermatrix.faa
    """
}




process species_tree {

    publishDir "results/species_tree", mode: 'copy'

    input:
    path supermatrix

    output:
    path "species_tree.treefile"

    script:
    """
    iqtree2 \
        -s ${supermatrix} \
        -st AA \
        -m MFP \
        -nt AUTO \
        -B 1000 \
        -pre species_tree
    """
}


process gene_tree_consensus {

    publishDir "results/species_tree/consensus", mode: "copy"

    input:
    path treefiles

    output:
    path "gene_tree_consensus.treefile"

    script:
    """
    # create list of gene trees
    for t in ${treefiles}; do
        echo \$t
    done > gene_trees.list

    # Majority-rule consensus (50%)
    iqtree2 \
        -t gene_trees.list \
        -con \
        -pre gene_tree_consensus \
        -quiet
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