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

    assemblies_ch.view { it.text }

    // //4. download proteomes
    // proteomes_ch = download_genomes(assemblies_tsv)
    
    // //5. rCreate gene to species mapping
    // create_mapping(proteomes_ch.sequences)

    // //6. run mmseqs2 to make 'BLAST' and create gene families
    // clusters_ch = run_mmseqs(proteomes_ch.sequences)

    // //7. filter orthologous genes (1-to-1)
    // orthologs_ch = filter_orthologs(clusters_ch.clusters)

    // //8 Extract cluster ids  from othologs
    // cluster_ids_ch = orthologs_ch
    //     .splitCsv(header: true, sep: '\t')
    //     .map { row -> row.cluster_id }
    //     .distinct()


    // //9. Make alignemnt (MAFFT) on clusters)
    // fasta_inputs_ch = cluster_ids_ch
    //     .combine(orthologs_ch)
    //     .combine(clusters_ch.allproteomes)
    //     .map { a, b, c, d -> tuple(a, b, c, d) }
    
    // fasta_ch = make_clusters_fasta(fasta_inputs_ch)

    // aln_ch = mafft_align(fasta_ch)

    // //10. Count ML trees on given alignments. One alignment per one script execution
    // trees_ch = gene_tree_ml(aln_ch)

    // //11. Filter trees only if they were made using bootstrap method. Otherwise return trees_ch (previous output)
    // trees_stream_ch = trees_ch.flatten()
    // trees_stream_ch.view { it.name }

    // filtered_trees_ch = params.bootstrap \
    //     ? filter_gene_trees(trees_stream_ch) \
    //     : trees_stream_ch


}



process select_genomes {

    tag { species }
    // Copy output to that path
    publishDir "data/proteomes", mode: 'copy'

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

    publishDir "data/proteomes", mode: 'copy'

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
    publishDir "data/proteomes", mode: 'copy'

    input:
    path assemblies_tsv

    
    output:
    path 'zipped', emit: zipped
    path 'sequences', emit: sequences

    script:
    """
    python3 ${projectDir}/scripts/download_genomes.py \
        --input ${assemblies_tsv} \
        --output_zipped zipped \
        --output_sequences sequences
    """
}

process create_mapping {
    publishDir "results/mapping", mode: 'copy'
    
    input:
    path proteome_sequences
    
    output:
    file 'gene_to_species.tsv'

    script:
    """
    python3 ${projectDir}/scripts/create_mapping_gene_to_specie.py \
        --input ${proteome_sequences} \
        --output gene_to_species.tsv
    """
}


process run_mmseqs {

    publishDir "results/clusters/mmseqs2", mode: 'copy'

    input:
    path proteome_sequences

    output:
    path "clusters_cluster.tsv", emit: clusters
    path "clusters_rep_seq.fasta", emit: repseq
    path "all_proteomes.faa", emit: allproteomes

    script:
    """
    mkdir -p tmp

    python3 ${projectDir}/scripts/make_BLAST-like_clusters.py \
        --input ${proteome_sequences} \
        --output clusters \
        --tmp tmp
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
    maxForks 20
    cpus 4

    input:
    tuple val(cluster_id), path(orthologs), path(proteomes)


    output:
    path "${cluster_id.replace('|','_')}.faa"

    script:
    """
    python3 ${projectDir}/scripts/create_fasta_from_clusters.py \
        --cluster_id '${cluster_id}' \
        --orthologs ${orthologs} \
        --proteomes ${proteomes} \
        --output ${cluster_id.replace('|','_')}.faa
    """
}


process mafft_align {
    publishDir "results/alignments", mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta.simpleName}.aln.faa"

    script:
    """
    mafft --thread ${task.cpus} ${fasta} > ${fasta.simpleName}.aln.faa
    """
}

// remember to make two executions - with and without bootstrap
process gene_tree_ml {
    // Adjust path depending on params.bootstrap. Returns folder path
    // ? means or. If params.bootstrap is True, use first value, else second one
    publishDir {
        params.bootstrap ?
        "results/gene_trees/bootstrap" :
        "results/gene_trees/no_bootstrap"
    }, mode: 'copy'

    cpus 4
    memory '4 GB'
    time '2h'

    input: //one alignemnt = one tree
    path aln
    
    output: //remove last extension and add your own A.aln.faa -> A.aln.treefile
    path "${aln.simpleName}.treefile"

    script:
    def bootstrap_flag = params.bootstrap ? "-B ${params.bootstrap_reps}" : ""

    """
    iqtree2 \
        -s ${aln} \
        -m MFP \
        -nt ${task.cpus} \
        ${bootstrap_flag} \
        -pre ${aln.simpleName} \
        -quiet
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


// process  {
//     input:
//     path
    
//     output:
//     path

//     script:
//     """
    
//     """
// }


// #TODO: te od serpula mają błąd w mazwie -> ujednolić jakoś 
// # TODO: usun podwójne podłogi __.
// TODO: zoptymalizuj wszystkie procesy tak by każdy z nich przyjmował JEDEN plik jako input a nie cały folder
// TODO: przeanalizuj każdy skrypt od początku do końca