#!/usr/bin/env nextflow

// nextflow run main.nf -resume

nextflow.enable.dsl=2

params.taxonomy = "data/short_taxonomy.csv"
params.scripts = "scripts"
params.results = 'results'

workflow {
    // 1. input taxonomy file
    taxonomy_ch = channel.fromPath(params.taxonomy)

    // //2. select genome assemblies (JNOS summaries)
    script1 = channel.fromPath("scripts/select_genomes.sh")
    json_summaries_dir_ch = select_genomes(taxonomy_ch, script1)
    
    // // //3. select best assemblies
    script2 = channel.fromPath("scripts/select_best_assemblies.py")
    assemblies_tsv = select_best_assemblies(json_summaries_dir_ch, script2)

    // // //4. download proteomes
    script3 = channel.fromPath("scripts/download_genomes.py")
    proteomes_ch = download_genomes(assemblies_tsv, script3)
    
    //5. rCreate gene to species mapping
    script4 = channel.fromPath('scripts/create_mapping_gene_to_specie.py')
    create_mapping(proteomes_ch.sequences, script4)

    //6. run mmseqs2 to make 'BLAST' and create gene families
    script5 = channel.fromPath('scripts/make_BLAST-like_clusters.py')
    clusters_ch = run_mmseqs(proteomes_ch.sequences, script5)

    //7. filter orthologous genes (1-to-1)
    script6 = channel.fromPath('scripts/filter_cluster.py')
    filter_orthologs(clusters_ch.clusters, script6)
}


process select_genomes {
    // Copy output to that path
    publishDir "data/proteomes", mode: 'copy'

    input:
    path taxonomy
    path script

    output:
    path "selection"


    script:
    """
    ls ./
    bash ${script} ${taxonomy} "selection"
    """
}


process select_best_assemblies {
    publishDir "data/proteomes", mode: 'copy'

    input:
    path json_files
    path script

    output:
    file "selected_assemblies.tsv"

    script:
    """
    python3 ${script} \
        --input ${json_files} \
        --output selected_assemblies.tsv
    """
}


process download_genomes {
    publishDir "data/proteomes", mode: 'copy'

    input:
    path assemblies_tsv
    path script

    
    output:
    path 'zipped', emit: zipped
    path 'sequences', emit: sequences

    script:
    """
    python3 ${script} \
        --input ${assemblies_tsv} \
        --output_zipped zipped \
        --output_sequences sequences
    """
}

process create_mapping {
    publishDir "results/mapping", mode: 'copy'
    
    input:
    path proteome_sequences
    path script
    
    output:
    file 'gene_to_species.tsv'

    script:
    """
    python3 ${script} \
        --input ${proteome_sequences} \
        --output gene_to_species.tsv
    """
}


process run_mmseqs {

    publishDir "results/clusters/mmseqs2", mode: 'copy'

    input:
    path proteome_sequences
    path script

    output:
    path "clusters_cluster.tsv", emit: clusters
    path "clusters_rep_seq.fasta", emit: repseq
    path "clusters_all_seqs.fasta", emit: allseqs
    path "all_proteomes.faa", emit: allproteomes

    script:
    """
    mkdir -p tmp

    python3 ${script} \
        --input ${proteome_sequences} \
        --output clusters \
        --tmp tmp
    """
}

process filter_orthologs {
    publishDir "results/clusters"
    
    input:
    path clusters
    path script

    output:
    file "orthologs1to1.tsv"

    script:
    """
    python3 ${script} \
        --input ${clusters} \
        --output orthologs1to1.tsv
    """
}



// // process  {
// //     input:
// //     path
    
// //     output:
// //     path

// //     script:
// //     """
    
// //     """
// // }