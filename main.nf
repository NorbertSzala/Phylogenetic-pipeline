#!/usr/bin/env nextflow

// nextflow run main.nf -resume

nextflow.enable.dsl=2

params.taxonomy = "data/short_taxonomy.csv"
params.scripts = "scripts"
params.results = 'results   '

workflow {
    // 1. input taxonomy file
    taxonomy_ch = Channel.fromPath(params.taxonomy)

    // //2. select genome assemblies (JNOS summaries)
    json_summaries_ch = select_genomes(taxonomy_ch)
    
    // //3. select best assemblies
    select_best_assemblies(json_summaries_ch)

    // //4. download proteomes
    // proteomes_ch = download_genomes(assemblies_ch)
    
    // //5. rCreate gene to species mapping
    // mapping_ch = create_mapping(proteomes_ch)

    // //6. run mmseqs2 to make 'BLAST' andcreate gene families
    // clusters_ch = run_mmseqs(proteomes_ch)

    // //7. filter orthologous genes (1-to-1)
    // filter_orthologs(clusters_ch)
}


process select_genomes {

    publishDir "data/proteomes/", mode: 'copy'

    input:
    path taxonomy

    output:
    path selection

    script:
    """
    ls ./
    bash ${projectDir}/scripts/select_genomes.sh \
        ${taxonomy} \
        selection
    """
}


process select_best_assemblies {
    publishDir "data/proteomes", mode: 'copy'

    input:
    path selection_dir

    output:
    path "nextflow_assemblies.tsv"

    script:
    """
    python3 ${projectDir}/scripts/select_best_assemblies.py \
        --input ${selection_dir} \
        --output nextflow_assemblies.tsv
    """
}


// process download_genomes {
//     input:
//     path  data/proteomes/selected_assemblies.tsv

    
//     output:
//     path "data/proteomes/sequences/*.faa"

//     script:
//     """
//     python3 ${projectDir}/scripts/download_genomes.py
//     """
// }

// process create_mapping {
//     input:
//     path "*.faa"
    
//     output:
//     path "results/clusters/gene_to_species.tsv"

//     script:
//     """
//     python3 ${projectDir}/scripts/create_mapping_gene_to_species.py
//     """
// }


// process run_mmseqs {

//     input:
//     path "*.faa"

//     output:
//     path "results/clusters/mmseqs2/*"

//     script:
//     """
//     python ${projectDir}/scripts/make_BLAST-like_clusters.py
//     """
// }

// process filter_orthologs {

//     input:
//     path "clusters_cluster.tsv"

//     output:
//     path "results/clusters/orthologs1to1.tsv"

//     script:
//     """
//     python ${projectDir}/scripts/filter_cluster.py
//     """
// }



// // process  {
// //     input:
// //     path
    
// //     output:
// //     path

// //     script:
// //     """
    
// //     """
// // }


// # TODO: przetestowac nextflow na wszystkich etapach do tej pory
// # TODO: ustawic logowanie we wszystkich skryptach do katalogu outputu bo nf tego nie lubi