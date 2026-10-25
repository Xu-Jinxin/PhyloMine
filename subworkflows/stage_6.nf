#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IQTREE_OG } from '../modules/iqtree_og'
include { ALERAX } from '../modules/alerax'

workflow stage_6 {
    take:
    filtered_og_ch
    species_root_tree
    gtdb_rep_genomes_list
    checkm2_report

    main:
    IQTREE_OG(filtered_og_ch)

    families_ch = filtered_og_ch
        .join(IQTREE_OG.out.iqtree_og_tree
            .map { tree -> tuple(tree.baseName, tree) })
        .toList()
        .map { items -> tuple(items.collect { item -> item[1] }, items.collect { item -> item[2] }) }

    ALERAX(families_ch,
           species_root_tree,
           gtdb_rep_genomes_list,
           checkm2_report)

}