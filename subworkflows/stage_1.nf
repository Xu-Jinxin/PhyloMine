#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKM2 } from '../modules/checkm2'
include { FILTER_GENOMES } from '../modules/checkm2'
include { SKDER } from '../modules/skder'
include { GTDBTK } from '../modules/gtdbtk'
include { GTDB_SELECT_REP } from '../modules/gtdb_select_representative'
include { ROOT_TREE } from '../modules/root_tree'

workflow stage_1 {
    take:
    genomes_ch

    main:
    CHECKM2(genomes_ch)
    FILTER_GENOMES(CHECKM2.out.checkm2_dir, genomes_ch)
    SKDER(FILTER_GENOMES.out.high_quality_genomes_dir)
    GTDBTK(SKDER.out.drep_genomes_dir)
    GTDB_SELECT_REP(SKDER.out.drep_genomes_dir, GTDBTK.out.gtdbtk_dir)
    ROOT_TREE(GTDBTK.out.gtdbtk_dir, GTDB_SELECT_REP.out.gtdb_rep_genomes_list)
    
    emit:
    rep_genomes = GTDB_SELECT_REP.out.rep_genomes
    gtdbtk_out  = GTDBTK.out
    species_root_tree = ROOT_TREE.out.species_root_tree
    gtdb_rep_genomes_list = GTDB_SELECT_REP.out.gtdb_rep_genomes_list
    checkm2_report = CHECKM2.out.quality_report
}
