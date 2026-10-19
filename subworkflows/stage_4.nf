#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CGCFINDER } from '../modules/cgcfinder'
include { ANTISMASH } from '../modules/antismash'
include { BIGSCAPE } from '../modules/bigscape'

workflow stage_4 {
    take:
    genome_assets_ch

    main:
    CGCFINDER(genome_assets_ch)
    ANTISMASH(genome_assets_ch)
    BIGSCAPE(ANTISMASH.out.antismash_dir.collect())
}