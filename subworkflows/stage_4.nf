#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CGCFINDER } from '../modules/cgcfinder'
include { ANTISMASH } from '../modules/antismash'
include { BIGSCAPE } from '../modules/bigscape'
include { DEFENSEFINDER } from '../modules/defensefinder'
include { CCTYPER } from '../modules/crisprcastyper'
include { HGTECTOR2 } from '../modules/hgtector2'
include { GENOMESPOT } from '../modules/genomespot'

workflow stage_4 {
    take:
    genome_assets_ch

    main:
    CGCFINDER(genome_assets_ch)
    ANTISMASH(genome_assets_ch)
    BIGSCAPE(ANTISMASH.out.antismash_dir.collect())
    DEFENSEFINDER(genome_assets_ch)
    CCTYPER(genome_assets_ch)
    HGTECTOR2(genome_assets_ch)
    GENOMESPOT(genome_assets_ch)

}