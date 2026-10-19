#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEFENSEFINDER } from '../modules/defensefinder'
include { CCTYPER } from '../modules/crisprcastyper'
include { HGTECTOR2 } from '../modules/hgtector2'

workflow stage_5 {
    take:
    genome_assets_ch

    main:
    DEFENSEFINDER(genome_assets_ch)
    CCTYPER(genome_assets_ch)
    HGTECTOR2(genome_assets_ch)

}