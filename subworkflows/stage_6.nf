#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IQTREE_OG } from '../modules/iqtree_og'

workflow stage_6 {
    take:
    filtered_og_ch

    main:
    IQTREE_OG(filtered_og_ch)

}