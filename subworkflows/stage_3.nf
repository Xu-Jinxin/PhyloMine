#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAKTA } from '../modules/bakta'
include { EGGNOG } from '../modules/eggNOG'
include { CAZYME } from '../modules/cazyme'
include { RGI } from '../modules/card'
include { VF } from '../modules/vf'
include { MEROPS } from '../modules/merops'

workflow stage_3 {
    take:
    clean_OG

    main:
    BAKTA(clean_OG)
    EGGNOG(clean_OG)
    CAZYME(clean_OG)
    RGI(clean_OG)
    VF(clean_OG)
    MEROPS(clean_OG)
}