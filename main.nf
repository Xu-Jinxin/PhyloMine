#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// if (!params.input || !params.outdir) {
//     exit 1, "Error ！\nUsage : nextflow run main.nf --input /path/to/genomes_dir --outdir /path/to/output\n"
// }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// import subworkflow
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

include { stage_1 } from './subworkflows/stage_1'
include { stage_2 } from './subworkflows/stage_2'
include { stage_3 } from './subworkflows/stage_3'
include { stage_4 } from './subworkflows/stage_4'
include { stage_5 } from './subworkflows/stage_5'
include { stage_6 } from './subworkflows/stage_6'

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// run main workflow
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

workflow {
    genomes_ch = Channel.fromPath(params.input)

    stage_1(genomes_ch)
    stage_2(stage_1.out.rep_genomes)
    stage_3(stage_2.out.clean_OG)
    stage_4(stage_2.out.assets_ch)
    stage_5(stage_2.out.assets_ch)
    stage_6(stage_2.out.filtered_og_ch)
}