#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// if (!params.input || !params.outdir) {
//     exit 1, "Error ！\nUsage : nextflow run main.nf --input /path/to/genomes_dir --outdir /path/to/output\n"
// }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// import functions 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

include { CHECKM2 } from './modules/local/checkm2'
include { FILTER_GENOMES } from './modules/local/checkm2'
include { SKDER } from './modules/local/skder'
include { GTDBTK } from './modules/local/gtdbtk'
include { PYRODIGAL } from './modules/local/pyrodigal'
include { ORTHOFINDER } from './modules/local/orthofinder'
include { CLUSTER_ORTHOLOGS } from './modules/local/orthofinder'
include { BAKTA } from './modules/local/bakta'
include { EGGNOG } from './modules/local/eggNOG'
include { CAZYME } from './modules/local/cazyme'
include { RGI } from './modules/local/card'
include { VF } from './modules/local/vf'
include { CGCFINDER } from './modules/local/cgcfinder'


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// run main workflow
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

workflow {
    genomes_dir = Channel.fromPath("${params.input}")
    CHECKM2(genomes_dir)
    FILTER_GENOMES(CHECKM2.out.checkm2_dir, genomes_dir)
    SKDER(FILTER_GENOMES.out.high_quality_genomes_dir)
    GTDBTK(SKDER.out.rep_genomes_dir)
    
    rep_genome = SKDER.out.rep_genomes
        .flatten()
        .map { genome -> tuple(genome.baseName, genome) }
    
    PYRODIGAL(rep_genome)
    ORTHOFINDER(PYRODIGAL.out.protein_file.collect())
    CLUSTER_ORTHOLOGS(ORTHOFINDER.out.orthofinder_dir)
    BAKTA(CLUSTER_ORTHOLOGS.out.clean_OG)
    EGGNOG(CLUSTER_ORTHOLOGS.out.clean_OG)
    CAZYME(CLUSTER_ORTHOLOGS.out.clean_OG)
    RGI(CLUSTER_ORTHOLOGS.out.clean_OG)
    VF(CLUSTER_ORTHOLOGS.out.clean_OG)

    protein_ch = PYRODIGAL.out.protein_file
        .map { faa -> tuple(faa.baseName, faa) }
    gff_ch = PYRODIGAL.out.gff_file
        .map { gff -> tuple(gff.baseName, gff) }

    cgc_input = protein_ch.join(gff_ch)
        .map { id, faa, gff -> tuple(id, faa, gff) }

    CGCFINDER(cgc_input)

}