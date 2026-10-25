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
    if (params.run_stage == null) {
        genomes_ch = Channel.fromPath(params.input)
        stage_1(genomes_ch)
        stage_2(stage_1.out.rep_genomes)
        stage_3(stage_2.out.clean_OG)
        stage_4(stage_2.out.assets_ch)
        stage_5(stage_2.out.assets_ch)
        stage_6(stage_2.out.filtered_og_ch, 
                stage_1.out.species_root_tree, 
                stage_1.out.gtdb_rep_genomes_list, 
                stage_1.out.checkm2_report)

    } else if (params.run_stage == 'stage2'){
        stage_2_only()
    } else if (params.run_stage == 'stage3') {
        stage_3_only()
    } else if (params.run_stage == 'stage4') {
        stage_4_only
    }
}

workflow stage_2_only {
    rep_genomes = Channel.fromPath(params.genome_dir)
    stage_2(rep_genomes)
}

workflow stage_3_only {
    clean_OG = Channel.fromPath(params.clean_og)
    stage_3(clean_OG)
}

workflow stage_4_only {
    genome_ch = Channel.fromPath("${params.genome_dir}/*.${params.genome_extension}")
        .map { genome -> tuple(genome.baseName, genome) }
    protein_ch = Channel.fromPath("${params.proteins_dir}/*.faa")
        .map { faa -> tuple(faa.name.replaceFirst(/\.faa$/, ''), faa) }
    gff_ch = Channel.fromPath("${params.gff_dir}/*.gff")
        .map { gff -> tuple(gff.name.replaceFirst(/\.gff$/, ''), gff) }
    gene_ch = Channel.fromPath("${params.genes_dir}/*.genes.fa")
        .map { gene -> tuple(gene.name.replaceFirst(/\.genes\.fa$/, ''), gene) }

    assets_ch = genome_ch
        .join(protein_ch)
        .join(gff_ch)
        .join(gene_ch)

    stage_4(assets_ch)
}

workflow stage_5_only {
    genome_ch = Channel.fromPath("${params.genome_dir}/*.${params.genome_extension}")
        .map { genome -> tuple(genome.baseName, genome) }
    protein_ch = Channel.fromPath("${params.proteins_dir}/*.faa")
        .map { faa -> tuple(faa.name.replaceFirst(/\.faa$/, ''), faa) }
    gff_ch = Channel.fromPath("${params.gff_dir}/*.gff")
        .map { gff -> tuple(gff.name.replaceFirst(/\.gff$/, ''), gff) }
    gene_ch = Channel.fromPath("${params.genes_dir}/*.genes.fa")
        .map { gene -> tuple(gene.name.replaceFirst(/\.genes\.fa$/, ''), gene) }

    assets_ch = genome_ch
        .join(protein_ch)
        .join(gff_ch)
        .join(gene_ch)

    stage_5(assets_ch)
}
