#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYRODIGAL } from '../modules/pyrodigal'
include { ORTHOFINDER } from '../modules/orthofinder'
include { CLUSTER_ORTHOLOGS } from '../modules/orthofinder'

workflow stage_2 {
    take:
    rep_genomes

    main:
    rep_genome_ch = rep_genomes
        .flatten()
        .map { genome -> tuple(genome.baseName, genome) }

    PYRODIGAL(rep_genome_ch)

    protein_ch = PYRODIGAL.out.protein_file
        .map { faa -> tuple(faa.name.replaceFirst(/\.faa$/, ''), faa) }
    gff_ch = PYRODIGAL.out.gff_file
        .map { gff -> tuple(gff.name.replaceFirst(/\.gff$/, ''), gff) }
    gene_ch = PYRODIGAL.out.gene_file
        .map { gene -> tuple(gene.name.replaceFirst(/\.genes\.fa$/, ''), gene) }

    assets_ch = rep_genome_ch
        .join(protein_ch)
        .join(gff_ch)
        .join(gene_ch)
    
    ORTHOFINDER(PYRODIGAL.out.protein_file.collect(), file("${projectDir}/bin/OGFilter.py"))
    
    filtered_og_ch = ORTHOFINDER.out.filtered_OG
        .flatten()
        .map { OG -> tuple(OG.baseName, OG) }

    CLUSTER_ORTHOLOGS(ORTHOFINDER.out.total_OG)

    emit:
    assets_ch = assets_ch
    clean_OG = CLUSTER_ORTHOLOGS.out.clean_OG
    filtered_og_ch = filtered_og_ch
}