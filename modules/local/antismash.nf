process ANTISMASH {
    tag "${id}"
    label "antismash"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple val(id), path(rep_genome)

    output:
    path "antismash/${id}", emit: antismash_dir

    script:
    """
    mkdir -p antismash/${id}
    antismash gencomp/skder/Dereplicated_Representative_Genomes/${id}.fa \
        --cpus 16 \
        --databases /share/data01/project/xujinxin/database/antismash \
        --genefinding-tool none \
        --genefinding-gff3 gencomp/pyrodigal/gff/${id}.gff \
        --output-dir gencomp/antismash/${id}
    """
    
    stub:
    """
    mkdir -p pyrodigal/{gff,genes,proteins}
    touch pyrodigal/gff/${id}.gff
    touch pyrodigal/genes/${id}.fa
    touch pyrodigal/proteins/${id}.faa
    """

}
