process PYRODIGAL {
    tag "${id}"
    label "pyrodigal"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple val(id), path(rep_genome)

    output:
    path "pyrodigal/gff/${id}.gff", emit: gff_file
    path "pyrodigal/genes/${id}.fa", emit: gene_file
    path "pyrodigal/proteins/${id}.faa", emit: protein_file

    script:
    """
    mkdir -p pyrodigal/{gff,genes,proteins}
    pyrodigal -p meta -f gff \
        -i ${rep_genome} \
        -o pyrodigal/gff/${id}.gff \
        -d pyrodigal/genes/${id}.fa \
        -a pyrodigal/proteins/${id}.faa \
        -j ${task.cpus}
    """
    
    stub:
    """
    mkdir -p pyrodigal/{gff,genes,proteins}
    touch pyrodigal/gff/${id}.gff
    touch pyrodigal/genes/${id}.fa
    touch pyrodigal/proteins/${id}.faa
    """

}