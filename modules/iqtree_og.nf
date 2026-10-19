process IQTREE_OG {
    tag "iqtree_og"
    label "iqtree_og"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple val(id), path(og)

    output:
    path "iqtree_og/${id}", emit: iqtree_og_dir
    path "iqtree_og/${id}/${id}.ufboot", emit: iqtree_og_tree

    script:
    """
    mkdir -p iqtree_og/${id}

    mafft --auto ${og} > ${id}.aln.fasta
    iqtree -m LG+I+G -s ${id}.aln.fasta --prefix iqtree_og/${id}/${id} -T ${task.cpus} --wbtl --ufboot 1000
    """

    stub:
    """
    mkdir -p iqtree_og
    """

}