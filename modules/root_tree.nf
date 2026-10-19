process ROOT_TREE {
    tag "root_tree"
    label "root_tree"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path gtdbtk_dir
    path ref_group

    output:
    path "gtdbtk.infer/gtdbtk.rooted.tree", emit: species_root_tree

    script:
    """
    export GTDBTK_DATA_PATH='/opt/gtdbtk_db'
    outgroup_taxon=\$(tail -n1 ${ref_group} | cut -f5)

    bash root_tree.sh -u ${gtdbtk_dir}/align/gtdbtk.*.user_msa.fasta.gz \
        -m ${gtdbtk_dir}/align/gtdbtk.*.msa.fasta.gz \
        -r ${ref_group} \
        -t \${outgroup_taxon}

    """

    stub:
    """
    mkdir -p gtdbtk.infer
    touch gtdbtk.infer/gtdbtk.rooted.tree
    """

}