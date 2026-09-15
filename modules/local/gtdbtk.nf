process GTDBTK {
    tag "gtdbtk"
    label "gtdbtk"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path rep_genomes_dir

    output:
    path "gtdbtk", emit: gtdbtk_dir

    script:
    """
    export GTDBTK_DATA_PATH='/opt/gtdbtk_db'
    gtdbtk classify_wf \
        --genome_dir ${rep_genomes_dir} \
        --out_dir  gtdbtk \
        --place_species \
        --extension "${params.genome_extension}" \
        --cpus "${task.cpus}" --force
    """

    stub:
    """
    mkdir -p gtdbtk
    """

}