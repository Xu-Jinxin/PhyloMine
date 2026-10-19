process EGGNOG {
    tag "eggNOG"
    label "eggNOG"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path clean_OG

    output:
    path "eggNOG", emit: eggNOG_dir

    script:
    def container_db_path = '/opt/eggNOG_db'
    """
    mkdir -p eggNOG
    emapper.py -i ${clean_OG} \
        -o eggNOG/OG \
        --data_dir ${container_db_path} \
        --cpu ${task.cpus}
    """

    stub:
    """
    mkdir -p eggNOG
    """
}