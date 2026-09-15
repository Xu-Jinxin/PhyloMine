 process CAZYME {
    tag "cazyme"
    label "cazyme"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path clean_OG

    output:
    path "CAZyme", emit: cazy_dir

    script:
    def container_db_path = '/opt/dbcan_db'
    """
    run_dbcan ${clean_OG} protein \
        --hmm_cpu ${task.cpus} \
        --dia_cpu ${task.cpus} \
        --tools hmmer diamond \
        --out_dir CAZyme  \
        --db_dir ${container_db_path}
    """

    stub:
    """
    mkdir -p CAZyme
    """
}