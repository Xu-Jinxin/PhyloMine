process HGTECTOR2 {
    tag "${id}"
    label "hgtector2"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple val(id), path(genome), path(faa), path(gff), path(gene)

    output:
    path "hgtector2/${id}", emit: hgtector2_dir

    script:
    def container_db_path = '/opt/hgtector2_db'

    """
    mkdir -p hgtector2/{search,analyze}

    hgtector search \
        -i ${faa} \
        -o hgtector2/${id}/search \
        -m diamond \
        -p ${task.cpus} \
        -d ${container_db_path}/db \
        -t ${container_db_path}/taxdump

    hgtector analyze \
        -i hgtector2/${id}/search \
        -o hgtector2/${id}/analyze \
        -t ${container_db_path}/taxdump

    """

    stub:
    """
    mkdir -p hgtector2/${id}
    touch hgtector2/${id}/stub.txt
    """
}