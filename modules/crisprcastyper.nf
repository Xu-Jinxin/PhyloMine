process CCTYPER {
    tag "${id}"
    label "crisprcastyper"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple val(id), path(genome), path(faa), path(gff), path(gene)

    output:
    path "crisprcastyper/${id}", emit: crisprcastyper_dir

    script:
    def container_db_path = '/opt/cctyper_db'
    """
    mkdir -p crisprcastyper
    
    cctyper ${genome} \
        crisprcastyper/${id} \
        --gff ${gff} \
        --prot ${faa} \
        --db ${container_db_path}
    """

    stub:
    """
    mkdir -p crisprcastyper/${id}
    touch crisprcastyper/${id}/stub.txt
    """
}