process BAKTA {
    tag "bakta"
    label "bakta"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path clean_OG

    output:
    path "bakta", emit: bakta_dir

    script:
    def container_db_path = '/opt/bakta_db'
    """
    bakta_proteins \
        --db ${container_db_path} \
        --output bakta \
        "${clean_OG}"
    """

    stub:
    """
    mkdir -p bakta
    """

}