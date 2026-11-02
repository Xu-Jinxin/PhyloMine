process BAKTFOLD {
    tag "baktfold"
    label "baktfold_gpu"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path bakta_json

    output:
    path "baktfold", emit: baktfold_dir

    script:
    def container_db_path = '/opt/baktfold_db'
    """
    baktfold proteins     \
        --input ${bakta_json} \
        --output baktfold \
        --threads 8 \
        --database ${container_db_path} \
        --foldseek-gpu

    """

    stub:
    """
    mkdir -p baktfold
    """

}