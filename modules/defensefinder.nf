process DEFENSEFINDER {
    tag "${id}"
    label "defensefinder"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple val(id), path(genome), path(faa), path(gff), path(gene)

    output:
    path "defensefinder/${id}", emit: defensefinder_dir

    script:
    def macsyfinder_models_path = '/opt/macsyfinder'

    """
    defense-finder run ${faa} \
        --out-dir defensefinder/${id} \
        --workers ${task.cpus} \
        --skip-model-version-check \
        --models-dir ${macsyfinder_models_path}

    """

    stub:
    """
    mkdir -p defensefinder/${id}
    """
}