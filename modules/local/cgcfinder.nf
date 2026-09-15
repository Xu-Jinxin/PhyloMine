process CGCFINDER {
    tag "${id}"
    label "cgc"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple val(id), path(faa), path(gff)

    output:
    path "cgc/${id}", emit: cgc_dir

    script:
    """
    mkdir -p cgc/${id}

    run_dbcan ${faa} protein \
        --cluster ${gff} \
        --out_dir cgc/${id}
    """

    stub:
    """
    mkdir -p cgc/${id}
    touch cgc/${id}/stub.txt
    """
}