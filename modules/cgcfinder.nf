process CGCFINDER {
    tag "${id}"
    label "cgc"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple val(id), path(genome), path(faa), path(gff), path(gene)

    output:
    path "cgc/${id}", emit: cgc_dir

    script:
    def container_db_path = '/opt/dbcan_db'
    """
    mkdir -p cgc/${id}

    run_dbcan CAZyme_annotation \
        --input_raw_data ${faa} \
        --mode protein \
        --output_dir cgc/${id}  \
        --db_dir ${container_db_path}

    run_dbcan gff_process \
    --input_gff ${gff} \
    --gff_type prodigal \
    --output_dir cgc/${id} \
    --db_dir ${container_db_path} 

    run_dbcan cgc_finder --output_dir cgc/${id}

    run_dbcan substrate_prediction \
        --input_raw_data ${faa} \
        --mode protein \
        --output_dir cgc/${id} \
        --db_dir ${container_db_path}

    run_dbcan cgc_circle_plot \
        --output_dir cgc/${id}
    """

    stub:
    """
    mkdir -p cgc/${id}
    touch cgc/${id}/stub.txt
    """
}