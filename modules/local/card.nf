process RGI {
    tag "rgi"
    label "rgi"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path clean_OG

    output:
    path "rgi", emit: rgi_dir

    script:
    def container_db_path = '/opt/card_db'
    
    """
    rgi main --input_type protein  \
        --input_sequence ${clean_OG} \
        --output_file rgi/card \
        --alignment_tool DIAMOND \
        --threads ${task.cpus} \
        --local ${container_db_path} \
        --clean 
    """

    stub:
    """
    mkdir -p rgi
    """
}

