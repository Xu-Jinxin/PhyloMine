process RGI {
    tag "rgi"
    label "rgi"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path clean_OG

    output:
    path "card", emit: card_dir

    script:
    def container_db_path = "/opt/card_db"
    
    """
    mkdir -p card 
    rgi load --card_json ${container_db_path}/card.json --local

    rgi main --input_type protein  \
        --input_sequence ${clean_OG} \
        --output_file card \
        --alignment_tool DIAMOND \
        --threads ${task.cpus} \
        --clean 
    """

    stub:
    """
    mkdir -p rgi
    """
}

