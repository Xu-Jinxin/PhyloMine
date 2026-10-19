process ANTISMASH {
    tag "${id}"
    label "antismash"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple val(id), path(genome), path(faa), path(gff), path(gene)

    output:
    path "antismash/${id}", emit: antismash_dir

    script:
    def container_db_path = '/opt/antismash_db'
    """
    mkdir -p antismash/${id}
    antismash ${genome} \
        --cpus ${task.cpus} \
        --databases ${container_db_path} \
        --genefinding-tool none \
        --genefinding-gff3 ${gff} \
        --output-dir antismash/${id}
    """
    
    stub:
    """
    mkdir -p antismash/${id}
    touch antismash/${id}/stub.txt
    """

}
