process GENOMESPOT {
    tag "genomespot"
    label "genomespot"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple val(id), path(genome), path(faa), path(gff), path(gene)

    output:
    path "genomespot/${id}.predictions.tsv", emit: genomespot_file

    script:
    def container_db_path = '/opt/genomespot_db'

    """
    mkdir -p genomespot
    python -m genome_spot.genome_spot \
        --models ${container_db_path} \
        --contigs  ${genome} \
        --proteins  ${faa} \
        --output genomespot/${id}

    """

    stub:
    """

    """

}