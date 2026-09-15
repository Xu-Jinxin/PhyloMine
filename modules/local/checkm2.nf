
process CHECKM2 {
    tag "checkm2"
    label "checkm2"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path genomes_dir

    output:
    path "checkm2", emit: checkm2_dir

    script:
    def container_db_path = '/opt/checkm2_db/uniref100.KO.1.dmnd'
    """
    checkm2 predict \
        --input "${genomes_dir}" \
        --output-directory checkm2 \
        --database_path "${container_db_path}" \
        --threads ${task.cpus} \
        -x ${params.genome_extension} \
        --force
    """

    stub:
    """
    mkdir -p checkm2
    touch checkm2/quality_report.tsv
    """
}

process FILTER_GENOMES {
    tag "filter_genomes"
    label "filter_genomes"
    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path checkm2_dir
    path genomes_dir

    output:
    path "high_quality_genomes", emit: high_quality_genomes_dir

    script:
    """
    mkdir -p high_quality_genomes
    awk -F '\t' '(\$2 >= ${params.min_completeness} && \$3 <= ${params.max_contamination})' \
        ${checkm2_dir}/quality_report.tsv | cut -f1 | while read genome; do
        cp ${genomes_dir}/\${genome}.${params.genome_extension} high_quality_genomes/
    done
    """

    stub:
    """
    mkdir -p high_quality_genomes
    """
}
