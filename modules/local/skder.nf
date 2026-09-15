process SKDER {
    tag "skder"
    label "skder"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path high_quality_genomes_dir

    output:
    path "skder", emit: skder_dir
    path "skder/Dereplicated_Representative_Genomes", emit: rep_genomes_dir
    path "skder/Dereplicated_Representative_Genomes/*", emit: rep_genomes

    script:
    """
    skder -g ${high_quality_genomes_dir} \
        -o skder \
        --threads ${task.cpus} \
        --max-memory ${task.memory.toGiga()} \
        -d dynamic
    """
    
    stub:
    """
    mkdir -p high_quality_genomes
    mkdir -p skder/Dereplicated_Representative_Genomes
    touch skder/Dereplicated_Representative_Genomes/genome{1..9}.fasta
    """
}