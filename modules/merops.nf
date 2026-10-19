process MEROPS {
    tag "merops"
    label "merops"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path clean_OG

    output:
    path "merops", emit: merops_dir

    script:
    def container_db_path = '/opt/merops_db'
    
    """
    mkdir -p merops

    diamond blastp --db ${container_db_path} \
        --query ${clean_OG} \
        --out merops/MEROPS.tsv \
        --outfmt 6 \
        --max-target-seqs 1 \
        --threads ${task.cpus} \
    """

    stub:
    """
    mkdir -p merops
    """
}