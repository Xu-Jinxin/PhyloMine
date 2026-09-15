process VF {
    tag "vf"
    label "vf"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path clean_OG

    output:
    path "vf", emit: vf_dir

    script:
    def container_db_path = '/opt/vf_db'
    
    """
    diamond blastp --db ${container_db_path} \
        --query ${clean_OG} \
        --out vf/VFDB.tsv \
        --outfmt 6 \
        --max-target-seqs 1 \
        --threads ${task.cpus} \
    """

    stub:
    """
    mkdir -p vf
    """
}