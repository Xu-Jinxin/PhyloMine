process BIGSCAPE {
    tag "bigscape"
    label "bigscape"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path antismash_dir, stageAs: 'antismash/*'

    output:
    path "bigscape", emit: bigscape_dir

    script:
    def container_db_path = '/opt/pfam_db'
    
    """
    bigscape cluster \
        -i antismash \
        -o bigscape \
        -p ${container_db_path}//Pfam-A.hmm \
        --mix
    """

    stub:
    """
    mkdir -p bigscape
    """

}