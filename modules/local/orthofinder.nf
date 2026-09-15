process ORTHOFINDER {
    tag "orthofinder"
    label "orthofinder"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path "input_proteins/*"

    output:
    path "orthofinder", emit: orthofinder_dir

    script:
    """
    orthofinder -f input_proteins -t ${task.cpus} \
        -M msa -S diamond -A mafft -T fasttree \
        -o orthofinder
    """

    stub:
    """
    mkdir -p orthofinder
    """

}

process CLUSTER_ORTHOLOGS{
    tag "cluster_OG"
    label "cluster_OG"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path orthofinder_dir

    output:
    path "total_OG.clean.70.fa", emit: clean_OG

    script:
    """
    results_dirs=( ${orthofinder_dir}/Results_* )
    if [ \${#results_dirs[@]} -ne 1 ]; then
        echo "ERROR: Expected exactly one Results_* directory, found \${#results_dirs[@]}" >&2
        exit 1
    fi

    cat "\${results_dirs}/Orthogroup_Sequences/"*.fa > total_OG.fa
    sed '/^>/! s/[UuXx*]//g' total_OG.fa > total_OG.clean.fa
    cd-hit -i total_OG.clean.fa -o total_OG.clean.70.fa -d 0 -c 0.7 -T 4 -n 5

    """
    stub:
    """
    touch total_OG.clean.70.fa
    """
}