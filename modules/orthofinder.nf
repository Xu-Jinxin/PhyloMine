process ORTHOFINDER {
    tag "orthofinder"
    label "orthofinder"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path "input_proteins/*"
    path "OGFilter.py", stageAs: "OGFilter.py"

    output:
    path "orthofinder", emit: orthofinder_dir
    path "orthofinder/OGfilter/*.fa", emit: filtered_OG
    path "total_OG.clean.fa", emit: total_OG 

    script:
    """
    orthofinder -f input_proteins -t ${task.cpus} \
        -M msa -S diamond -A mafft -T fasttree -os \
        -o orthofinder

    python OGFilter.py \
        -g orthofinder/Results_*/Orthogroups/Orthogroups.GeneCount.tsv \
        -s orthofinder/Results_*/Orthogroup_Sequences \
        -o orthofinder/OGfilter \
        -max_copies 1 \
        -min_species 8
    
    cat orthofinder/OGfilter/OG*.fa > total_OG.fa
    sed '/^>/! s/[UuXx*]//g' total_OG.fa > total_OG.clean.fa

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
    path total_OG

    output:
    path "total_OG_clean_70_rep_seq.fasta", emit: clean_OG

    script:
    """
    mmseqs easy-cluster ${total_OG} total_OG_clean_70 tmp_dir \
        --min-seq-id 0.7 \
        --threads ${task.cpus}

    """
    stub:
    """
    touch total_OG_clean_70_rep_seq.fasta
    """
}