process ALERAX {
    tag "alerax"
    label "alerax"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    tuple path(og_fastas, stageAs: "filtered_OG/*"), path(gene_trees, stageAs: "gene_trees/*")
    path species_root_tree
    path representative_genome_tsv
    path checkm2_report

    output:
    path "alerax", emit: alerax_out
    path "family.txt", emit: family_file
    path "missing.txt", emit: missing_file

    script:
    """
    mkdir -p mapping

    {
        awk -F '\t' -v OFS='\t' '{name = \$1; sub(/^(RS|GB)_/, "", name); print name, (100 - \$4) / 100}' ${representative_genome_tsv}
        awk -F '\t' -v OFS='\t' 'NR > 1 && \$2 >= ${params.min_completeness} && \$3 <= ${params.max_contamination} {print \$1, (100 - \$2) / 100}' ${checkm2_report}
    } > missing.txt

    echo "[FAMILIES]" > family.txt

    for fa in filtered_OG/*.fa; do
        OG=\$(basename "\${fa}" .fa)

        # gene headers look like >GCA_xxx__seq1_1; the species id is the part before "__"
        grep "^>" "\${fa}" | sed "s/^>//" | \
            awk -F "__" -v OFS="\t" '{print \$0, \$1}' > mapping/\${OG}.mapping

        echo "- \${OG}" >> family.txt
        echo "gene_tree = gene_trees/\${OG}.ufboot" >> family.txt
        echo "mapping = mapping/\${OG}.mapping" >> family.txt
    done

    mpiexec --use-hwthread-cpus --oversubscribe -np ${task.cpus} \
        alerax -f family.txt \
        -s ${species_root_tree} \
        -p alerax \
        --fraction-missing-file missing.txt \
        --gene-tree-samples 100 \
        --gene-tree-rooting MAD \
        --species-tree-search SKIP \
        --no-tl \
        --transfer-constraint PARENTS \
        --model-parametrization GLOBAL \
        --origination LCA
    """

    stub:
    """
    mkdir -p alerax
    touch family.txt missing.txt
    """
}