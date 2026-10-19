process ALERAX {
    tag "alerax"
    label "alerax"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path "iqtree_og/*"
    path "species_root_tree"

    output:
  
    script:
    
    """
    cat <(cut -f1,4 gencomp/gtdb_rep_genomes/representative_genome.tsv) \
        <(awk -F '\t' '($2 >= 70 && $3 <= 5)' gencomp/checkm2/quality_report.tsv | cut -f1,2) | \
        awk 'BEGIN{FS=OFS="\t"}{print $1, 100-$2}' > missing.txt
    
    ls gencomp/orthofinder/OGfilter/*.fa | while read file; do 
        OG=$(basename ${file} | sed "s/.fa//")
        grep "^>" gencomp/orthofinder/OGfilter/${OG}.fa | sed "s/>//" | \
            awk -F "__" '{print $0 "\t" $1 }' \
            > test/${OG}.mapping
    done

    echo "[FAMILIES]" > family.txt
    ls gencomp/orthofinder/OGfilter/*.fa | while read file; do 
        OG=$(basename ${file} | sed "s/.fa//")
        echo "- ${OG}" >> family.txt
        echo "gene_tree = gencomp/orthofinder/OGfilter/${OG}.ufboot" >> family.txt
        echo "mapping = test/${OG}.mapping" >> family.txt
    done
    
    mpiexec -np ${task.cpus} \
        alerax -f family.txt \
        -s ${species_root_tree} \
        -p output_3_per_family \
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

    """
}