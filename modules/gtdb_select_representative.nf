process GTDB_SELECT_REP {
    tag "gtdb_select_rep"
    label "gtdb_select_rep"

    publishDir "${params.outdir}/gencomp/", mode: 'symlink'

    input:
    path drep_genomes_dir
    path gtdbtk_dir

    output:
    path "rep_genomes/*", emit: rep_genomes
    path "gtdb_rep_genomes/representative_genome.tsv", emit: gtdb_rep_genomes_list

    script:
    """
    taxon=\$(find_common_taxon.sh -p ${gtdbtk_dir})
    
    fetch_rep_genome.sh -m ${params.gtdb_auxillary_metadata} \
        -g ${params.gtdbtk_genome} \
        -t \${taxon} \
        -r ${params.level} \
        -o gtdb_rep_genomes
    
    mkdir -p rep_genomes

    # cp ${drep_genomes_dir}/* rep_genomes
    # ls gtdb_rep_genomes/*_genomic.fna.gz | while read file; do
    #     accession=\$(basename \${file} | sed "s/_genomic.fna.gz//")
    #     zcat \${file} > rep_genomes/\${accession}.${params.genome_extension}
    # done
    # seqkit replace -p '.+' -r "GCA_964619305.1__seq{nr}"  gencomp/rep_genomes/GCA_964619305.1.fa

    for file in ${drep_genomes_dir}/*; do
        base=\$(basename "\${file}")
        name=\${base%.fa}
        seqkit replace -p '.+' -r "\${name}__seq{nr}" "\${file}" \
            > rep_genomes/\${name}.${params.genome_extension}
    done

    for file in gtdb_rep_genomes/*_genomic.fna.gz; do
        accession=\$(basename "\${file}" | sed 's/_genomic.fna.gz//')
        zcat "\${file}" | \
            seqkit replace -p '.+' -r "\${accession}__seq{nr}" \
            > rep_genomes/\${accession}.${params.genome_extension}
    done

    """

    stub:
    """
    mkdir -p gtdbtk
    """

}