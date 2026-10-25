#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 [ -u user-msa ] [ -m msa ] [ -r ref-group ] [ -h help ]"
    echo "  -u, --user-msa             GTDB user msa file"
    echo "  -m, --msa                  GTDB msa file"
    echo "  -t, --outgroup-taxonomy    outgroup taxonomy"
    echo "  -r, --ref-group            Reference group"
    echo "  -h, --help                 Show this help"
    exit 1
}

msa="gtdbtk/align/gtdbtk.bac120.msa.fasta.gz"
user_msa="gtdbtk/align/gtdbtk.bac120.user_msa.fasta.gz"
ref_group="gtdb_rep_genomes/representative_genome.tsv"
outgroup_taxon="o__Flavobacteriales_H"

OPTS=$(getopt -o u:m:r:t:h -l user-msa:,msa:,ref-group:,outgroup-taxonomy:,help -n "$0" -- "$@") || usage
eval set -- "$OPTS"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u|--user-msa)    user_msa="$2";   shift 2 ;;
        -m|--msa)         msa="$2";        shift 2 ;;
        -r|--ref-group)   ref_group="$2";  shift 2 ;;
        -t|--outgroup-taxonomy)    outgroup_taxon="$2";  shift 2 ;;
        -h|--help)        usage ;;
        --) shift;        break ;;
        *)  echo "Internal error: $1" >&2; exit 1 ;;
    esac
done

zcat ${msa} | \
    grep -w -A 1 -f <(cut -f1 "$ref_group") | \
    sed "/^--/d" | \
    gzip -c > gtdbtk.temp.msa.fasta.gz

cat gtdbtk.temp.msa.fasta.gz \
    ${user_msa} \
    > gtdbtk.combine_msa.fasta.gz

gtdbtk infer --msa_file gtdbtk.combine_msa.fasta.gz \
    --out_dir gtdbtk.infer

gtdbtk root --input_tree gtdbtk.infer/gtdbtk.unrooted.tree \
    --outgroup_taxon "${outgroup_taxon}" \
    --output_tree gtdbtk.infer/gtdbtk.rooted.tree

sed -i -E 's/(RS|GB)_(GC[AF]_)/\2/g' gtdbtk.infer/gtdbtk.rooted.tree