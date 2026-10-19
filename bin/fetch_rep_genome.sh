#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 [-m metadata] [-g genomic_reps] [-t taxon] [-r rank] [-o outpath]"
    echo "  -m, --metadata       GTDB metadata TSV"
    echo "  -g, --genomic-reps   GTDB genomic_fasta_reps directory"
    echo "  -t, --taxon          target taxon, e.g. o__Mycobacteriales"
    echo "  -r, --rank           rank, one of: phylum class order family genus species"
    echo "  -o, --outpath        output directory"
    exit 1
}

metadata="/share/data01/project/xujinxin/database/gtdb_r232/auxillary_files/metadata_r232.tsv"
genomic_reps="/home/database/GTDB_R232/genomic_fasta_reps"
taxon="o__Flavobacteriales"
rank="order"
outpath="gtdb_rep_genomes"

# 输出
rep_tsv="${outpath}/representative_genome.tsv"
acc_list="${outpath}/accessions.txt"

OPTS=$(getopt -o m:g:t:r:o:h -l metadata:,genomic-reps:,taxon:,rank:,outpath:,help -n "$0" -- "$@") || usage
eval set -- "$OPTS"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--metadata)       metadata="$2";         shift 2 ;;
        -g|--genomic-reps)   genomic_reps="$2";     shift 2 ;;
        -t|--taxon)          taxon="$2";            shift 2 ;;
        -r|--rank)           rank="$2";             shift 2 ;;
        -o|--outpath)        outpath="$2";          shift 2 ;;
        -h|--help)           usage ;;
        --) shift; break ;;
         *)  echo "Internal error: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$metadata" &&  -n "$genomic_reps" && -n "$taxon" && -n "$rank" && -n "$outpath" ]] || usage

[[ -f "$metadata" ]] || { echo "Error: metadata not found: $metadata" >&2; exit 1; }
[[ -d "$genomic_reps" ]] || { echo "Error: genomic_reps not found: $genomic_reps" >&2; exit 1; }

if [ ! -d "${outpath}" ]; then
    mkdir -p "${outpath}"
fi


declare -A levels=([phylum]="p__" [class]="c__" [order]="o__"
                   [family]="f__" [genus]="g__" [species]="s__") 
  
filter_and_score() {
    local metadata="$1" taxon="$2"
    awk -v target="$taxon" '
        BEGIN { FS = OFS = "\t"
                checkm2_completeness = 3 
                checkm2_contamination = 4
                checkm_completeness = 6
                checkm_contamination = 7
                gtdb_genome_representative = 18
                gtdb_representative = 19
                gtdb_taxonomy = 20 } 
        $gtdb_representative == "t" && index($gtdb_taxonomy, target) {
            score = ($checkm2_completeness - 5*$checkm2_contamination) + ($checkm_completeness - 5*$checkm_contamination)
            print $gtdb_genome_representative, $gtdb_taxonomy, score, $checkm2_completeness
        }' "$metadata"
}

dedup_by_taxon() {
    local rank_prefix="$1"
    awk -v rank="$rank_prefix" '
        BEGIN { FS = OFS = "\t" }
        match($2, rank "[^;]*") {
            taxon_name = substr($2, RSTART + length(rank), RLENGTH - length(rank))
            if (!seen[taxon_name]++) print $0, rank taxon_name }'
}

copy_genomes() {
    local acc_list="$1" reps="$2" out="$3"
    while IFS= read -r accession; do
        [[ -z "$accession" ]] && continue
        local L1="${accession:0:3}"   # GCF / GCA
        local L2="${accession:4:3}"   # 数字部分前三位，如 965 / 000
        local src="${reps}/${L1}/${L2}/${accession}_genomic.fna.gz"
        if [[ -f "$src" ]]; then
            cp -v "$src" "${out}/"
        else
            echo "WARN: not found: $src" >&2
        fi
    done < "$acc_list"
}

# main

filter_and_score "$metadata" "$taxon" \
    | sort -nr -k3,3 \
    | dedup_by_taxon "${levels[$rank]}" \
    > "$rep_tsv"

grep -w -v "${taxon}" "$rep_tsv" > "${rep_tsv}.tmp" && mv "${rep_tsv}.tmp" "$rep_tsv"

cut -f1 "$rep_tsv" | sed -E 's/^(RS|GB)_//' > "$acc_list"
copy_genomes "$acc_list" "$genomic_reps" "$outpath"
