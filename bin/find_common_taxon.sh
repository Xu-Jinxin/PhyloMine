#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 [ -p path ] [ -s|--strict ]"
    echo "  -p, --path       GTDB result path"
    echo "  -s, --strict     Suffix sensitive: treat o__Flavobacteriales_D and o__Flavobacteriales as different"
    echo "  -h, --help       Show this help"
    exit 1
}

path=""
strict=0

OPTS=$(getopt -o p:sh -l path:,strict,help -n "$0" -- "$@") || usage
eval set -- "$OPTS"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--path)   path="$2";  shift 2 ;;
        -s|--strict) strict=1;   shift 1 ;;
        -h|--help)   usage ;;
        --) shift; break ;;
        *)  echo "Internal error: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$path" ]]; then
    echo "Error: -p/--path is required" >&2
    usage
fi

bac="${path}/gtdbtk.bac120.summary.tsv"
ar="${path}/gtdbtk.ar53.summary.tsv"

if [[ -s "$bac" && -s "$ar" ]]; then
    echo "Error: both bac120 and ar53 summaries exist in $path" >&2
    echo "       PhyloMine cannot handle genomes spanning multiple domains in a single run." >&2
    exit 1
elif [[ -s "$bac" ]]; then
    gtdb_summary="$bac"
elif [[ -s "$ar" ]]; then
    gtdb_summary="$ar"
else
    echo "Error: neither gtdbtk.bac120.summary.tsv nor gtdbtk.ar53.summary.tsv found in $path" >&2
    exit 1
fi

#   strict=0 时去掉 "__" 后名字里第一个 "_" 起的后缀
#   o__Flavobacteriales_D -> o__Flavobacteriales
#   c__Bacteroidia        -> c__Bacteroidia（无后缀，不变）

cut -f2 "${gtdb_summary}" | sed '1d' | awk -F';' -v strict="$strict" '
function norm(s) {
    if (strict) return s
    pos = index(s, "__")
    if (pos == 0) return s
    prefix = substr(s, 1, pos + 1)
    rest   = substr(s, pos + 2)
    p = index(rest, "_")
    if (p == 0) return s
    return prefix substr(rest, 1, p - 1)
}
NR == 1 {
    for (i = 1; i <= NF; i++) common[i] = norm($i)
    n = NF
    next }
{
    for (i = 1; i <= n; i++)
        if (norm($i) != common[i]) {
            for (j = i; j <= n; j++) common[j] = ""
            break
        }}
END {
    if (n == "") exit
    for (i = n; i >= 1; i--)
        if (common[i] != "" && common[i] !~ /^[a-z]__$/) {
            print common[i]
            break
        }
}'