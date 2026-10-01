#!/bin/bash
# bash commands to make phylogenetic trees ultrametric with treetime
# bash help_treetime.sh input_tree
# by Hannah Goetsch


# check dependency

if ! command -v treetime >/dev/null 2>&1; then
    error "TreeTime is not installed or is not available.
Install it with:
    pip install phylo-treetime"
fi


# default values

input_tree="$1" # given through command line
sequence_length=1000000 # only influences "resolution" (overall scaling) not actual branch length
clock_rate=0.00001 # only influences tree height
# ignore outliers: --clock-filter 0
outdir="./treetime_output/"
output_tree="./treetime_tree.nwk"

mkdir $outdir


dates_file="${outdir%/}/strains.csv"

# produce dates_file
{
    echo "strain,date"

    tr -d '\n\r' < "$input_tree" |
        grep -oE '[(,][^():,;]+' |
        sed -E 's/^[,(][[:space:]]*//; s/[[:space:]]*$//' |
        awk 'NF { print $0 ",2026" }' # same year for all strains required to ensure ultrametric tree
} > $dates_file

sed -i 's/"//g' $dates_file

# confirm that dates_file  has been produced

[[ -s "$dates_file" ]] ||
    error "Something went wrong in producing the strain-date-file."

# run treetime https://treetime.readthedocs.io/en/latest/index.html

treetime --tree $input_tree --dates $dates_file --sequence-length $sequence_length --clock-rate $clock_rate --keep-root --clock-filter 0 --outdir $outdir

nexus_tree="${outdir%/}/timetree.nexus"

# confirm that treetime produced the expected Nexus tree

[[ -s "$nexus_tree" ]] ||
    error "TreeTime completed, but not successfull. Check the error message above."

echo "The Nexus ultrametric tree can be found here: $nexus_tree ."


# convert nexus-output of treetime to nwk-file with  R (ape)

Rscript --vanilla - "$nexus_tree" "$output_tree" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
nexus_file  <- args[1]
newick_file <- args[2]
tree <- ape::read.nexus(nexus_file)
ape::write.tree(tree, file = newick_file)
RSCRIPT

echo "TreeTime done. Ultrametric tree: $output_tree ."
