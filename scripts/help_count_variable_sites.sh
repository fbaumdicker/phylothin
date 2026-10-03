#!/bin/bash
# bash commands to count variable sites in a panX SNP alignment.
# bash help_count_variable_sites.sh input_aln
# Copyright (C) 2026 Hannah Goetsch

input_aln="$1"

aln_file2="./SNP.aln.2"
tr -d '\n' < $input_aln > $aln_file2

aln_file3="./SNP.aln.3"
sed 's/ <unknown description>/,/g' $aln_file2 > $aln_file3
rm $aln_file2

aln_file1="./SNP.aln"
#sed 's/>/\\n/g' $aln_file3 > $aln_file1 # for GNU sed
sed 's/>/\
/g' $aln_file3 > $aln_file1 # for macOS/BSD sed
rm $aln_file3

# read the first non-empty record.
while IFS=',' read -r sequence_name sequence; do
    if [[ -n "$sequence" ]]; then
        break
    fi
done < "$aln_file1"

# remove a possible carriage return from Windows-formatted input.
sequence=${sequence//$'\r'/}

# count the characters in the sequence.
num_snps=${#sequence}

echo "Number of variable sites: $num_snps"
