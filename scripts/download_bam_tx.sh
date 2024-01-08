#!/bin/bash

REMOTE_FOLDER="sharepoint:Chen Lab Backups/Upd_sc"

metadata=$1
condition_lower="${2%.*}"
batch_number="${2#*.}"
grep_string=$3
output_path=$4

cb_file=`mktemp`
bash "${0%download_bam_tx.sh}"bam_generate_cell_barcode_filter.sh $1 $2 $3 > $cb_file

if [[ $condition_lower == nos ]]; then
    dir_name=Nos-Upd_H3-GFP_Rep${batch_number}
else
    dir_name=Tj-Upd_H3-GFP_Rep${batch_number}
fi

rclone cat "$REMOTE_FOLDER/$dir_name/outs/possorted_genome_bam.bam" | \
    samtools view | grep -F xf:i:25 | grep -F -f $cb_file | \
    perl -ne 'for my $m (/(?:TX:Z:|\G;)(FBtr\d+)[^;]+/g) { print "$m\n" }' | \
    sort | uniq -c > $output_path