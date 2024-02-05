#!/bin/bash

REMOTE_FOLDER="sharepoint:Chen Lab Backups/Upd_sc"

# Path to scRNA-seq-Metadata.csv
metadata=$1
# condition/driver: Nos or Tj. However, the argument is expected as lowercase.
condition_lower="${2%.*}"
# Batch: Replicate of the Nos-driven or Tj-driven experiment (numbered 1 or 2).
batch_number="${2#*.}"
# Grep the csv file for this string, then get the record containing the cell
# barcode and print that to a file for searching SAM output.
grep_string=$3
# Output TSV: transcript id and pseudobulk count.
output_path=$4

cb_file=`mktemp`
bash "${0%download_bam_tx.sh}"bam_generate_cell_barcode_filter.sh $1 $2 $3 > $cb_file

if [[ $condition_lower == nos ]]; then
    dir_name=Nos-Upd_H3-GFP_Rep${batch_number}
else
    dir_name=Tj-Upd_H3-GFP_Rep${batch_number}
fi

mkdir -p `dirname "$output_path"`

rclone cat "$REMOTE_FOLDER/$dir_name/outs/possorted_genome_bam.bam" | \
    samtools view | grep -F xf:i:25 | grep -F -f $cb_file | \
    # The match follows a TX:Z: tag name, or a previous match \G followed by ;.
    perl -ne 'for my $m (/(?:TX:Z:|\G;)(FBtr\d+)[^;]+/g) { print "$m\n" }' | \
    sort | uniq -c > $output_path