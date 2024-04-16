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

cb_file=`mktemp`
bash "${0%estimate_library_size_bam_tx.sh}"metadata_to_samtools_tag_values_10x.sh $1 $2 $3 > $cb_file

if [[ $condition_lower == nos ]]; then
    dir_name=Nos-Upd_H3-GFP_Rep${batch_number}
else
    dir_name=Tj-Upd_H3-GFP_Rep${batch_number}
fi

rclone cat "$REMOTE_FOLDER/$dir_name/outs/possorted_genome_bam.bam" | \
    samtools view -@ 8 -D CB:$cb_file | grep -F xf:i:25 | grep -c -F RE:A:E