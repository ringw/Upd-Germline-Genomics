#!/bin/bash

set -e

READS_1_PATH=$1
READS_2_PATH=$2
OUTPUT_PATH=$3

mkdir -p `dirname $OUTPUT_PATH`

OUTPUT_1=`mktemp -u --suffix=.fifo`
mkfifo $OUTPUT_1
OUTPUT_2=`mktemp -u --suffix=.fifo`
mkfifo $OUTPUT_2

(
    rclone cat "sharepoint:Bio Data/$READS_1_PATH" | gunzip -c
) > $OUTPUT_1 &
(
    rclone cat "sharepoint:Bio Data/$READS_2_PATH" | gunzip -c
) > $OUTPUT_2 &

TEMP_PATH=$(mktemp -p /tmp --dry-run $(echo $OUTPUT_PATH | tr '/' '.').XXXXXXXXXX)
bowtie2 --threads 12 -x ~/OneDrive/Documents/Dmel_r6.47_bowtie2/Dmel_r6.47 \
    --no-discordant -1 $OUTPUT_1 -2 $OUTPUT_2 \
    2> >(tee ${OUTPUT_PATH%.bam}.log >&2) | \
    samtools view -b -L scripts/drosophila_genomic_regions.tab -q 20 | \
    # samtools collate -O - $TEMP_PATH.collate | \
    samtools fixmate -m - - | \
    samtools sort -T /tmp -@ 8 - | \
    samtools markdup -@ 12 -r -s - $OUTPUT_PATH \
        2> ${OUTPUT_PATH%.bam}.markdup.log
    # Suppress sixth column (base quality).
    # samtools mpileup -B - | \
    # gzip -9 -c > $OUTPUT_PATH 
wait
samtools index $OUTPUT_PATH

rm $OUTPUT_1 $OUTPUT_2