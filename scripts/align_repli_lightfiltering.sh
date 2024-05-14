#!/bin/bash

set -eo pipefail

BOWTIE_REFERENCE=$1
READS_PATH=$2
OUTPUT_PATH=$3

mkdir -p `dirname $OUTPUT_PATH`

rclone cat "sharepoint:Bio Data/$READS_PATH" | \
  gunzip -c | \
  bowtie2 --threads "${BOWTIE_THREADS:-12}" -x "$BOWTIE_REFERENCE" \
    -U - 2> >(tee ${OUTPUT_PATH%.bam}.log >&2) | \
  samtools view -b | \
  samtools sort -T /tmp -@ 8 | \
  samtools markdup -@ 8 --duplicate-count -r -s - $OUTPUT_PATH \
    2> ${OUTPUT_PATH%.bam}.markdup.log
samtools index -@ 12 $OUTPUT_PATH