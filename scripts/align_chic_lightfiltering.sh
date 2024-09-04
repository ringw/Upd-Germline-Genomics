#!/bin/bash

set -eo pipefail

BOWTIE_REFERENCE=$1
READS_1_PATH=$2
READS_2_PATH=$3
OUTPUT_PATH=$4

mkdir -p `dirname $OUTPUT_PATH`

OUTPUT_1=`mktemp -u --suffix=.fifo`
mkfifo $OUTPUT_1
OUTPUT_2=`mktemp -u --suffix=.fifo`
mkfifo $OUTPUT_2

(
    rclone cat --include="`basename $READS_1_PATH`" "sharepoint:Bio Data/`dirname $READS_1_PATH`" | gunzip -c
) > $OUTPUT_1 &
(
    rclone cat --include="`basename $READS_2_PATH`" "sharepoint:Bio Data/`dirname $READS_2_PATH`" | gunzip -c
) > $OUTPUT_2 &

TEMP_PATH=$(mktemp -p /tmp --dry-run $(echo $OUTPUT_PATH | tr '/' '.').XXXXXXXXXX)
bowtie2 --threads "${BOWTIE_THREADS:-12}" -x "$BOWTIE_REFERENCE" \
    --no-discordant -1 $OUTPUT_1 -2 $OUTPUT_2 \
    2> >(tee ${OUTPUT_PATH%.bam}.log >&2) | \
    samtools view -b -f 0x02 | \
    samtools fixmate -m - - | \
    samtools sort -T /tmp -@ 8 - | \
    samtools markdup -@ 12 --duplicate-count -r -s - $OUTPUT_PATH \
        2> ${OUTPUT_PATH%.bam}.markdup.log
wait
samtools index -@ 12 $OUTPUT_PATH

rm $OUTPUT_1 $OUTPUT_2
