#!/bin/bash

set -e
set -x

TJ_FILE_PATTERNS=(
  KZ2357/GC3772016_S1_L002
  KZ2357/GC3772017_S2_L002
  KZ2357/GC3772018_S3_L002
  KZ2504_KZ2505/GC3768010_S10_L001
  KZ2504_KZ2505/GC3768011_S11_L001
  KZ2504_KZ2505/GC3768012_S12_L001
  KZ2504_KZ2505/GC3768013_S13_L001
  KZ2504_KZ2505/GC3768014_S14_L001
  KZ2613_KZ1614/GC4651004_S4_L001
  KZ2613_KZ1614/GC4651005_S5_L001
  KZ2613_KZ1614/GC4651006_S6_L001
)

TJ_FILE_RCLONES=(
  "${TJ_FILE_PATTERNS[@]/#/sharepoint:Bio Data/}"
)

SAMPLE_RATE=0.05

OUTPUT_1=`mktemp -u --suffix=.fifo`
mkfifo $OUTPUT_1
OUTPUT_2=`mktemp -u --suffix=.fifo`
mkfifo $OUTPUT_2

for (( i=0; i<${#TJ_FILE_PATTERNS[@]}; i++ )); do
    f=${TJ_FILE_PATTERNS[i]}
    (
        rclone cat "sharepoint:Bio Data/$f"_R1_001.fastq.gz | gunzip -c | \
            perl -ne 'BEGIN{srand('$i')} $p = rand() < '$SAMPLE_RATE' if (/\@/); print if ($p)'
    ) > $OUTPUT_1 &
    (
        rclone cat "sharepoint:Bio Data/$f"_R2_001.fastq.gz | gunzip -c | \
            perl -ne 'BEGIN{srand('$i')} $p = rand() < '$SAMPLE_RATE' if (/\@/); print if ($p)'
    ) > $OUTPUT_2 &

    TEMP_PATH=chic_tj_sample_reads.`echo -n $f | tr / .`
    bowtie2 --threads 8 -x ~/OneDrive/Documents/Dmel_r6.47_bowtie2/Dmel_r6.47 \
        -1 $OUTPUT_1 -2 $OUTPUT_2 \
        2> >(tee $TEMP_PATH.log >&2) | \
        samtools view -b -L scripts/drosophila_genomic_regions.tab | \
        samtools collate -O - /tmp/$TEMP_PATH.collate | \
        samtools fixmate -m - - | \
        samtools sort -T /tmp/$TEMP_PATH -@ 8 - | \
        samtools markdup -@ 8 -r -s - $TEMP_PATH.markdup.bam 2> $TEMP_PATH.markdup.log
    wait
done

rm $OUTPUT_1 $OUTPUT_2

for f in ${TJ_FILE_PATTERNS[@]}; do
    echo -n chic_tj_sample_reads.
    echo -n $f | tr '/' '.'
    echo .markdup.bam
done | xargs samtools merge -f chic_tj_sample_reads.markdup.bam

samtools index chic_tj_sample_reads.markdup.bam
# rm
echo "${TEMP_BAMS[@]}"