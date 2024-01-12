rclone cat --include='Tj_Early-Mid-*.fastq.gz' sharepoint:"Chen Lab Backups"/Upd_Repli/Validation | \
  gunzip -c | \
  perl -ne 'BEGIN{srand(1)} $p = rand() < 0.05 if (/\@/); print if ($p)' | \
  bowtie2 --threads 8 -x ~/OneDrive/Documents/Dmel_r6.47_bowtie2/Dmel_r6.47 -U - 2> \
  >(tee repli_tj_sample_reads.log >&2) | \
  samtools sort - | \
  samtools markdup -@ 8 -r -s - repli_tj_sample_reads.markdup.bam 2> repli_tj_sample_reads.markdup.log
samtools index repli_tj_sample_reads.markdup.bam