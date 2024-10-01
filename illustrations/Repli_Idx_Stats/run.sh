#!/bin/bash

mkdir -p illustrations/Repli_Idx_Stats/bam
for f in nos_E2 nos_E4 nos_G2 nos_G4 nos_J2 nos_J4 nos_L2 nos_L4 tj_E1 tj_E2 \
    tj_G2 tj_J2 tj_L1 tj_L2; do
for f in nos_L2 nos_L4; do 
  samtools view -bq 20 -o illustrations/Repli_Idx_Stats/bam/$f.bam \
      repli/chr/$f.bam
  samtools index -@ 8 illustrations/Repli_Idx_Stats/bam/$f.bam
done
