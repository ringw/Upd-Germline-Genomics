set -e

REFERENCE_FILE=references/chic_reference.fa.fai
EUCHROMATIN_REFS=$'2L\n2R\n3L\n3R\n4\nX\nY'
# Include multi-mapping reads from the locus (Histone locus).
SUPPLEMENTAL_LOCUS_2L="2L:21398000-21549000"

# Read the BAM file with our decided-upen filtering (0x2 mate mapped in proper
# pair, and the alignment of the fragment is the optimal alignment with a MAPQ
# >= 20).
function read_ref_with_mapq() {
  BAM_FILE=$1
  REF_NAME=$2
  samtools view -f 0x2 -F 0x100 -q 20 $BAM_FILE $REF_NAME
}

# For each fragment on the reference sequence, generate a new pair of reads
# extended all the way to the ends of the reference sequence. The seq and
# base-pair quality are reset to fake values, and the read length is kept at 50.
# This creates a nice small BAM file even though there may be a great deal of
# multi-mapping reads in the alignment. This is being done for short reference
# sequences (scaffolds, TEs) where we are later going to summarize the number of
# fragments aligned to it, and there is a performance gain to simplifying the
# contents of the BAM file as part of our sequence alignment target.
function read_short_scaffold() {
  BAM_FILE=$1
  REF_NAME=$2
  REF_LENGTH=$3
  READ_LENGTH=50
  READ_REVERSE_POS=$(( REF_LENGTH - READ_LENGTH + 1 ))
  READ_CIGAR='"'${READ_LENGTH}M'"'
  READ_PNEXT=1
  SEQ='"'`printf %${READ_LENGTH}s | tr " " "N"`'"'
  QUAL='"'`printf %${READ_LENGTH}s | tr " " "F"`'"'

  samtools view -f 0x2 -F 0x10 $BAM_FILE $REF_NAME | awk '
    {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, 1, $5, '$READ_CIGAR\
', $7, '$READ_REVERSE_POS', '$REF_LENGTH', '$SEQ', '$QUAL'}'
  samtools view -f 0x12 $BAM_FILE $REF_NAME | awk '
    {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, '$READ_REVERSE_POS', $5, '$READ_CIGAR\
', $7, '$READ_PNEXT', '-$REF_LENGTH', '$SEQ', '$QUAL'}'
}

(
  # Print SAM header to the new file.
  samtools view -H $1
  while read ind; do
    REF_NAME=`echo "$ind" | awk '{print $1}'`
    REF_LENGTH=`echo "$ind" | awk '{print $2}'`
    if [[ $REF_NAME == 2L ]]
    then
      # In subprocess: print header, print strictly mapping reads, and also
      # print all reads in the supplemental locus.
      # In pipe: Sort the two sets of reads, fix overlap between strictly
      # mapping reads/Histone locus reads using makrdup, and apply "view" (which
      # removes the SAM header).
      (
        samtools view -H $1
        read_ref_with_mapq $1 $REF_NAME
        samtools view $1 $SUPPLEMENTAL_LOCUS_2L
      ) | \
        samtools sort -T /tmp -@ 8 - | \
        samtools markdup -@ 12 -r -s -O SAM - - | \
        samtools view -O SAM
    elif echo "$EUCHROMATIN_REFS" | grep -w -q $REF_NAME
    then
      read_ref_with_mapq $1 $REF_NAME
    else
      read_short_scaffold $1 $REF_NAME $REF_LENGTH
    fi
  done < $REFERENCE_FILE
) | samtools view -b -o $2

samtools index -@ 12 $2
