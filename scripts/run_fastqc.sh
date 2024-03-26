#!/bin/bash

set -e
set -x # undo this later

TEMP_DIR=`mktemp -d`
# ##*/ is the basename operator
TEMP_FILE=$TEMP_DIR/${1##*/}

rclone cat "sharepoint:Bio Data/$1" > $TEMP_FILE
fastqc $TEMP_FILE
cp ${TEMP_FILE%.fastq.gz}_fastqc.zip $2
rm -r $TEMP_DIR