#!/bin/bash

metadata=$1
sample_name=$2
grep_string=$3


grep '"'$sample_name'_[ATGC]*-1",.*"'$grep_string'"' $metadata | \
  # awk -F '["_]' '{print "xf:i:25.*CB:Z:" $3; print "CB:Z:" $3 ".*xf:i:25"}'
  awk -F '["_]' '{print "CB:Z:" $3}'