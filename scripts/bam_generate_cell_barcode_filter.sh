#!/bin/bash

metadata=$1
sample_name=$2
grep_string=$3


grep '"'$sample_name'_[ATGC]*-1",.*"'$grep_string'"' $metadata | \
  awk -F '["_]' '{print "CB:Z:" $3}'