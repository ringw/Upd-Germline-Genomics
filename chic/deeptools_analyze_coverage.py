#!/usr/bin/python

import numpy as np
import pandas as pd
import sys

chic_samples = pd.read_csv("chic/chic_samples.csv")
chic_samples = chic_samples.loc[(chic_samples.rejected != True).fillna(True)]

for unused_index, chic_sample in chic_samples.iterrows():
  if "ChIP" in chic_sample["group"]:
    sys.stdout.write(
      "bamCoverage --extendReads 200 --minMappingQuality 20 --normalizeUsing RPKM -bs 20"
    )
  else:
    sys.stdout.write(
      "bamCoverage --minFragmentLength 150 --maxFragmentLength 200 --extendReads 200 --minMappingQuality 20 --normalizeUsing RPKM -bs 20"
    )
  sys.stdout.write(
    " -b chic/masked/{}/{}.bam".format(chic_sample["group"], chic_sample["sample"])
  )
  sys.stdout.write(
    " -o chic/deepTools_masked/{}/{}.bw\n".format(chic_sample["group"], chic_sample["sample"])
  )