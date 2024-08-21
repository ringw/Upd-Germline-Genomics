#!/usr/bin/python

import numpy as np
import pandas as pd
import sys

chic_samples = pd.read_csv("chic/chic_samples.csv")
chic_samples = chic_samples.loc[(chic_samples.rejected != True).fillna(True)]

celltype = dict(nos="Germline", tj="Somatic")

for unused_index, chic_sample in chic_samples.iterrows():
  if "ChIP" in chic_sample["group"] or "GC4651004_S4_L001" == chic_sample["sample"]:
    continue
  sys.stdout.write(
    "Rscript chic/median_coverage.R"
  )
  sys.stdout.write(
    " {} chic/wholeFragments_chr/{}/{}.bw\n".format(
      chic_sample["sample"],
      celltype[chic_sample["driver"]],
      chic_sample["sample"]
    )
  )