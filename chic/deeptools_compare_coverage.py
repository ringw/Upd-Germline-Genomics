#!/usr/bin/python

import numpy as np
import pandas as pd
import sys

chic_samples = pd.read_csv("chic/chic_samples.csv")
chic_samples = chic_samples.loc[(chic_samples.rejected != True).fillna(True)]

for unused_index, chic_sample in chic_samples.iterrows():
  if chic_sample["molecule"] != "H3":
    continue
  
  loctreat, = np.where(
    (chic_samples.group == chic_sample["group"])
    & (chic_samples.driver == chic_sample["driver"])
    & (chic_samples.rep == chic_sample["rep"])
    & (chic_samples.molecule.str.match("H3K"))
  )[0]
  treat_sample = chic_samples.iloc[loctreat]
  f1 = "chic/deepTools_masked/{}/{}.bw".format(
    treat_sample["group"], treat_sample["sample"]
  )
  f2 = "chic/deepTools_masked/{}/{}.bw".format(
    chic_sample["group"], chic_sample["sample"]
  )
  print(
    "bigwigCompare -b1 {} -b2 {} -o chic/deepTools_masked/{}/L2FC_{}_Rep{}.bw".format(
      f1, f2,
      chic_sample["group"],
      dict(nos="Germline", tj="Somatic")[chic_sample["driver"]],
      chic_sample["rep"]
    )
  )