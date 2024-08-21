#!/usr/bin/python

import numpy as np
import pandas as pd
import sys

chic_samples = pd.read_csv("chic/chic_samples.csv")
chic_samples = chic_samples.loc[(chic_samples.rejected != True).fillna(True)]

celltype = dict(nos="Germline", tj="Somatic")

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
  f1 = "chic/wholeFragments_chr/{}/{}.bw".format(
    celltype[chic_sample["driver"]], chic_sample["sample"]
  )
  f2 = "chic/wholeFragments_chr/{}/{}.bw".format(
    celltype[chic_sample["driver"]], treat_sample["sample"]
  )
  print(
    "Rscript chic/median_log2fc.R {} {} chic/wholeFragments_chr/{}/L2FC_{}_{}_Rep{}.bw".format(
      f1, f2,
      celltype[chic_sample["driver"]],
      chic_sample["group"],
      celltype[chic_sample["driver"]],
      chic_sample["rep"]
    )
  )