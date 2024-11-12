# chr lengths from dmel r6.47
chr.lengths <- c(
  `2L`=23513712, `2R`=25286936, `3L`=28110227, `3R`=32079331, `4`=1348131, X=23542271, Y=3667352
)
scaffolds.length <- 6178042

chr.colors <- setNames(
  hcl(seq(0, 330, length.out=8)[-8] + 100, 100, 65),
  names(chr.lengths)
)

arm.colors <- c(
  "2L" = hcl(126, 100, 65),
  "2LC" = hcl(126, 80, 35),
  "2RC" = hcl(173, 80, 35),
  "2R" = hcl(173, 100, 65),
  "3L" = hcl(31, 100, 65),
  "3LC" = hcl(31, 80, 35),
  "3RC" = hcl(79, 80, 35),
  "3R" = hcl(79, 100, 65),
  "4" = hcl(297, 100, 65),
  "X" = hcl(344, 100, 65),
  "Y" = hcl(250, 100, 65)
)

misc.lengths <- c(
  # `2Cen_mapped_Scaffold_10_D1684` = 19956, `2Cen_mapped_Scaffold_43_D1668` = 44411, 
  # `3Cen_mapped_Scaffold_1_D1896_D1895` = 76224, `3Cen_mapped_Scaffold_27_D1777` = 11983, 
  # `3Cen_mapped_Scaffold_31_D1643_D1653_D1791` = 87365, `3Cen_mapped_Scaffold_36_D1605` = 36913, 
  # `3Cen_mapped_Scaffold_41_D1641` = 22604, `3Cen_mapped_Scaffold_50_D1686` = 23238, 
  rDNA = 76973
)

masked.lengths = c(chr.lengths, `2L_Histone_Repeat_Unit`=5061)

masked.feature.lengths <- c(chr.lengths, misc.lengths, `2L_Histone_Repeat_Unit`=5061)

experiment.driver <- tribble(
  ~driver, ~celltype,
  "nos", "Germline",
  "tj", "Somatic"
)

repli_level_colors <- list(
  E="#2FEBF0", # hcl(196, 65, 85)
  EM="#AC7CDA", # hcl(284, 72, 60)
  ML="#9F3B5C", # hcl(356, 65, 40)
  L="#650812" # hcl(10, 59, 20)
)

repli_early_late_background <- list(
  E="#75E5EA", # hcl(196, 50, 85)
  L="#F2A2B7" # hcl(356, 50, 75)
)
