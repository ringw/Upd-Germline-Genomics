# chr lengths from dmel r6.47
chr.lengths = c(
  `2L`=23513712, `2R`=25286936, `3L`=28110227, `3R`=32079331, `4`=1348131, X=23542271, Y=3667352
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