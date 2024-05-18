# chr lengths from dmel r6.47
chr.lengths = c(
  `2L`=23513712, `2R`=25286936, `3L`=28110227, `3R`=32079331, `4`=1348131, X=23542271, Y=3667352
)

masked.lengths = c(chr.lengths, `2L_Histone_Repeat_Unit`=5061)

experiment.driver <- tribble(
  ~driver, ~celltype,
  "nos", "Germline",
  "tj", "Somatic"
)