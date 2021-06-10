# FR and CDR bounds in IMGT alignment
# http://www.imgt.org/IMGTrepertoire/2D-3Dstruct/FClength/human/IGH/IGHV/Hu_IGHVnber.html
# (Accessed Dec. 15 2019)
# Converted to nucleotide numbering
imgt_region_bounds <- tibble(
  region = c('FR1','CDR1','FR2','CDR2','FR3'),
  start = c(1, 79, 115,166,196),
  end = c(78,114,165,195,312)
)
# CDRs found by looking for conserved cysteine and tryp. positions (provided by partis)