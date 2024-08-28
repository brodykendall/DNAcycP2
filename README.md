#TODO:

# TEMPORARY README FILE

# Usage:

```r
devtools::install_github("brodykendall/dnacycpv2-R")
# Substitute path_to_data:
dat_tiling <- read.csv("/path_to_data/cycle5.txt")
dat_tiling_pred <- dnacycpv2::cycle(dat_tiling$Sequence)
```
