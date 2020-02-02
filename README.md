# _talklr_ uncovers ligand-receptor mediated intercellular crosstalk 
talklr uses Kullback-Leibler divergence, a concept from inforamtion theory, to uncover interesting ligand-receptor interactions among multiple cell types in single cell RNA-seq data.  
This is an R package for talklr, so that you can integrate talklr analysis into your any existing scRNA-seq bioinformatic workflow.  
To install talklr R package, you need to first install the _devltools_ R package, then run:
```r
library(devtools)
install_github("yuliangwang/talklr")
library(talklr)
```
Now talklr is ready for use!  
Alternatively, you can use the talklr Shiny App that does not require any programming experience. talklr App is available at:  
https://yuliangwang.shinyapps.io/talklr/