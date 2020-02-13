# _talklr_ uncovers ligand-receptor mediated intercellular crosstalk 
talklr (intercellular cross__TALK__ using __L__igand-__R__eceptor)uses Kullback-Leibler divergence, a concept from inforamtion theory, to uncover interesting ligand-receptor interactions among multiple cell types in single cell RNA-seq data.  
This is an R package for talklr, so that you can integrate talklr analysis into your any existing scRNA-seq bioinformatic workflow.  
To install talklr R package, you need to first install the _devltools_ R package, then run:
```r
library(devtools)
install_github("yuliangwang/talklr")
library(talklr)
```
You can follow the tutorial below to get started using talklr on your single cell RNA-seq data (or bulk RNA-seq of purified cells):  
https://rpubs.com/wang341/570900  
Now talklr is ready for use!  
Alternatively, you can use the talklr Shiny App that does not require any programming experience. talklr App is available at:  
https://yuliangwang.shinyapps.io/talklr/  
If you use talklr in your work, please cite the paper below:  
Yuliang Wang. (2020). talklr uncovers ligand-receptor mediated intercellular crosstalk. bioRxiv, 2020.02.01.930602. doi:10.1101/2020.02.01.930602
