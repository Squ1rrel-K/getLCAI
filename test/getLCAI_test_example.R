# install
setwd("E:/Coding/getLCAI/test")
# install.packages("./getLCAI_1.0.0.zip", repos = NULL, type = "source")      # Linux
install.packages("./getLCAI_1.0.0.zip", repos = NULL, type = "win.binary") 
library(getLCAI)
data(exp_example)
data(pheno_example)

outlist = getlcai(exp = exp_example,
                  pheno = pheno_example,
                  control = "Control",
                  experimental = "experimental",
                  type = 'Array',
                  plotPCA = TRUE)

