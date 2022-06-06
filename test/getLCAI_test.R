# configure your own working directory
setwd("F:/Code/getLCAI/test")

install.packages(
  "./getLCAI_1.0.0.zip",
  repos = NULL,
  type = "win.binary"
)
library(getLCAI)

# configure relative path to test the code
exp_test_path = "./2/exp_GSE165843.txt"
pheno_test_path = "./2/GSE165843_phe.txt"

exp_test <- read.table(exp_test_path)
pheno_test <- read.table(pheno_test_path)

#data(exp_example)
#data(pheno_example)

outlist = getlcai(
  exp = exp_test_path,
  pheno = pheno_test_path,
  control = "treated",
  experimental = "Control",
  type = 'Array',
  plotPCA = TRUE
)

