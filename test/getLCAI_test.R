# !!! 2 steps are needed to test your data
# !!! please make sure you have the right configuration

# step 1: configure your own testing directory
setwd("E:/Coding/getLCAI/test")

# step 2: configure 5 necessary values
# Your won exp path and pheno path, relatively
exp_test_path = "./2/exp_GSE165843.txt"
pheno_test_path = "./2/GSE165843_phe.txt"
control_type = "Control"
experimental_type = "experimental"
data_type = "Array" # data_type can only be 'Array' or 'RNA-seq'

install.packages("./getLCAI_1.0.0.zip",
                 repos = NULL,
                 type = "win.binary")
library(getLCAI)

exp_test <- read.table(exp_test_path, header = TRUE)
pheno_test <- read.table(pheno_test_path, header = TRUE)

#data(exp_example)
#data(pheno_example)

outlist = getlcai(
  exp = exp_test_path,
  pheno = pheno_test_path,
  control = control_type,
  experimental = experimental_type,
  type = data_type,
  plotPCA = TRUE
)

