# !!! 2 steps are needed to test your own data.
# !!! You can also open getLCAI_test_example.R to see built-in results.
# !!! please make sure you have the right configuration!

# step 1: configure your own testing directory
setwd("E:/Coding/getLCAI/test")

# step 2: configure 5 necessary values
# relative exp path and pheno path
exp_test_path = "./test_data/exp_GSE165843.txt"
pheno_test_path = "./test_data/GSE165843_phe.txt"
control_type = "shAMPKa"
experimental_type = "shCTL"
data_type = "Array" # data_type can only be 'Array' or 'RNA-seq'

install.packages("./getLCAI_1.0.0.zip",
                 repos = NULL,
                 type = "win.binary")
library(getLCAI)

exp_test <- read.table(exp_test_path, header = TRUE)
pheno_test <- read.table(pheno_test_path, header = TRUE)

outlist = getlcai(
  exp = exp_test_path,
  pheno = pheno_test_path,
  control = control_type,
  experimental = experimental_type,
  type = data_type,
  plotPCA = TRUE
)

