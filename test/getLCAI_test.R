install.packages(
  "G:\\Coding\\getLcai\\getLCAI_test\\getLCAI_1.0.0.zip",
  repos = NULL,
  type = "win.binary"
)
library(getLCAI)

# configure absolute path to test the code
exp_test_repo = "G:\\Coding\\getLcai\\2\\exp_GSE165843.txt"
pheno_test_repo = "G:\\Coding\\getLcai\\2\\GSE165843_phe.txt"

exp_test <- read.table(exp_test_repo, header = TRUE)
pheno_test <- read.table(pheno_test_repo, header = TRUE)

data(exp_test)
data(pheno_test)
#data(exp_example)
#data(pheno_example)

outlist = getlcai(
  exp = "G:\\Coding\\getLcai\\2\\exp_GSE165843.txt",
  pheno = "G:\\Coding\\getLcai\\2\\GSE165843_phe.txt",
  control = "shAMPKa",
  experimental = "shCTL",
  type = 'Array',
  plotPCA = TRUE
)
