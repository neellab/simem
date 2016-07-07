rm(list=ls())

### To install required packages, uncomment and run this
# source("http://www.bioconductor.org/biocLite.R")
# biocLite(c("Biobase", "preprocessCore", "genefilter"))
# install.packages(c("blme", "doParallel", "ggplot2", "locfit", "MASS", "plyr", "reshape"))

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("blme"))
suppressPackageStartupMessages(library("doParallel"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("locfit"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("preprocessCore"))
suppressPackageStartupMessages(library("reshape"))


source("../../R/data_format_lib.R")
source("../../R/model_lib.R")
source("../../R/simem_lib.R")

load("../../data/shrna/achilles_screens_with_weights.eset")
fdat = fData(achilles_screens)
pheno = pData(achilles_screens)

achilles_screens

hp = read.delim("../../data/annotations/achilles_hairpin_annotations.txt", header=T, as.is=T, check.names=F)
hpWeights = hp[,c("probeset_id", "gene_id", "weight")]

expr = read.delim("../../data/expression/ccle_achilles_expression.txt", header=T, as.is=T, check.names=F)

# Perform the analysis for ids contained in both the expression matrix and the shRNA data.
commonIds = intersect(expr$gene_id, fdat$gene_id)

# Subset the shRNA data to only include the Achilles screens profiled by CCLE expressions
samplesToInclude = which(pheno$cell_line %in% colnames(expr))
achilles = achilles_screens[,samplesToInclude]

testIds = c(598, #BCL2L1
            3845, #KRAS
            7849, #PAX8
            6663 #SOX10
)

t(getDefaultVariableMap())

vars = getDefaultVariableMap()
vars$reagentId = "probeset_id"
vars$timeNum = "timeNum"
vars$timeGroup = "timepointIndex"

results = simem(achilles,
               geneIds=testIds,
               reagentWeights=hpWeights,
               annotationsPerId=expr, 
               analyzeReagents=TRUE,
               inverseVarianceWeights=TRUE,
               endPoint=TRUE,
               parallelNodes=1,
               variableMap=vars)

results = simem(achilles,
               geneIds=testIds,
               reagentWeights=hpWeights,
               annotationsPerId=expr, 
               analyzeReagents=TRUE,
               inverseVarianceWeights=TRUE,
               endPoint=TRUE,
               parallelNodes=2,
               variableMap=vars)

results$gene

reagent = results$reagent
reagent[order(reagent$symbol, reagent$probeset_id), ]
