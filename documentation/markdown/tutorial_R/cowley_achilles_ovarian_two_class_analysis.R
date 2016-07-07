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

load("../../data/shrna/achilles_screens_expanded_with_weights.eset")
fdat = fData(achilles_screens_expanded)
pheno = pData(achilles_screens_expanded)

achilles_screens_expanded

cancerTypes = unique(pheno[,c("cell_line", "type")])
cancerTypes$tissue = ifelse(cancerTypes$type == "ovarian", "ovarian", "other")
# Remove the types column, since it also exists in the phenoData annotations of the Achilles annotations
cancerTypes = cancerTypes[,-2]


testIds = c(7849, #PAX8
            6656 #SOX1
)

t(getDefaultVariableMap())

vars = getDefaultVariableMap()
vars$reagentId = "reagent"
# These values currently need to be specified, but are not used in the single time-point Achilles analysis
# The need for these will hopefully be removed in the case of single time-point analyses in a future update
vars$timeNum = "cell_doublings"
vars$timeGroup = "doubling_time_hrs"

results = simem(achilles_screens_expanded,
               geneIds=testIds,
               covariate="tissue",
               annotationsPerCellLine=cancerTypes, 
               analyzeReagents=TRUE,
               inverseVarianceWeights=TRUE,
               endPoint=TRUE,
               parallelNodes=1,
               variableMap=vars)

# siMEM run on 2 processor cores
results = simem(achilles_screens_expanded,
               geneIds=testIds,
               covariate="tissue",
               annotationsPerCellLine=cancerTypes, 
               analyzeReagents=TRUE,
               inverseVarianceWeights=TRUE,
               endPoint=TRUE,
               parallelNodes=2,
               variableMap=vars)

results$gene

reagent = results$reagent
reagent[order(reagent$symbol, reagent$reagent), ]
