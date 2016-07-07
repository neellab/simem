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


load("../../data/shrna/breast_screens_with_weights.eset")
breast_screens

hp = read.delim("../../data/annotations/hairpin_annotations.txt", header=T, as.is=T, check.names=F)
hpWeights = hp[,c("trcn_id", "gene_id", "weight")]

subtypes = read.delim("../../data/annotations/cell_line_subtypes.txt", header=T, as.is=T, check.names=F)
status = subtypes[,c("cell_line", "subtype_neve")]

status$erbb2 = ifelse(status$subtype_neve == "her2", "her2", "other")
status[1:10,]

signalingPathway = c(207,   #AKT1
                    11140, #CDC37
                    2064,  #ERBB2
                    55914, #ERBB2IP
                    2065,  #ERBB3
                    2475,  #MTOR
                    5290,  #PIK3CA
                    6009,  #RHEB
                    25803, #SPDEF
                    7022) #TFAP2C

results = simem(screens = breast_screens,
                geneIds = signalingPathway,
                covariate = "erbb2",
                reagentWeights = hpWeights,
                annotationsPerCellLine = status,
                inverseVarianceWeights = TRUE,
                signalProbWeights = TRUE,
                analyzeReagents = TRUE,
                covariateFactorOrder = c("other", "her2"), 
                parallelNodes = 1
                )

results = simem(screens = breast_screens,
                geneIds = signalingPathway,
                covariate = "erbb2",
                reagentWeights = hpWeights,
                annotationsPerCellLine = status,
                inverseVarianceWeights = TRUE,
                signalProbWeights = TRUE,
                analyzeReagents = TRUE,
                covariateFactorOrder = c("other", "her2"), 
                parallelNodes = 3
                )

results$gene

str(results$gene_detailed)

reagent = results$reagent
reagent[order(reagent$symbol, reagent$trcn_id), ]

str(results$reagent_detailed)
