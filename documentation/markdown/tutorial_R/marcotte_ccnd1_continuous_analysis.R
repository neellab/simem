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
#breast_screens

hp = read.delim("../../data/annotations/hairpin_annotations.txt", header=T, as.is=T, check.names=F)
hpWeights = hp[,c("trcn_id", "gene_id", "weight")]

rnaseq = read.delim("../../data/rnaseq/breast_rnaseq_qn.txt", header=T, as.is=T, check.names=F)
# Get the row index for CCND1 expression
ccnd1Index = which(rnaseq$symbol == "CCND1")
# Remove all identifier columns, and drop the columns dimension, resulting in a cell-line-named vector of log2-FPKM values
ccnd1Expr = unlist(rnaseq[ccnd1Index,-grep("gene_id|symbol|ensembl_id", colnames(rnaseq))])
ccnd1Expr

ccnd1Values = cbind.data.frame(cell_line=names(ccnd1Expr),
                               ccnd1=ccnd1Expr,
                               stringsAsFactors=FALSE)
ccnd1Values[1:5,]

genesOfInterest = c(595,   #CCND1
                    1019,  #CDK4
                    1021)  #CDK6

results = simem(screens = breast_screens,
                geneIds = genesOfInterest,
                covariate = "ccnd1",
                reagentWeights = hpWeights,
                annotationsPerCellLine = ccnd1Values,
                inverseVarianceWeights = TRUE,
                signalProbWeights = TRUE,
                analyzeReagents = TRUE,
                parallelNodes = 1
                )

results$gene

reagent = results$reagent
reagent[order(reagent$symbol, reagent$trcn_id), ]
