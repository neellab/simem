library("Biobase")
library("preprocessCore")
library("locfit")
library("blme")
library("reshape")
# Contains the empty() function
library("plyr")
library("genefilter")
library("ggplot2")
library("MASS")

# Cut out 3 legacy parameters
simemIsogenic = function(screens, 
                        covariate=NULL, 
                        geneIds=NULL, 
                        reagentWeights=NULL, 
                        cellLineWeights=NULL, 
                        annotationsPerCellLine=NULL, 
                        covariateFactorOrder=NULL,
                        inverseVarianceWeights=TRUE, 
                        signalProbWeights=FALSE, 
                        analyzeReagents=FALSE, 
                        calculateRelativeDropout=FALSE,
                        # default is "gaussian", with "t" producing p-values as reported in manuscript
                        pvalueType="gaussian",
                        # Fit siMEM models using the lme4 or blme packages as specified
                        modelEngine="blme",
                        variableMap=getDefaultVariableMap(),
                        formulas=NULL, 
                        printFreq=20, 
                        parallelNodes=1
			   ) {

	results = simemAnalysis(screens=screens,
                          covariate=covariate,
                          geneIds=geneIds,
                          analysisType="isogenic",
                          reagentWeights=reagentWeights,
                          cellLineWeights=cellLineWeights,
                          annotationsPerCellLine=annotationsPerCellLine,
                          annotationsPerId=NULL,
                          covariateFactorOrder=covariateFactorOrder,
                          inverseVarianceWeights=inverseVarianceWeights,
                          signalProbWeights=signalProbWeights,
                          signalProb=signalProb,
                          analyzeReagents=analyzeReagents,
                          calculateRelativeDropout=calculateRelativeDropout,
                          pvalueType=pvalueType,
                          modelEngine=modelEngine,
                          variableMap=variableMap,
                          formulas=formulas,
                          printFreq=printFreq,
                          endPoint=FALSE,
                          parallelNodes=parallelNodes)

	return(results)
}

simem = function(screens, 
                geneIds=NULL, 
                covariate=NULL, 
                reagentWeights=NULL, 
                cellLineWeights=NULL, 
                annotationsPerCellLine=NULL, 
                annotationsPerId=NULL, 
                covariateFactorOrder=NULL, 
                inverseVarianceWeights=TRUE, 
                signalProbWeights=FALSE, 
                analyzeReagents=FALSE,
                calculateRelativeDropout=FALSE,
                # default is "gaussian", with "t" producing p-values as reported in manuscript
                pvalueType="gaussian",
                # Fit siMEM models using the lme4 or blme packages as specified
                modelEngine="blme",
                variableMap=getDefaultVariableMap(), 
                formulas=NULL, 
                printFreq=20, 
                endPoint=FALSE,
                parallelNodes=1
				   ) {

	results = simemAnalysis(screens=screens,
          							geneIds=geneIds,
          							covariate=covariate,
          							analysisType="panel",
          							reagentWeights=reagentWeights,
          							cellLineWeights=cellLineWeights,
          							annotationsPerCellLine=annotationsPerCellLine,
          							annotationsPerId=annotationsPerId,
          							covariateFactorOrder=covariateFactorOrder,
          							inverseVarianceWeights=inverseVarianceWeights,
          							signalProbWeights=signalProbWeights,
          							analyzeReagents = analyzeReagents,	
                        calculateRelativeDropout=calculateRelativeDropout,
          							pvalueType=pvalueType,
          							modelEngine=modelEngine,
          							variableMap = variableMap,
          							formulas=formulas,
          							printFreq=printFreq,
          							endPoint=endPoint,
          							parallelNodes=parallelNodes
          							)

	return(results)
}


##
# The simemAnalysis() function does basic type and input sanity checking, then calls the simemChunk() function
# to perform the main data extraction and model fitting loop. In cases where more than one processor is specified
# this function calls the doMC functions to register multiple processor cores, split the gene ids into equal sized
# subsets, and calls simemChunk() functions in parallel with each gene id subset. Once all concurrent jobs complete,
# the function concatenates the fitted models from multiple processors (if necessary), and applies Benjamini-Hochberg
# FDR correction to model parameter P-values.
# 
##
simemAnalysis = function(screens, 
                        geneIds=NULL, 
                        covariate=NULL,
                        analysisType="panel", 
						            reagentWeights=NULL,
						            cellLineWeights=NULL, 
						            annotationsPerCellLine=NULL,
						            annotationsPerId=NULL,
						            covariateFactorOrder=NULL,
						            inverseVarianceWeights=TRUE,
						            signalProbWeights=FALSE,
						            analyzeReagents=TRUE,
                        calculateRelativeDropout=FALSE,
						            # default is "gaussian", with "t" producing p-values as reported in manuscript
						            pvalueType="gaussian",
						            # Fit siMEM models using the lme4 or blme packages as specified
						            modelEngine="blme",
            						variableMap=getDefaultVariableMap(),
#						            formulas=getModelFormulas(),
                        formulas=NULL, 
						            printFreq=20, 
						            endPoint=FALSE, 
						            parallelNodes=1
				   ) {

	if(parallelNodes > 1) {
		library(doParallel)
	  cl = makePSOCKcluster(parallelNodes)
	  registerDoParallel(cl)
#	  registerDoParallel(cores=parallelNodes)
	}

	##################################################################
	# BEGIN - Input variable sanity check
	##################################################################
	allowedEsets = c("ExpressionSet")

	if(length(intersect(allowedEsets, class(screens))) < 1) {
		expected = paste(allowedEsets, sep="", collapse=", ")
		observed = paste(class(screens), sep="", collapse=", ")

		msg = paste("ERROR in simemAnalysis():\n\n",
					"'screens' input paramter must be an object instantiating one of the following classes: ", expected,
					"\n\nHowever, the following classes were associated with the 'screens' paramter: ", observed,
					"\n\n\n\n", sep="")
		stop(msg)
	}

	if(!(analysisType == "panel" || analysisType == "isogenic")) {
		msg = paste("ERROR in simemAnalysis():\n\n",
					"Incorrect 'analysisType' parameter (either 'panel', 'panel_two_covariate' or 'isogenic' allowed): ", analysisType,
					"\n\n\n\n", sep="")
		stop(msg)
	}

	if(analysisType == "isogenic" && is.data.frame(annotationsPerId)) {
		msg = paste("ERROR in simemAnalysis():\n\n",
					"Specifying an 'annotationsPerId' parameter value is incompatible with an 'isogenic' analysis type.",
					"\nThis parameter is intended to specify per-gene annotations, such as copy number status, across screens.",
					"\nIt is thus meant to be used in conjunction with an across-screens (or 'panel') analysis type.",
					"\nPlease do not specify the 'annotationsPerId' parameter in the simemAnalysis() function call.", 
					"\n\n\n\n", sep="")
		stop(msg)
	}


	# Load screens data for sanity-checking
	pheno = pData(screens)
	fdat = fData(screens)

	# After we verify that the input is in ExpressionSet format, extract the gene ids if none specified
	if(!is.null(geneIds)) {
	  geneIds = unique(geneIds)
	} else {
	  fdat = fData(screens)
	  # Extract the geneId column's name from the variableMap
	  geneIdParam = variableMap$geneId
	  geneIds = unique(fdat[,geneIdParam])
	}
	
	
	##################################################################
	# BEGIN - Input variable sanity check
	##################################################################
	# TEST that key columns are present in screens phenoData
	phenoColsToKeep = c(variableMap$sampleId,
          						variableMap$cellLineId,
          						variableMap$replicateId,
          						variableMap$timeGroup,
          						variableMap$timeNum)

	if(length(intersect(phenoColsToKeep, colnames(pheno))) != length(phenoColsToKeep)) {
		expected = paste(phenoColsToKeep, sep="", collapse=", ")
		observed = paste(intersect(phenoColsToKeep, colnames(pheno)), sep="", collapse=", ")

		msg = paste("ERROR in simemAnalysis():\n\n",
					"screens phenoData data frame must have, at a minimum, the following columns: ", expected,
					"\n\nHowever, only the following were observed: ", observed,
					"\n\nRequired column names can be customized using the 'variableMap' parameter",
					"\n\n\n\n", sep="")
		stop(msg)
	}

	if(!is.data.frame(variableMap) || nrow(variableMap) != 1) {
		msg = paste("ERROR in simemAnalysis():\n\n",
    					"Incorrect 'variableMap' input (data frame with single row required)",
    					"\n\n\n\n", sep="")
		stop(msg)
	}

  # Ensure that if if reagentWeights is a data frame,
	# It contains the reagentId as its first column header, 
	# and it contains the weight column
	if(is.data.frame(reagentWeights) &&
	   (colnames(reagentWeights)[1] != variableMap$reagentId ||
	   length(grep("weight", colnames(reagentWeights))) == 0)) {
		msg = paste("ERROR in simemAnalysis():\n\n",
    					"Reagent weights table must have as first column: ", variableMap$reagentId,
    					"\n\nThe data frame must also contain a 'weight' column.",
    					"\n\nReagent identifier can be customized using the 'variableMap' parameter",
    					"\n\n\n\n", sep="")
		stop(msg)
	}

	# Ensure that the annotations per cell line, if specified, have the cell line as first 
	# column
	if(is.data.frame(annotationsPerCellLine) && 
	   colnames(annotationsPerCellLine)[1] != variableMap$cellLineId) {
		msg = paste("ERROR in simemAnalysis():\n\n",
					"Additional 'annotationsPerCellLine' data frame must have as first column: ", variableMap$cellLineId,
					"\n\nRequired column names can be customized using the 'variableMap' parameter",
					"\n\n\n\n", sep="")
		stop(msg)
	}

	# If we're specifying per gene annotations (eg: Copy Number status for each gene)
	# ensure that the first column of annotationsById is "gene_id"
	if(is.data.frame(annotationsPerId) && colnames(annotationsPerId)[1] != variableMap$geneId) {
		msg = paste("ERROR in simemAnalysis():\n\n",
					"Additional 'annotationsPerId' data frame must have as first column:  ", variableMap$geneId,
					"\n\nRequired column names can be customized using the 'variableMap' parameter",
					"\n\n\n\n", sep="")
		stop(msg)
	}

	# TEST sampleId labeling for the phenoData of the screens
	if(length(intersect(colnames(pheno), variableMap$sampleId)) != 1) {
		observed = paste(colnames(pheno), sep="", collapse=", ")
		msg = paste("ERROR in simemAnalysis():\n\n",
					"At least one of the screens phenoData columns must be: ", variableMap$sampleId,
					"\n\nWe observe the following columns in the phenoData: ", observed,		
					"\n\nRequired column names can be customized using the 'variableMap' parameter", "\n\n\n\n", sep="")
		stop(msg)
	}

	# TEST gene/reagent ids
	# test that gene and reagent ids are properly set in the fData of the screens
	idVars = c(variableMap$reagentId, variableMap$geneId, variableMap$symbol)
	idMatch = intersect(colnames(fdat), idVars)

	if(length(idMatch) != length(idVars)) {
		expected = paste(idVars, sep="", collapse=", ")
		observed = paste(colnames(fdat), sep="", collapse=", ")
		msg = paste("ERROR in simemAnalysis():\n\n",
					"Cannot find reagent id, gene id or gene symbol in screens featureData.",
					"\n\nWe expect to find the following columns: ", expected,
					"\n\nWe observe the following columns in the featureData: ", observed,
			"\n\nId variable columns can be customized using the 'variableMap' parameter", "\n\n\n\n", sep="")
		stop(msg)
	}

	#############################################
	# MERGE ADDITIONAL ANNOTATIONS into the pheno matrix if theye're provided
	# This merged data frame will be used to sanity-check other input paramters
	#############################################
	if(is.data.frame(annotationsPerCellLine)) {
		# add additional columns to the pheno on the fly
		pheno = merge(pheno, annotationsPerCellLine)
		# Merge may not preserve row names, and these are used later on
		rownames(pheno) = pheno[,variableMap$sampleId]
	}

	# Ensure covariate variable exists in either, provided it hasn't been extracted
	# from annotationsPerId just above
	if(is.character(covariate) && !is.data.frame(annotationsPerId) && !(length(intersect(colnames(pheno), covariate)) %in% c(1,2))) {
		msg = paste("ERROR in simemAnalysis():\n\n",
					"Cannot find the following 'covariate' column(s) in the screens phenoData (or annotationsPerCellLine, if specified): ", covariate,
					"\n\n\n\n", sep="")
		stop(msg)
	}

	# Test whether filter variables match those found in the data (pheno+annotations)
# 	if(is.list(filters)) {
# 		filterVars = names(filters)
# 		matchedVars = intersect(filterVars, colnames(pheno))
# 		mismatchedVars = setdiff(filterVars, colnames(pheno))
# 
# 		if(length(mismatchedVars) > 0) {
# 			expected = paste(matchedVars, sep="", collapse=", ")
# 			observed = paste(mismatchedVars, sep="", collapse=", ")
# 			msg = paste("ERROR in simemAnalysis():\n\n",
# 						"Some specified filter variables not found in phenoData or additional annotations.",
# 						"\n\nFilter variables matching existing phenoData variables: ", expected,
# 						"\n\nFilter variables missing from phenoData: ", observed, "\n\n\n\n", sep="")
# 			stop(msg)
# 		}
# 	}

	if(!is.null(covariate) && length(covariate) > 1 && is.data.frame(annotationsPerId)) {
		msg = paste("WARNING in simemAnalysis(): Both multiple 'covariate' values and 'annotationsPerId' parameters were specified.",
					"\n\nAs a result, 'covariate' parameter values beyond the first will be ignored, and the per gene id, per cell line covariate values extracted from annotationsPerId will be used as one of two allowed covariates.",
					"\n\n\n\n", sep="")
		warning(msg)
	}
	##################################################################
	# END - Input variable sanity check
	##################################################################

	# If we're performing an end-point analysis, adjust the time-course formulas used
	# by for marcotte data
	if(is.null(formulas) & endPoint == TRUE) {
		formulas = getModelFormulas()
		formulas$panelGene = as.formula("expr ~ 1 + COVARIATE + (1|reagentId/cellLineId)")
		formulas$panelRg = as.formula("expr ~ 1 + COVARIATE + (1|cellLineId)")
	} else if(is.null(formulas)) {
		formulas = getModelFormulas()
	}

	##################################################################
  # Check whether the weight matrices are included, if a weighted regression is specified
	##################################################################
	if(inverseVarianceWeights & is.null(assayData(screens)[["varWeights"]])) {
	  msg = paste("WARNING in simemAnalysis(): 'varWeights' entry not found in assayData(screens).",
	              "This occurs if the 'screens = addPrecisionWeights(screens)' command, which calculates mean-variance weights, has not been run",
                "Setting inverseVarianceWeights=FALSE for this analysis", 
	              "\n\n\n\n", sep="")
	  warning(msg)
	  inverseVarianceWeights = FALSE
		#screens = addPrecisionWeights(screens)
	}

	if(signalProbWeights & is.null(assayData(screens)[["signalProbWeights"]])) {
	  msg = paste("WARNING in simemAnalysis(): 'signalProbWeights' entry not found in assayData(screens).",
	              "This occurs if the 'screens = addSignalProbWeights(screens, signalProb)' command, which calculates signal-noise weights, has not been run",
	              "Setting signalProbWeights=FALSE for this analysis", 
	              "\n\n\n\n", sep="")
	  warning(msg)

	  signalProbWeights = FALSE
		#screens = addSignalProbWeights(screens, signalProb)
	}


	if(parallelNodes > 1) {

		# Transform gene id vector into a list with each entry containing a subset of the gene ids
		# to be passed to a different processing core.
		chunkSize = ceiling(length(geneIds)/parallelNodes)
		chunkIndex = rep(1:parallelNodes, each=chunkSize)
		chunkIndex = chunkIndex[1:length(geneIds)]
		## Split Ids into chunks and do parallel processing
		geneIdsList = by(geneIds, INDICES=list(chunkIndex), function(v) v)

		packages = c("Biobase", "preprocessCore", "blme", "reshape", "genefilter", "plyr")

# For whatever reason, the first solution doesn't work, so hard-code the vector for now
#		fns = sapply(ls(), function(element) is.function(eval(as.name(element))))
#    fnsToExport = names(fns[fns == TRUE])
#     The code below only passes the 3 hard-coded functions to each of the processes, and ignores the rest of the vector
#     Not hard-coding anything leads to no function definitions being passed to the parallel processes, resulting in errors
#    fnsToExport = c("simemChunk", "extractData", "addPrecisionWeights", as.character(fnsToExport))
    fnsToExport = c("addPrecisionWeights",
                    "addSignalProbWeights",
                    "applyFDR",
                    "assignWeightsToObservations",
                    "combineResults",
                    "extractData",
                    "fitAnova",
                    "fitIsogenicCovariateGene",
                    "fitIsogenicCovariateRg",
                    "fitModel",
                    "fitModelAndFormatOutput",
                    "fitPanelCovariateGene",
                    "fitPanelCovariateRg",
                    "fitSimplerRandomEffects",
                    "formatOutput",
                    "getDefaultVariableMap",
                    "getMeanSd",
                    "getMeanVar",
                    "getMinimalSummary",
                    "getModelFormulas",
                    "getModelPValues",
                    "getRelativeDropout",
                    "getRelativeDropoutRates",
                    "getReplicateMeanVar",
                    "getSignalProb",
                    "getVarianceWeights",
                    "loadMeanVar",
                    "pvaluesNorm",
                    "pvaluesT",
                    "resetFactorLevels",
                    "simem",
                    "simemAnalysis",
                    "simemChunk",
                    "simemIsogenic",
                    "standardizeColumnNames",
                    "summaryAnova",
                    "summaryFitEffectCoefs",
                    "summaryFitLogLik",
                    "tryCatchWarningsErrors",
                    "withinBetweenDF")
    
    # Split siMEM processing to multiple processor cores, with a subset of gene ids processed by each core.
		final = foreach(geneIds=iter(geneIdsList), .combine="combineResults", .packages=packages, .export=fnsToExport, .verbose=T) %dopar% {
					results <- simemChunk(screens=screens,
					                      geneIds=geneIds,
					                      analysisType=analysisType,
					                      covariate=covariate,
					                      reagentWeights=reagentWeights,
					                      cellLineWeights=cellLineWeights,
					                      annotationsPerCellLine=annotationsPerCellLine,
					                      annotationsPerId=annotationsPerId,
					                      covariateFactorOrder=covariateFactorOrder,
					                      inverseVarianceWeights=inverseVarianceWeights,
					                      signalProbWeights=signalProbWeights,
					                      analyzeReagents=analyzeReagents,	
                                calculateRelativeDropout=calculateRelativeDropout,
					                      pvalueType=pvalueType,
					                      modelEngine=modelEngine,
					                      variableMap = variableMap,
					                      formulas=formulas,
					                      endPoint=endPoint,
					                      printFreq=NULL);
		}

		# Clean up once parallel processing complete
    stopCluster(cl)

	} else {
		final = simemChunk(screens=screens,
		                   covariate=covariate,
		                   geneIds=geneIds,
		                   analysisType=analysisType,
		                   reagentWeights=reagentWeights,
		                   cellLineWeights=cellLineWeights,
		                   annotationsPerCellLine=annotationsPerCellLine,
		                   annotationsPerId=annotationsPerId,
		                   covariateFactorOrder=covariateFactorOrder,
		                   inverseVarianceWeights=inverseVarianceWeights,
		                   signalProbWeights=signalProbWeights,
		                   analyzeReagents = analyzeReagents,	
                       calculateRelativeDropout=calculateRelativeDropout,
		                   pvalueType=pvalueType,
		                   modelEngine=modelEngine,
		                   variableMap = variableMap,
		                   formulas=formulas,
		                   endPoint=endPoint,
		                   printFreq=printFreq);
	}

	# Correct pvalues for multiple testing
	final[["gene_detailed"]] = applyFDR(final[["gene_detailed"]])
	final[["gene"]] = applyFDR(final[["gene"]])
	
	final[["gene_detailed"]] = formatOutput(final[["gene_detailed"]])
	final[["gene"]] = formatOutput(final[["gene"]])
	
	if(analyzeReagents) {
		final[["reagent_detailed"]] = applyFDR(final[["reagent_detailed"]])
		final[["reagent"]] = applyFDR(final[["reagent"]])
		
		final[["reagent_detailed"]] = formatOutput(final[["reagent_detailed"]])
		final[["reagent"]] = formatOutput(final[["reagent"]])
	}

	return(final)
}



combineResults = function(list1, list2) {

	combined = list()
	
#	sink("~/Desktop/work/data_processing_shrna/2013-05-09_simem_parallel/debug-parallel.txt", append=T)
#	print(str(list1))
#	print(str(list2))
#	sink()

	combined[["gene"]] = rbind(list1[["gene"]], list2[["gene"]])
	combined[["gene_detailed"]] = rbind(list1[["gene_detailed"]], list2[["gene_detailed"]])

	if(!is.null(list1[["reagent"]])) {
	  combined[["reagent"]] = rbind(list1[["reagent"]], list2[["reagent"]])	  
	  combined[["reagent_detailed"]] = rbind(list1[["reagent_detailed"]], list2[["reagent_detailed"]])
	}

	return(combined)
}


####
# Fit gene and reagent level models to a subset of the gene ids, typically performing this analysis
# on a separate processor core for parallelization...
#### 
simemChunk = function(screens, 
                      geneIds, 
                      analysisType="panel", 
                      covariate=NULL, 
          						reagentWeights=NULL, 
          						cellLineWeights=NULL, 
          						annotationsPerCellLine=NULL, 
          						annotationsPerId=NULL, 
          						covariateFactorOrder=NULL,
          						inverseVarianceWeights=TRUE, 
          						signalProbWeights=FALSE, 
          						analyzeReagents=TRUE,
                      calculateRelativeDropout=FALSE,
          						# default is "gaussian", with "t" producing p-values as reported in manuscript
          						pvalueType="gaussian",
          						# Fit siMEM models using the lme4 or blme packages as specified
          						modelEngine="blme",
          						variableMap=getDefaultVariableMap(), 
          						formulas=NULL, 
          						endPoint=FALSE, 
                      printFreq=20
				   ) {

	fitsSummary = list()
	fitsWarnings = vector()
	fitsErrors = vector()

	fitsSummaryRg = list()
	fitsWarningsRg = vector()
	fitsErrorsRg = vector()

	# signalProbWeights are only relevant for 3+ timepoint experiments.
	if(endPoint == TRUE) {
		signalProbWeights=FALSE
	}

	counter = 1

	for(gene_id in geneIds) {

		geneKey = as.character(gene_id)

		if(!is.null(printFreq)) {
		  if(counter %% printFreq == 0) {
		    print(paste("Fitting gene models: ", counter))		    
		  }
		}

		# Extract all measurements associated with the specified gene, and match each measurement
		# to its annotations and genomic covariate. The data frame returned is in "long" format
		# with one measurement per row.
		datgene = extractData(screens, 
		                      idType="gene", 
		                      id=gene_id, 
		                      covariate=covariate, 
            							reagentWeights=reagentWeights, 
                          cellLineWeights=cellLineWeights,
            							annotationsPerCellLine=annotationsPerCellLine,
            							annotationsPerId=annotationsPerId, 
                          covariateFactorOrder=covariateFactorOrder,
            							inverseVarianceWeights=inverseVarianceWeights, 
                          signalProbWeights=signalProbWeights,
                          variableMap=variableMap)


		
		# It may be that, after filtering out reagents frequently starting below noise threshold
		# no reagents for the gene are left. A classic example from the Marcotte et al. 2012 screens is ALK.
		# In that case, just skip storing the empty dataframe and fitting the model
		if(!empty(datgene)) {

			if(length(unique(datgene$reagentId)) >= 3) {

				if(analysisType == "isogenic") {
					fitgene = fitIsogenicCovariateGene(datgene, formulas=formulas, pvalueType=pvalueType, modelEngine=modelEngine, outputLevel="summary")
				} else {
				  #
					fitgene = fitPanelCovariateGene(datgene, formulas=formulas, pvalueType=pvalueType, modelEngine=modelEngine, outputLevel="summary", endPoint=endPoint)
				}

				fitsSummary[[geneKey]] = fitgene[["summary"]]
				fitsWarnings[geneKey] = fitgene[["warning"]]
				fitsErrors[geneKey] = fitgene[["error"]]
			}

			# Only fit reagent-level models if specified (TRUE by default)
			if(analyzeReagents==TRUE) {
				## Fit reagent models at same time, to avoid storing huge data list
				rgs = unique(datgene$reagentId)

				for(reagent_id in rgs) {

					datReagent = datgene[datgene$reagentId == reagent_id,]
					
					# Even if we subset the data frame, some variables will be factors where some factor levels will 
					# have no associated measurements. Having factor levels with no associated measurements messes up
					# the model fitting, since the model will try to fit data for -each- factor level. The reset
					# resolves this issue.
					datReagent = resetFactorLevels(datReagent)

					reagentKey = as.character(reagent_id)

					if(nrow(datReagent) > 0) {

						if(analysisType == "isogenic") {
							fitReagent = fitIsogenicCovariateRg(datReagent, formulas=formulas, pvalueType=pvalueType, modelEngine=modelEngine, outputLevel="summary")
						} else {
							fitReagent = fitPanelCovariateRg(datReagent, formulas=formulas, pvalueType=pvalueType, modelEngine=modelEngine, outputLevel="summary", endPoint=endPoint)
						}

						fitsSummaryRg[[reagentKey]] = fitReagent[["summary"]]
						fitsWarningsRg[reagentKey] = fitReagent[["warning"]]
						fitsErrorsRg[reagentKey] = fitReagent[["error"]]
					}
				}
			}
		}

		counter = counter+1
	}

	fdat = fData(screens)
	annot = unique(fdat[,c(variableMap$geneId, variableMap$symbol)])
	annotRg = unique(fdat[,c(variableMap$reagentId, variableMap$geneId, variableMap$symbol)])

	allfitsSummary = do.call(rbind, fitsSummary)

	# For end-point models, we don't have enough info to calculate the relative dropout rate.
  # For continuous covariates, the RDR makes no sense, so need to exclude
	# 

	minRowsForRDR = 100

	if(calculateRelativeDropout == TRUE & endPoint == FALSE & nrow(allfitsSummary) >= minRowsForRDR) {
		allfitsSummary = getRelativeDropoutRates(allfitsSummary, geneOrReagent="gene")
	}
	
	allfitsMinimal = getMinimalSummary(allfitsSummary, 
	                                   calculateRelativeDropout=calculateRelativeDropout, 
	                                   endPoint=endPoint)

	allfitsSummary = cbind.data.frame(rownames(allfitsSummary), allfitsSummary, stringsAsFactors=F)
	colnames(allfitsSummary)[1] = variableMap$geneId
	allfitsSummary[,variableMap$geneId] = as.integer(allfitsSummary[,variableMap$geneId])
	allfitsSummary = cbind.data.frame(allfitsSummary, warnings=fitsWarnings, errors=fitsErrors, stringsAsFactors=F)
	allfitsSummary = merge(annot, allfitsSummary)

	allfitsMinimal = cbind.data.frame(rownames(allfitsMinimal), allfitsMinimal, stringsAsFactors=F)
	colnames(allfitsMinimal)[1] = variableMap$geneId
	allfitsMinimal[,variableMap$geneId] = as.integer(allfitsMinimal[,variableMap$geneId])
	allfitsMinimal = cbind.data.frame(allfitsMinimal, warnings=fitsWarnings, errors=fitsErrors, stringsAsFactors=F)
	allfitsMinimal = merge(annot, allfitsMinimal)
	
	# Only fit reagent-level models if specified (TRUE by default)
	if(analyzeReagents==TRUE) {

		modeledRgIds = names(fitsSummaryRg)

		allfitsSummaryRg = do.call(rbind, fitsSummaryRg)

		# For end-point models, we don't have enough info to calculate the relative dropout rate.
		if(calculateRelativeDropout == TRUE & endPoint == FALSE & nrow(allfitsSummaryRg) >= minRowsForRDR) {
			allfitsSummaryRg = getRelativeDropoutRates(allfitsSummaryRg, geneOrReagent="reagent")
		}

		allfitsMinimalRg = getMinimalSummary(allfitsSummaryRg, 
		                                     calculateRelativeDropout=calculateRelativeDropout, 
		                                     endPoint=endPoint)

		allfitsSummaryRg = cbind.data.frame(rownames(allfitsSummaryRg), allfitsSummaryRg, stringsAsFactors=F)
		colnames(allfitsSummaryRg)[1] = variableMap$reagentId
		allfitsSummaryRg = cbind.data.frame(allfitsSummaryRg, warnings=fitsWarningsRg, errors=fitsErrorsRg, stringsAsFactors=F)
		allfitsSummaryRg = merge(annotRg, allfitsSummaryRg)

		allfitsMinimalRg = cbind.data.frame(rownames(allfitsMinimalRg), allfitsMinimalRg, stringsAsFactors=F)
		colnames(allfitsMinimalRg)[1] = variableMap$reagentId
		allfitsMinimalRg = cbind.data.frame(allfitsMinimalRg, warnings=fitsWarningsRg, errors=fitsErrorsRg, stringsAsFactors=F)
		allfitsMinimalRg = merge(annotRg, allfitsMinimalRg)
	}

	final = list()
	final[["gene_detailed"]] = allfitsSummary
	final[["gene"]] = allfitsMinimal

	# Only fit reagent-level models if specified (TRUE by default)
	if(analyzeReagents==TRUE) {
		final[["reagent_detailed"]] = allfitsSummaryRg
		final[["reagent"]] = allfitsMinimalRg
	}

	return(final)
}

