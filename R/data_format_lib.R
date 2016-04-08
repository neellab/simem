##################################################################
# Requirements for generic data extraction function:
# - extract arbitrary subset of screens, identified by:
#   cellLineId, condition, any combination of filtering criteria organized as a list (variable => vars)
#
# - optionally: add data column(s) to meta-data extracted from eSet pData()
#   these will be merged with the metadata using unique "sampleId"
#   filters can apply to these as well
#   COVARIATE column may or may not be added using this additional column
#
#   QUESTION: how to reconcile entering an "COVARIATE" matrix such as per-gene CN vs.
#   adding multi-covariate matrix?
#
#   IMPORTANT: check that
#
# - optionally: specify list of reagentIds with weights [0;1].
# Reagents with weight 0 are to be excluded.
# Total number of non-zero-weight reagents stored.
# Data weights renormalized to sum of non-zero weights, after T_i-1, variance weighting
# 
# - optionally: specify a data frame of per-gene or per-reagent values, that will be slotted
# as the COVARIATE column
#
#
# INPUT ARGUMENTS:
# eset -	eset containing array or sequencing-based quantification of screens

# idType -	"gene" or "reagent"

# id -	a variable with an id of the appropriate type

# covariate -	variable column, contained in either the eset phenoData matrix or the
#			"annotations" variable, to include as additional covariate in the
# 			model.

# reagentWeights -	a data frame, which MUST have as 1st column the unique reagent
#					identifier, with header the value of variableMap$reagentId
#					(by default, "trcn_id"). Additionally, at least one column
#					of the input data frame must be labeled "weight".
#					At its simplest, labeling each reagent with 0 or 1 indicates
#					which reagents to exclude from the analysis by some previous
#					filtering or QC analysis. Data points with a weight of 0
#					are ignored by the mixed-effect model. By default, all data
#					points are assumed to have equal weight. This can be changed,
#					for example, by using inverse-variance weights, to incorporate
#					additional information.

# annotationsPerCellLine - additional columns that associate values with the different
#						Cell Lines in the experiment.
#						IMPORTANT: this data frame MUST have as 1st column
#						the value of variableMap$sampleId (by default, "sample_id")

# annotationsPerId - 	additional matrix, indexed by samples for columns, and genes/reagents
#						on rows, that provides gene/reagent specific values per sample.
#						for example, these can be copy number status, or expression value
#						IMPORTANT: this data frame MUST start with either: 
#						- 	the value of variableMap$geneId (by default, "gene_id")
#							if the value of idType parameter is "gene"
#						- 	the value of variableMap$reagentId (by default, "trcn_id")
#							if the value of idType parameter is "reagent"

# variableMap - identifies the column names in the eSet phenoData and annotations
#				table that correspond to the data format required for the Mixed-Effect,
#				and other, R modeling packages.
##################################################################
extractData = function(eset, 
                      idType, 
                      id, 
                      covariate=NULL, 
                      annotationsPerCellLine=NULL,
                      annotationsPerId=NULL,
                      covariateFactorOrder=NULL,
                      reagentWeights=NULL,
                      cellLineWeights=NULL,
                      inverseVarianceWeights=TRUE,
                      signalProbWeights=FALSE,
                      variableMap = getDefaultVariableMap() 
    						) {

	# Load eset data
	shrna = exprs(eset)
	pheno = pData(eset)
	fdat = fData(eset)

	if(inverseVarianceWeights) {
		varWeights = assayData(eset)[["varWeights"]]
	}

	if(signalProbWeights) {
		signalWeights = assayData(eset)[["signalProbWeights"]]
	}

	##################################################################
	# BEGIN - Input variable sanity check
	##################################################################
	# TEST that key columns are present in eset phenoData
	phenoColsToKeep = c(variableMap$sampleId,
						variableMap$cellLineId,
						variableMap$replicateId,
						variableMap$timeGroup,
						variableMap$timeNum)

	if(length(intersect(phenoColsToKeep, colnames(pheno))) != length(phenoColsToKeep)) {

		expected = paste(phenoColsToKeep, sep="", collapse=", ")
		observed = paste(intersect(phenoColsToKeep, colnames(pheno)), sep="", collapse=", ")

		msg = paste("ERROR in extractData():\n\n",
					"eset phenoData data frame must have, at a minimum, the following columns: ", expected,
					"However, only the following were observed: ", observed,
					"\n\nRequired column names can be customized using the 'variableMap' parameter",
					"\n\n\n\n", sep="")
		stop(msg)
	}

	if(!(idType %in% c("gene", "reagent"))) {
		msg = paste("ERROR in extractData():\n\n",
					"Incorrect 'idType' value (either 'gene' or 'reagent' allowed): ", idType,
					"\n\n\n\n", sep="")
		stop(msg)
	}

	if(!is.data.frame(variableMap) || nrow(variableMap) != 1) {
		msg = paste("ERROR in extractData():\n\n",
					"Incorrect 'variableMap' input (data frame with single row required)",
					"\n\n\n\n", sep="")
		stop(msg)
	}

	# TEST sampleId labeling
	# check that annotations table, if specified, has sampleId as first column
	if(is.data.frame(reagentWeights) &&
	   (colnames(reagentWeights)[1] != variableMap$reagentId ||
	   length(grep("weight", colnames(reagentWeights))) == 0)) {
		msg = paste("ERROR in extractData():\n\n",
					"Reagent weights table must have as first column: ",
					variableMap$reagentId,
					"\n\nThe data frame must also contain a 'weight' column.",
					"\n\nReagent identifier can be customized using the 'variableMap' parameter", "\n\n\n\n", sep="")
		stop(msg)
	}

	# TEST sampleId labeling
	# check that annotations table, if specified, has sampleId as first column
	if(is.data.frame(annotationsPerCellLine) && colnames(annotationsPerCellLine)[1] != variableMap$cellLineId) {
		msg = paste("ERROR in extractData():\n\n",
					"Additional 'annotationsPerCellLine' table must have as first column ", variableMap$cellLineId,
					"\n\nRequired column names can be customized using the 'variableMap' parameter",
					"\n\n\n\n", sep="")
		stop(msg)
	}

	# If the additional annotations per cell line table has non-id columns that overlap
	# with the pheno column, warn of this and drop the additional columns from the pheno matrix
	# TEST sampleId labeling
	# check that annotations table, if specified, has sampleId as first column
	if(is.data.frame(annotationsPerCellLine)) {
		overlap = setdiff(
						  intersect(colnames(pheno), colnames(annotationsPerCellLine)),
						  variableMap$cellLineId)
		# Warn of duplicate column names, and drop them from the pheno matrix.
		if(length(overlap) > 0) {

			colsToRemove = paste(overlap, sep="", collapse="|")
			superceded = paste(overlap, sep="", collapse=", ")

			pheno = pheno[,-grep(colsToRemove, colnames(pheno))]
		
			msg = paste0("WARNING in extractData():\n\n",
						"Duplicate columns other than cell line identifier shared between eset phenoData and 'annotationsPerCellLine'.\n",
						"The following columns of phenoData were superceded by 'annotationsPerCellLine' input to avoid data merging issues:\n",
						superceded,
						"\n\n\n\n")
			warning(msg)
		}
	}


	# If we're specifying a per gene/reagent set of annotations (eg: Copy Number status)
	# ensure that, if we're looking at
	# gene-level data, the first column of annotationsById is "gene_id"
	# reagent-level data, the first column of annotationsById is "reagent_id"
	####if(is.data.frame(annotationsPerId) &&
	####	(	((idType == "gene") && (colnames(annotationsPerId)[1] != variableMap$geneId)) ||
	####		((idType == "reagent") && (colnames(annotationsPerId)[1] != variableMap$reagentId))	)) {
	####	msg = paste("ERROR in extractData():\n\n",
	####				"Additional 'annotationsPerId' table must have as first column either ",
	####				variableMap$geneId, " or ", variableMap$reagentId,
	####				"\n\nRequired column names can be customized using the 'variableMap' parameter", "\n\n\n\n", sep="")
	####	stop(msg)
	####}
	if(is.data.frame(annotationsPerId) && (colnames(annotationsPerId)[1] != variableMap$geneId)) {
		msg = paste("ERROR in extractData():\n\n",
					"Additional 'annotationsPerId' table must have as first column ", variableMap$geneId, 
					"\n\nRequired column names can be customized using the 'variableMap' parameter", "\n\n\n\n", sep="")
		stop(msg)
	}

	# Ensure that the id we're extracting data for has a unique entry in annotationsPerId
	if(is.data.frame(annotationsPerId)) {

		idVar = variableMap$geneId

		if(idType == "reagent") {
			gene = fdat[fdat[,variableMap$reagentId] == id, variableMap$geneId]
		} else {
			gene = id
		}

		checkIdMatch = which(annotationsPerId[,variableMap$geneId] == gene)
		matchCount = length(checkIdMatch)

		if(matchCount != 1) {
			msg = paste("ERROR in extractData():\n\n",
						"While checking 'annotationsPerId' data frame for unique matching gene or reagent id.", 
						"\n\nid specified: ", id,
						"\n\nNumber of matching rows found for id in 'annotationsPerId': ", matchCount,
						"\n\nPlease ensure that the first column of 'annotationsPerId' contains an entry for the specified id.",
						"\n\n\n\n", sep="")
			stop(msg)
		}
	}

	# TEST sampleId labeling for the phenoData of the eset
	if(length(intersect(colnames(pheno), variableMap$sampleId)) != 1) {
		observed = paste(colnames(pheno), sep="", collapse=", ")
		msg = paste("ERROR in extractData():\n\n",
					"At least one of the eset phenoData columns must be: ", variableMap$sampleId,
					"\n\nWe observe the following columns in the phenoData: ", observed,		
					"\n\nRequired column names can be customized using the 'variableMap' parameter", "\n\n\n\n", sep="")
		stop(msg)
	}

	# TEST gene/reagent ids
	# test that gene and reagent ids are properly set in the fData of the eset
	idVars = c(variableMap$reagentId, variableMap$geneId, variableMap$symbol)
	idMatch = intersect(colnames(fdat), idVars)

	if(length(idMatch) != length(idVars)) {
		expected = paste(idVars, sep="", collapse=", ")
		observed = paste(colnames(fdat), sep="", collapse=", ")
		msg = paste("ERROR in extractData():\n\n",
					"Cannot find reagent id, gene id or gene symbol in eset featureData.",
					"\n\nWe expect to find the following columns: ", expected,
					"\n\nWe observe the following columns in the featureData: ", observed,
			"\n\nId variable columns can be customized using the 'variableMap' parameter", "\n\n\n\n", sep="")
		stop(msg)
	}

	#############################################
	# MERGE ADDITIONAL ANNOTATIONS into the pheno matrix if theye're provided
	#############################################
	if(is.data.frame(annotationsPerCellLine)) {
		# add additional columns to the pheno on the fly
		pheno = merge(pheno, annotationsPerCellLine)
		# Merge may not preserve row names, and these are used later on
		rownames(pheno) = pheno[,variableMap$sampleId]
	}

	# If a per-id annotation table is supplied, pull out the annotation matching
	# the specified id, and put the "per-id" values into the "COVARIATE" column
	if(is.data.frame(annotationsPerId)) {
		idVar = variableMap$geneId

		if(idType == "reagent") {
			gene = fdat[fdat[,variableMap$reagentId] == id, variableMap$geneId]
		} else {
			gene = id
		}

		values = annotationsPerId[annotationsPerId[,idVar] == gene,-1,drop=F]
		values = t(values)
		values = cbind.data.frame(rownames(values), values, stringsAsFactors=F)
		colnames(values)[1:2] = c(variableMap$cellLineId, "COVARIATE")

		# add additional columns to the pheno on the fly
		pheno = merge(pheno, values)
		rownames(pheno) = pheno[,variableMap$sampleId]
	}

	# Ensure covariate variable exists in either, provided it hasn't been extracted
	# from annotationsPerId just above
	if(is.character(covariate) && !is.data.frame(annotationsPerId) && length(intersect(colnames(pheno), covariate)) != length(covariate)) {
		msg = paste("ERROR in extractData():\n\n",
					"Cannot find this 'covariate' column in the eset phenoData (or annotationsPerCellLine, if specified): ", covariate,
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
# 			msg = paste("ERROR in extractData():\n\n",
# 						"Some specified filter variables not found in phenoData or additional annotations.",
# 						"\n\nFilter variables matching existing phenoData variables: ", expected,
# 						"\n\nFilter variables missing from phenoData: ", observed, "\n\n\n\n", sep="")
# 			stop(msg)
# 		}
# 	}

	if(!is.null(covariate) && is.data.frame(annotationsPerId)) {
		msg = paste("WARNING in extractData(): Both 'covariate' and 'annotationsPerId' parameters were specified.",
					"\n\nAs a result, 'covariate' parameter will be ignored and the id-specific values pulled from annotationsPerId will be used instead.",
					"\n\n\n\n", sep="")
		warning(msg) 
	}
	##################################################################
	# END - Input variable sanity check
	##################################################################


	##################################################################
	# BEGIN - Subset samples to analyze using filters, keep only relevant columns
	##################################################################

#	# include all samples by default
#	sampleFilter = rep(TRUE, nrow(pheno))

#	# APPLY FILTER - get the samples satisfying all filtering criteria
#	if(is.list(filters)) {
#		for(filterVar in names(filters)) {
#			# AND each new filter condition with previous ones
#			sampleFilter = sampleFilter & (!is.na(pheno[,filterVar]) & (pheno[,filterVar] %in% filters[[filterVar]]))
#		}
#	}

# 	# Restrict to a subset of samples, and to relevant columns
# 	# Keep basic information, plus any variables filtered on, as well as the
# 	# COVARIATE variable (the covariate, such as ERBB2 amp)
# 	if(is.list(filters)) {
# 		phenoColsToKeep = union(phenoColsToKeep, names(filters))
# 	}

	# We took the row from the annotationsPerId table corresponding to our entry
	if(is.data.frame(annotationsPerCellLine)) {
		phenoColsToKeep = union(phenoColsToKeep, colnames(annotationsPerCellLine))
	}

	# We took the row from the annotationsPerId table corresponding to our entry
	if(is.data.frame(annotationsPerId)) {
		phenoColsToKeep = union(phenoColsToKeep, "COVARIATE")
	}

	# We preferentially extract the covariate column from the perId annotations table
	# if specified, otherwise use the "covariate" variable if specified
#	if(is.character(covariate) && !is.data.frame(annotationsPerId)) {
#		phenoColsToKeep = union(phenoColsToKeep, covariate)
#	}

	# Keep only the columns required for the analysis
	pheno = pheno[,phenoColsToKeep]

	# keep only a subset of samples, and columns relevant to the analysis
#	pheno = pheno[sampleFilter, phenoColsToKeep]
	##################################################################
	# END - Subset samples to analyze using filters, keep only relevant columns
	##################################################################


	##################################################################
	# BEGIN - Subset shRNA intensity/count data to gene/reagent of interest
	##################################################################
	# extract rows corresponding to gene or reagent of interest
	if(idType == "reagent") {
		filtered = which(fdat[,variableMap$reagentId] == id)
	} else if(idType == "gene") {
		filtered = which(fdat[,variableMap$geneId] == id)
	} else {
		msg = paste("ERROR in extractData():\n\n",
					idType, " id value specified does not match any in the eset featureData: ", id,
					"\n\n\n\n", sep="")
		stop(msg)
	}

	gene_id = unique(fdat[filtered,variableMap$geneId])
	symbol = unique(fdat[filtered,variableMap$symbol])

	# extract reagent/gene ids for subset of screens matching our filters
	# before merging data with annotations, ensure that properly labeled sample id
	# column exists
	tempWide = shrna[filtered,rownames(pheno),drop=F]
	temp = t(tempWide)
	temp = cbind.data.frame(rownames(temp), temp, stringsAsFactors=F)
	colnames(temp)[1] = variableMap$sampleId

	##################################################################
	# END - Subset shRNA intensity/count data to extract gene/reagent of interest
	##################################################################


	##################################################################
	# BEGIN - reshape data frame to associate annotations with each intensity measurement
	##################################################################
	# merge on sampleId, assuming column/sample names in intensity matrix don't match
	# column names in phenoData (sample names can reasonably be assumed distinct from
	# "doublings", "cell_line" and other such phenoData column names)
	dat = merge(pheno, temp)
	melted = melt(dat, id.vars=phenoColsToKeep, variable_name=variableMap$reagentId)

	# Make the variable name more generic, in case we want to use TRPS
	melted$geneId = rep(gene_id, nrow(melted))
	melted$symbol = rep(symbol, nrow(melted))

	# Standardize column names so that model formulae can be specified using these
	# consistent column names
	cols = standardizeColumnNames(colnames(melted), variableMap=variableMap)

	# If we require an covariate column (covariate input param is set)
	# and the COVARIATE column isn't already specified (eg: from annotationsPerId)
	if(is.character(covariate) && length(covariate) == 1 && length(grep("COVARIATE", cols)) == 0) {		
		cols = gsub(paste0("^",covariate,"$"), "COVARIATE", cols)

	# Alternatively, if 2+ covariates were specified, and no annotationsPerId value was specified, 
	} else if(is.character(covariate) && length(covariate) > 1 && length(grep("COVARIATE", cols)) == 0) {
		cols = gsub(paste0("^",covariate[1],"$"), "COVARIATE", cols)
		cols = gsub(paste0("^",covariate[2],"$"), "COVARIATE2", cols)

	# Alternatively, if annotationsPerId value was specified, the per-gene values will already
	# have been placed into the COVARIATE column, if an additional covariate has been specified
	# rename that column as COVARIATE2. If length(covariate) argument > 1, ignore further values
	} else if(is.character(covariate) && length(grep("COVARIATE", cols)) == 1) {
		cols = gsub(paste0("^",covariate[1],"$"), "COVARIATE2", cols)

	} else if(is.character(covariate) && length(grep("COVARIATE", cols)) == 1) {
		cols = gsub(paste0("^",covariate[1],"$"), "COVARIATE2", cols)

	}
#	else {
#		# WARNING? Or just remove this clause?
#		warning("\n\nUnhandled condition in if-else block that sets COVARIATE columns...")
#	}

	# the melt() function puts the intensity into a column named "value"
	cols = gsub("^value$", "expr", cols)

	colnames(melted) = cols

	melted$cellLineIdUnique = paste(melted$reagentId, melted$cellLineId, sep="-")
	melted$replicateIdUnique = paste(melted$reagentId, melted$cellLineId, melted$replicateId, sep="-")

	##################################################################
	# END - reshape data frame to associate annotations with each intensity measurement
	##################################################################

	meltedVarWts = NA
	meltedNoiseWts = NA

	if(inverseVarianceWeights) {
		varWeights = assayData(eset)[["varWeights"]]
		# extract reagent/gene ids for subset of screens matching our filters
		# before merging data weights with annotations, ensure that properly labeled 
		# sample id column exists
		varWts = varWeights[filtered,rownames(pheno),drop=F]
		varWts = t(varWts)
		varWts = cbind.data.frame(rownames(varWts), varWts, stringsAsFactors=F)
		colnames(varWts)[1] = variableMap$sampleId
		datVarWts = merge(pheno, varWts)
		meltedVarWts = melt(datVarWts, id.vars=phenoColsToKeep, variable_name=variableMap$reagentId)
		meltedVarWts = standardizeColumnNames(meltedVarWts, variableMap=variableMap)
		colnames(meltedVarWts) = gsub("^value$", "variance_weights", colnames(meltedVarWts))
	}

	if(signalProbWeights) {
		signalWeights = assayData(eset)[["signalProbWeights"]]
		# extract reagent/gene ids for subset of screens matching our filters
		# before merging data weights with annotations, ensure that properly labeled 
		# sample id column exists
		signalWts = signalWeights[filtered,rownames(pheno),drop=F]
		signalWts = t(signalWts)
		signalWts = cbind.data.frame(rownames(signalWts), signalWts, stringsAsFactors=F)
		colnames(signalWts)[1] = variableMap$sampleId
		datNoiseWts = merge(pheno, signalWts)
		meltedNoiseWts = melt(datNoiseWts, id.vars=phenoColsToKeep, variable_name=variableMap$reagentId)
		meltedNoiseWts = standardizeColumnNames(meltedNoiseWts, variableMap=variableMap)
		colnames(meltedNoiseWts) = gsub("^value$", "signal_weights", colnames(meltedNoiseWts))
	}

	# CALCULATE WEIGHTS, using a combination of reagent-specific, cell line-specific
	# inverse-variance weights and microarray noise (signal probability) weights
	melted = assignWeightsToObservations(dat=melted,
										reagentWeights=reagentWeights,
										cellLineWeights=cellLineWeights,
										inverseVarianceWeights=inverseVarianceWeights,
										signalProbWeights=signalProbWeights,
										varWeightsDF=meltedVarWts,
										signalWeightsDF=meltedNoiseWts,
										variableMap=variableMap)

	colSort = sort(colnames(melted))
	melted = melted[,colSort]

	# If the covariate is categorical, cast it as a Factor to facilitate modeling
	if(!is.null(melted$COVARIATE) && is.character(melted$COVARIATE)) {

		observedCovariateValues = unique(melted$COVARIATE)

		#if necessary, reorder the factor levels in COVARIATE, which are by default alphabetical
		# However, both vectors must be characters, and all observed COVARIATE values must be
		# accounted for in the covariateFactorOrder
		if(is.character(covariateFactorOrder) && is.character(melted$COVARIATE)
		   && length(intersect(covariateFactorOrder, melted$COVARIATE)) == length(observedCovariateValues)) {

			subsetOrder = covariateFactorOrder[covariateFactorOrder%in% observedCovariateValues]
			# Reset factor levels, for example, if the first level should be "NON-AMP"
			# and we want to compare the covariate of "AMP"
			melted$COVARIATE = factor(melted$COVARIATE, levels=subsetOrder)
		} else {
			melted$COVARIATE = factor(melted$COVARIATE)
		}
	}

	# If we have 2 covariates, for now the second covariate will be handled alphabetically
	if(!is.null(melted$COVARIATE2) && is.character(melted$COVARIATE2)) {
		melted$COVARIATE2 = factor(melted$COVARIATE2)
	}

	return(melted)
}

# Map column names in data to standard ones to be used in model fitting process
standardizeColumnNames = function(dat, variableMap=getDefaultVariableMap()) {

	if(is.vector(dat)) {
		cols = dat
	} else if(is.matrix(dat) | is.data.frame(dat)) {
		cols = colnames(dat)
	}

	cols = gsub(paste0("^",variableMap$sampleId,"$"), "sampleId", cols)
	cols = gsub(paste0("^",variableMap$cellLineId,"$"), "cellLineId", cols)
	cols = gsub(paste0("^",variableMap$timeGroup,"$"), "timeGroup", cols)
	cols = gsub(paste0("^",variableMap$timeNum,"$"), "timeNum", cols)
	cols = gsub(paste0("^",variableMap$replicateId,"$"), "replicateId", cols)
	cols = gsub(paste0("^",variableMap$reagentId,"$"), "reagentId", cols)

	if(is.vector(dat)) {
		dat = cols
	} else if(is.matrix(dat) | is.data.frame(dat)) {
		colnames(dat) = cols
	}

	return(dat)
}



# The reagent/cell line weights are extracted from input data frames if required
# inverse variance and signal weights have already been extracted and added as
# columns to 'dat', if the options were specified
assignWeightsToObservations = function(dat,
                  									   reagentWeights,
                  									   cellLineWeights,
                  									   inverseVarianceWeights,
                  									   varWeightsDF,
                  									   signalProbWeights,
                  									   signalWeightsDF,
                  									   variableMap) {
	##################################################################
	# BEGIN - Calculate final weight of each data point in the dataset
	# Each type of weight should be renormalized so the total weight
	# of all observations is identical to the number of observations.
	# The alternative is to assig absolute weights, which could lead
	# to implausible results.
	##################################################################
	weightRG = is.data.frame(reagentWeights)
	weightCL = is.data.frame(cellLineWeights)
	weightVar = is.data.frame(varWeightsDF)
	weightSignal = is.data.frame(signalWeightsDF)

	if((weightRG | weightCL | weightVar | weightSignal) & nrow(dat) > 0) {
		# Default all data points to equal value
		dat$W=1

		#If specified, attach per-reagent weights to the data
		if(weightRG) {
			# By this point, the reagent column is identified by the string "reagentId"
			# regardless of the original column label
			rgs = unique(dat$reagentId)
			colnames(reagentWeights)[1] = "reagentId"
			rgWeights = reagentWeights[reagentWeights$reagentId %in% rgs, c("reagentId", "weight"), drop=F]
			colnames(rgWeights)[2] = "reagent_weights"

			#append reagent-specific weights
			dat = merge(dat, rgWeights, by="reagentId", all.x=T)
	
			#remove measurements associated with reagents of weight 0.
			toKeepRG = !(is.na(dat$reagent_weights)) & dat$reagent_weights > 0
			dat = dat[toKeepRG,]
		}

		#If specified, attach per-reagent weights to the data
		if(weightCL) {
			# By this point, the reagent column is identified by the string "reagentId"
			# regardless of the original column label
			cellLines = unique(dat$cellLineId)
			colnames(cellLineWeights)[1] = "cellLineId"
			clWeights = cellLineWeights[cellLineWeights$cellLineId %in% cellLines, c("cellLineId", "weight"), drop=F]
			colnames(clWeights)[2] = "cell_line_weights"

			#append reagent-specific weights
			dat = merge(dat, clWeights, by="cellLineId", all.x=T)

			#remove measurements associated with cell lines of weight 0.
			toKeepCL = !(is.na(dat$cell_line_weights)) & dat$cell_line_weights > 0
			dat = dat[toKeepCL,]
		}

		if(weightVar) {
			dat = merge(dat, varWeightsDF)
			dat = dat[!is.na(dat$variance_weights) & dat$variance_weights > 0,]
		}

		if(weightSignal) {
			dat = merge(dat, signalWeightsDF)
			dat = dat[!is.na(dat$signal_weights) & dat$signal_weights > 0,]
		}

		# Total number of non-zero observations
		totalObservations = nrow(dat)
    
		if(weightRG) {
			totalRG = sum(dat$reagent_weights)
			rescalingRG = totalObservations/totalRG
			dat$reagent_weights = dat$reagent_weights * rescalingRG
			dat$W = dat$W * dat$reagent_weights
		}

		if(weightCL) {
			totalCL = sum(dat$cell_line_weights)
			rescalingCL = totalObservations/totalCL
			dat$cell_line_weights = dat$cell_line_weights * rescalingCL
			dat$W = dat$W * dat$cell_line_weights
		}

		if(weightVar) {
			totalVar = sum(dat$variance_weights)
			rescalingVar = totalObservations/totalVar
			dat$variance_weights = dat$variance_weights * rescalingVar
			dat$W = dat$W * dat$variance_weights
		}

		if(weightSignal) {
			totalSignal = sum(dat$signal_weights)
			rescalingSignal = totalObservations/totalSignal
			dat$signal_weights = dat$signal_weights * rescalingSignal
			dat$W = dat$W * dat$signal_weights
		}

		# Rescale total weight sum to be equal to original number of data points
		# Otherwise, one can arbitrarily increase significance of regression 
		# parameters simply by multiplying all weightings by 2, or 100, for example.
		# if total = sum(weights)
		# then multiplying all weights by desired/total ensures
		# that the data point weights sum to the desired value,
		# in this case the original number of data points.
		total = sum(dat$W)
		rescaling = totalObservations/total
		dat$W = dat$W * rescaling

		# If certain cell lines/reagents are entirely excluded from the data
		# reset the factor information so that factor levels with no data are excluded
		dat = resetFactorLevels(dat)
	}

	return(dat)
}





###
# If we subset the model data's data frame, some variables will be factors where some factor levels 
# will have no associated measurements. Having factor levels with no associated measurements messes up
# the model fitting, since the model will try to fit data for -each- factor level. The reset
# resolves this issue.
###
resetFactorLevels = function(dat) {
	# Make sure no non-represented factor values exist
	dat$reagentId = factor(as.character(dat$reagentId))
	dat$cellLineId = factor(as.character(dat$cellLineId))
	dat$cellLineIdUnique = factor(as.character(dat$cellLineIdUnique))
	dat$replicateIdUnique = factor(as.character(dat$replicateIdUnique))
	return(dat)
}


getDefaultVariableMap = function() {

	# screenId - a unique identifier for each screen. Cell lines screened under different
	# 			perturbations must have distinct screenIds (eg: cell line + perturbation)
	# cellLineId - cell line name. The same cell line can be perturbed in different
	# 				ways (Eg: drug screen), leading to multiple different screenId
	#				associated with the same cell line

	variableMap = data.frame(geneId = "gene_id",
                          symbol="symbol",
                          reagentId="trcn_id", 
                          cellLineId="cell_line",
                          sampleId="sample_id",
                          replicateId="replicate",
                          timeGroup="timepointStd",
                          timeNum="timepointIndex",
                          condition="condition",
                          screenId="screen_id", 
                          replicateGroupId="replicate_group",
                          stringsAsFactors=F)
	return(variableMap)
}

getMeanSd = function(dat, meanSdTable, variableMap=getDefaultVariableMap()) {
	# check that 2 columns are passed - means, and timepoints associated with means
	colsDat = c("expr", variableMap$timeGroup)
	cols = c("mu", "sigma", variableMap$timeGroup)

#	colsDat = c("expr", "timepoint")
#	cols = c("mu", "sigma", "timepoint")

	if(length(intersect(colnames(dat), colsDat)) != 2) {
		expected = paste(colsDat, sep="", collapse=", ")
		observed = paste(colnames(dat), sep="", collapse=", ")
		msg = paste("ERROR in getMeanSd():\n\n",
					"The following 2 columns are expected in the 'dat' parameter: ", expected,
					"\n\nHowever, we observe the following columns in 'dat': ", observed,
					"\n\nThe 'expr' column must be named as such and is required.",			
					"\n\nThe ", variableMap$timeGroup ," column name can be changed using the 'variableMap' parameter.",
					"\n\n\n\n", sep="")
		stop(msg)
	}

	if(length(intersect(colnames(meanSdTable), cols)) != 3) {
		expected = paste(cols, sep="", collapse=", ")
		observed = paste(colnames(meanSdTable), sep="", collapse=", ")
		msg = paste("ERROR in getMeanSd():\n\n",
					"The following 3 columns are expected in the meanSdTable: ", expected,
					"\n\nHowever, we observe the following columns in the meanSdTable: ", observed,
					"\n\nThe 'mu' and 'sigma' columns must be named as such and are required.",			
					"\n\nThe ", variableMap$timeGroup ," column name can be changed using the 'variableMap' parameter.",
					"\n\n\n\n", sep="")
		stop(msg)
	}

	# Get the timepoint values closest to those available in the specified table
	datTime = sort(unique(dat[,variableMap$timeGroup]))
	meanSdTime = sort(unique(meanSdTable[,variableMap$timeGroup]))

	timeMap = sapply(datTime, function(v, times) {times[which.min(abs(times-v))]}, times=meanSdTime)
	timeMap = cbind.data.frame(datTime=datTime, meanSdTime=meanSdTime)
	dat2 = merge(dat, timeMap, all.x=T, by.x=variableMap$timeGroup, by.y=datTime, sort=F)

	print(table(dat$expr == dat2$expr, useNA="always"))

	mapped = list()

	# perform mean-SD mapping for each timepoint separately
	# eg: if mean-SD function provided for times 0,1,2, perform 3 mapping iterations
	for(onetime in meanSdTime) {
		# get values to map, and mapping table, for the given time value
		onetimeDat = dat2[dat2$datTime == onetime,]
		onetimeTab = meanSdTime[variableMap$timeGroup == onetime,]
		# get table lookup indices for the mean expression values by mapping each
		# to the closest 'mean' value entry in the mean-SD mapping table
		meanIndices = sapply(onetimeDat$expr,
							 function(v, increments) {which.min(abs(increments-x))},
							 increments=unlist(onetimeTab[,variableMap$timeGroup]))
		# for the mapped indices, return the matching Std. Deviation
		# This sd value is presumably calculated as a smoothed version of the
		# individual probe replicate mean-SD values
		onetimeDat$sd = onetimeTab[meanIndices, "sigma"]
		
		# Append the dat table, with Std. Dev. values added, to the list
		mapped = c(mapped, list(onetimeDat))
	}

	final = do.call(rbind, mapped)
	return(final)
}

# Load empirical mean-variance tables and convert them to functions using spline interpolation
loadMeanVar = function(meanVarFile, minVar = 0.1) {
	meanVar = read.delim(meanVarFile, header=T, as.is=T, check.names=F)

	# Limit how low the median variance goes
	minSd = minVar^(1/2)
	meanVar$variance = ifelse(meanVar$variance < minVar, minVar, meanVar$variance)
	meanVar$sigma = ifelse(meanVar$sigma < minSd, minSd, meanVar$sigma)

	timeGroups = unique(meanVar$timeGroup)

	fns = list()
	meanVarFns = list()
	meanSdFns = list()

	for(timept in timeGroups) {
		dat = meanVar[meanVar$timeGroup == timept, c("mu", "sigma", "variance")]
		meanVarFns[[as.character(timept)]] = splinefun(x=dat$mu, y=dat$variance, method="natural")
		meanSdFns[[as.character(timept)]] = splinefun(x=dat$mu, y=dat$sigma, method="natural")
	}

	fns[["sd"]] = meanSdFns
	fns[["var"]] = meanVarFns

	return(fns)
}


formatOutput = function(outputDF) {

  # rename dropout_rate_estimate as baseline_rate
  # rename dropout_rate_diff_X_estimate as difference_X
  # rename relative_dropout_rate_X as rdr_
  # rename_dropout_rate_diff_X_[pvalue|fdr] as [pvalue|fdr]_X
  colnames(outputDF) = gsub("dropout_rate_estimate", "baseline_trend", colnames(outputDF))
  colnames(outputDF) = gsub("dropout_rate_(se|df|t|pvalue|fdr)$", "baseline_\\1", colnames(outputDF))
  colnames(outputDF) = gsub("dropout_rate_diff_(.*)_estimate", "difference_\\1", colnames(outputDF))
  colnames(outputDF) = gsub("dropout_rate_diff_(.*)_(se|df|t|pvalue|fdr)$", "\\2_\\1", colnames(outputDF))

  # In the case of end-point columns, we have columns of the format X_[estimate|pvalue|fdr]
  colnames(outputDF) = gsub("intercept_estimate", "intercept", colnames(outputDF))
  colnames(outputDF) = gsub("([a-zA-Z0-9]*)_estimate", "difference_\\1", colnames(outputDF))
  colnames(outputDF) = gsub("([a-zA-Z0-9]*)_(se|df|t|pvalue|fdr)$", "\\2_\\1", colnames(outputDF))

  # If column name ends with _, remove it. (This occurs if the covariate is continuous, e.g., expression)
  colnames(outputDF) = gsub("_$", "", colnames(outputDF))

  newOrder = c(
                setdiff(colnames(outputDF), c("warnings", "errors")),
                c("warnings", "errors")
               )

  final = outputDF[,newOrder]
  return(final)
}
