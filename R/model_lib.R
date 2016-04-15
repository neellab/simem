## Azin's note: error handling code adapted from R core team/Martin Maechler
## as noted below.
##================================================================##
###  In longer simulations, aka computer experiments,		 ###
###  you may want to						 ###
###  1) catch all errors and warnings (and continue)		 ###
###  2) store the error or warning messages			 ###
###								 ###
###  Here's a solution	(see R-help mailing list, Dec 9, 2010):	 ###
##================================================================##
##' We want to catch *and* save both errors and warnings, and in the case of
##' a warning, also keep the computed result.
##'
##' @title tryCatch both warnings and errors
##' @param expr
##' @return a list with 'value' and 'warning', where
##'   'value' may be an error caught.
##' @author Martin Maechler
#  Copyright (C) 2010 The R Core Team
tryCatchWarningsErrors <- function(expr)
{
	result = NULL
    W = NULL
	E = NULL
	# Warning Handler
    warningHandler = function(w){W <<- w; invokeRestart("muffleWarning")}
	executed = withCallingHandlers(tryCatch(expr, error=function(e) e), warning = warningHandler)

	if("error" %in% class(executed)) {
		E = executed
	} else {
		result = executed
	}

    final = list(
				result = result,
				error = E,
				warning = W)
	return(final)
}


fitModel = function(dat, modelFormula, modelEngine="blme", prior="wishart") {

	model = NULL

	if(modelEngine == "blme") {
  	# if weights are specified, fit a weighted model, otherwise fit an unweighted model.
  	if(is.null(dat$W)) {
  
  		model = tryCatchWarningsErrors(
  					blmer(modelFormula, data=dat, REML=FALSE,
  					      control=lmerControl(optimizer="Nelder_Mead"),
  							fixef.prior = NULL,
  							var.prior = NULL,
  							cov.prior = prior)
  					)
  	} else {

  		# Fit a data-weighted hierarchical model with weakly informative priors.
  		model = tryCatchWarningsErrors(
  					blmer(modelFormula, data=dat, weights=W, REML=FALSE,
  					      control=lmerControl(optimizer="Nelder_Mead"),
  							fixef.prior = NULL,
  							var.prior = NULL,
  							cov.prior = prior)
  					)
  	}
	} else {
	    ### Fit linear mixed model using the lme4 package, without any weakly informative priors
	
  	  # if weights are specified, fit a weighted model, otherwise fit an unweighted model.
  	  if(is.null(dat$W)) {
  	    
  	    # Fit a data-weighted hierarchical model with no (ie: uniform) priors.
  	    model = tryCatchWarningsErrors(
  	      lmer(modelFormula, 
  	           data=dat, 
  	           control=lmerControl(optimizer="Nelder_Mead"),
  	           REML=FALSE)
  	    )
  	  } else {

  	    # Fit a data-weighted hierarchical model with no (ie: uniform) priors.
  	    model = tryCatchWarningsErrors(
  	      lmer(modelFormula, 
  	           data=dat, 
  	           weights=W, 
  	           control=lmerControl(optimizer="Nelder_Mead"),
  	           REML=FALSE)
  	    )
  	  }
	  }

	return(model)
}

###
# Compare 2 nested models using Likelihood Ratio Test, with error handling
### 
fitAnova = function(fit0, fit1) {
	anov = NULL
	anov = tryCatchWarningsErrors(anova(fit0, fit1))
	return(anov)
}

###
# Extract and return measures of model fit to data (AIC, log-Likelihood)
# These are useful to compare alternate, possibly non-nested, model structures
###
summaryFitLogLik = function(fit) {

	if(is.null(fit)) {
		lst = data.frame(AIC=NA, logLik=NA)
	} else {
		aic = extractAIC(fit)[2]
		logL = as.numeric(logLik(fit))
		lst = data.frame(AIC=aic, logLik=logL)
	}

	return(lst)
}


summaryFitEffectCoefs = function(fit) {

	colcount = 8

	if(is.null(fit)) {

		lst = data.frame(matrix(NA, nrow=1, ncol=colcount))

	} else {

		coefs = coef(summary(fit))

		rownames(coefs) = gsub("I\\(.*\\^2\\)", "timeNum2", rownames(coefs))
		rownames(coefs) = gsub("[\\(\\)]", "", rownames(coefs))
		rownames(coefs) = gsub(":", "_", rownames(coefs))
		rownames(coefs) = gsub("COVARIATE2", "", rownames(coefs))
		rownames(coefs) = gsub("COVARIATE", "", rownames(coefs))
		colnames(coefs) = c("estimate", "SE", "t")

		flat = list()

		for(i in 1:nrow(coefs)) {
				currRow = coefs[i,,drop=F]
				colnames(currRow) = paste(rownames(currRow), colnames(currRow), sep="_")
				flat[[i]] = currRow
		}

		lst = do.call(cbind,flat)

		cols = tolower(colnames(lst))

		cols = gsub("timenum_(estimate|t|se)$", "dropout_rate_\\1", cols)
		# All other estimates are differences relative to the baseline covariate level
		cols = gsub("timenum_(.*)_(estimate|t|se)$", "dropout_rate_diff_\\1_\\2", cols)

		cols = gsub("timenum2_(estimate|t|se)$", "dropout_quadratic_\\1", cols)
		# All other estimates are differences relative to the baseline covariate level
		cols = gsub("timenum2_(.*)_(estimate|t|se)$", "dropout_quadratic_diff_\\1_\\2", cols)

		colnames(lst) = cols
	}

	return(lst)
}


summaryAnova = function(anov) {

	colcount = 3

	if(is.null(anov)) {
		lst = data.frame(matrix(NA, nrow=1, ncol=colcount))
	} else {
		chisq = round(anov$`Chisq`[2], digits=2)
		dfs = anov$`Chi Df`[2]
		pvalue_lrt = anov$`Pr(>Chisq)`[2]

		lst = data.frame(chisq, dfs, pvalue_lrt)
	}

	colnames(lst) = c("chisq", "chisq_df", "pvalue_lrt")
	return(lst)
}

###
# Common model fitting function called for any model formula.
# Other functions are responsible for defining the type of model
# to fit, and calling this function to do the grunt work.
###
fitModelAndFormatOutput = function(dat, 
                                  modelFormula,
                                  #modelFormulaNoCovariate,
                                  modelType,
                                  pvalueType="gaussian",
                                  modelEngine="blme",
                                  prior="wishart",
                                  outputLevel=c("summary")) {

	allLevels = c("minimal", "summary", "full")
	outputLevel = intersect(outputLevel, allLevels)

	if(length(outputLevel) < 1) {
		specified = paste(outputLevel, sep="", collapse=", ")
		allowed = paste(allLevels, sep="", collapse=", ")
		msg = paste("ERROR in fitModelAndFormatOutput():",
					"\n\nNo valid output levels specified. Specified levels: ", specified, 
					"\n\nOutput evels must be some combination of: ", allowed,
					"\n\n\n\n", paste = "")

		stop(msg)
	}

	modelTypes = c("isogenicRg", "isogenicGene", "panelGene", "panelRg")

	if(length(outputLevel) < 1) {
		allowed = paste(modelTypes, sep="", collapse=", ")
		msg = paste("ERROR in fitModelAndFormatOutput():",
					"\n\nInvalid model type specified: ", modelType, 
					"\n\nModel type must be one of: ", allowed,
					"\n\n\n\n", paste = "")

		stop(msg)
	}

	# Evaluate model fit with and without the timeNum parameter to evaluate
	# dropout significance
#	f0 = modelFormulaNoCovariate
	f1 = modelFormula

#	m0 = NA
	m1 = NA
	minimalSummary = NA
	finalSummary = NA

	# Weakly informative priors are only supported by the blme package, not lme4, so remove any prior specification
	# in case the model engine is different.
	if(modelEngine != "blme") {
	  prior = NULL
	}

	# While currently disabled, we also have the option of fitting models with and without the covariate
	# and perform a LRT to determine significance.
#	m0 = fitModel(dat, f0, prior)
	m1 = fitModel(dat, f1, modelEngine=modelEngine, prior=prior)

	# Extract and capture any warning/error messages that were returned during model fitting process.
	warningMsg = ifelse(is.null(m1$warning), "", m1$warning$message)
	errorMsg = ifelse(is.null(m1$error), "", m1$error$message)

	if(!(is.null(m1$result)) & !(is.null(m1$result))) {

#		anov01 = fitAnova(m0$result, m1$result)
	  summaryFitMetrics = summaryFitLogLik(m1$result)
	  summaryFixedEffects = summaryFitEffectCoefs(m1$result)

	  # P-value calculations based on the t distribution require providing Denominator Degrees of Freedom (ddf)
	  # values, which are calculated using a heuristic and valid only for the default siMEM model structures
	  # For analyses containing hundreds or thousands of measurements, both calculations produce very similar results
	  # Since the t with high ddf values (> 30 to 50) is well-approximated by a gaussian, at least for the vast bulk
	  # of non-extreme p-values in the analysis, which determine the FDR adjustment used to set significance thresholds. 
	  # The extreme p-values will rank at the head of the list and survive FDR adjustment regardless of approach used.
	  if(pvalueType=="t") {
	    ddf = withinBetweenDF(dat, modelType=modelType)
	    summaryFixedEffectsPValues = getModelPValues(summaryFixedEffects, pvalueType, ddf)
	  } else {
	    summaryFixedEffectsPValues = getModelPValues(summaryFixedEffects, pvalueType)
	  }

#		summaryAOV = summaryAnova(anov01$result)

		finalSummary = cbind.data.frame(summaryFitMetrics,
										summaryFixedEffects,
#										summaryAOV,
										summaryFixedEffectsPValues)

		if("minimal" %in% outputLevel) {
			minimalSummary = getMinimalSummary(finalSummary)
		}
	}

	final = list()

	if("minimal" %in% outputLevel) {
		final[["minimal"]] = minimalSummary
	}

	if("summary" %in% outputLevel) {
		final[["summary"]] = finalSummary
	}

	if("full" %in% outputLevel) {
		final[["full"]] = m1
	}

	final[["warning"]] = warningMsg
	final[["error"]] = errorMsg

	return(final)
}


###
# 
# Fit a single hairpin model for sets of isogenic screens. Isogenic screens are screens performed on the same
# underlying cell line, e.g., cell line +/- drug, or cell line +/- introduced mutation. The presumption being
# that observed differences are due to the perturbation. We don't group by cell line in an isogenic analysis,
# since the underlying cell line is the same.
#
###
fitIsogenicCovariateRg = function(dat, 
                                  formulas=getModelFormulas(), 
                                  pvalueType="gaussian", 
                                  modelEngine="blme",
                                  outputLevel=c("summary")) {

	f1 = formulas$isogenicRg

	# look for the "timeNum" parameter in the model formula to determine whether we're doing an end-point
	# or short time-course analysis
	equationHasTime = grep("timeNum", as.character(f1))
	# If the parameter is not found in the equation, integer(0) is returned, and has length 0
	equationHasTime = length(equationHasTime) > 0

	if(endPoint == TRUE | !equationHasTime) {
	  prior = "gamma"
	} else {
	  prior = "wishart"
	}

	final = fitModelAndFormatOutput(dat=dat,
                									modelFormula=f1,
                									modelType="isogenicRg",
                                  pvalueType=pvalueType,
                									modelEngine=modelEngine,
                									prior=prior,
                									outputLevel=outputLevel)

	return(final)
}


fitIsogenicCovariateGene = function(dat, 
                                    formulas=getModelFormulas(), 
                                    pvalueType="gaussian", 
                                    modelEngine="blme",
                                    outputLevel=c("summary")) {

	f1 = formulas$isogenicGene

	# look for the "timeNum" parameter in the model formula to determine whether we're doing an end-point
	# or short time-course analysis
	equationHasTime = grep("timeNum", as.character(f1))
	# If the parameter is not found in the equation, integer(0) is returned, and has length 0
	equationHasTime = length(equationHasTime) > 0

	if(endPoint == TRUE | !equationHasTime) {
	  prior = "gamma"
	} else {
	  prior = "wishart"
	}

	final = fitModelAndFormatOutput(dat=dat,
                									modelFormula=f1,
                									modelType="isogenicGene",
                                  pvalueType=pvalueType,
                									modelEngine=modelEngine,
                									prior=prior,
                									outputLevel=outputLevel)

	return(final)
}

fitPanelCovariateGene = function(dat, 
                                 formulas=getModelFormulas(), 
                                 pvalueType="gaussian", 
                                 modelEngine="blme",
                                 outputLevel=c("summary"), 
                                 endPoint=FALSE) {

  f1 = formulas$panelGene

  # look for the "timeNum" parameter in the model formula to determine whether we're doing an end-point
  # or short time-course analysis
  equationHasTime = grep("timeNum", as.character(f1))
  # If the parameter is not found in the equation, integer(0) is returned, and has length 0
  equationHasTime = length(equationHasTime) > 0

	if(endPoint==TRUE | !equationHasTime) {
		prior = "gamma"
	} else {
		prior = "wishart"
	}

	final = fitModelAndFormatOutput(dat=dat,
                									modelFormula=f1,
                									modelType="panelGene",
                                  pvalueType=pvalueType, 
                									modelEngine=modelEngine,
                									prior=prior,
                									outputLevel=outputLevel)

	return(final)
}


###
# Fit a single reagent differential essentiality model assayed in multiple cell lines with heterogenous background
# This produces hairpin-level results 
###
fitPanelCovariateRg = function(dat, 
                               formulas=getModelFormulas(), 
                               pvalueType="gaussian", 
                               modelEngine="blme",
                               outputLevel=c("summary"), 
                               endPoint=FALSE) {

  f1 = formulas$panelRg

  # look for the "timeNum" parameter in the model formula to determine whether we're doing an end-point
  # or short time-course analysis
  equationHasTime = grep("timeNum", as.character(f1))
  # If the parameter is not found in the equation, integer(0) is returned, and has length 0
  equationHasTime = length(equationHasTime) > 0

  if(endPoint == TRUE | !equationHasTime) {
    prior = "gamma"
  } else {
    prior = "wishart"
  }

	final = fitModelAndFormatOutput(dat=dat,
                									modelFormula=f1,
                									modelType="panelRg",
                                  pvalueType=pvalueType, 
                									modelEngine=modelEngine,
                									prior=prior,
                									outputLevel=outputLevel)

	return(final)
}


###
# Specify different (simpler) random effects structures to the fixed effects
# and fit all the alternate model structures using the data.
# Pick the best fitting result that didn't produce a model fitting error or warning
# This should eventually use an anova-based step function, like the LMERConvenienceFunctions
# package, but this package seems to add a lot of unnecessary junk to the model formula 
# (redundant ranef statements, etc...) as a result of its forward fitting procedure
###
fitSimplerRandomEffects = function( dat,
                								    fixef1,
                								    randomVector,
                								    modelType,
                								    pvalueType="gaussian", 
                								    modelEngine="blme",
                  									prior="wishart",
                  									outputLevel) {

	formula1 = paste0(fixef1, " + ", randomVector)
	models = cbind.data.frame(formula1, randomEffect=randomVector, stringsAsFactors=F)

	allFits = list()
	allLogL = list()

	for(i in 1:nrow(models)) {

		f1 = as.formula(models[i,"formula1"])

		currFit = fitModelAndFormatOutput(dat=dat,
                											modelFormula=f1,
                											modelType=modelType,
                											modelEngine=modelEngine,
                											prior=prior,
                											outputLevel="summary")

		allFits[[i]] = currFit

		logLikelihood = ifelse(sum(is.na(currFit$summary)) > 0, NA, currFit[["summary"]][1,"logLik"])

		allLogL[[i]] = cbind.data.frame(index=i,
                										logL=logLikelihood,
                										warn=currFit$warning,
                										error=currFit$error,
                										stringsAsFactors=F)
	}

	allLogL = do.call(rbind, allLogL)

	noErrors = allLogL[allLogL$warn == "" & allLogL$error == "", ]

	if(nrow(noErrors) > 0) {
		bestModel = noErrors[which.max(noErrors$logL), "index"]
		final = allFits[[bestModel]]
		randomEffect = models[bestModel, "randomEffect"]
		final$warning = paste0("NOTE: Due to error/warning resulting from full model fit, model with simpler random effect structure used: ", randomEffect)
	} else {
		final = NA
	}

	return(final)
}

##
# Extract a subset of the columns from the full output, namely:
# - basic identifiers
# - for difference(s) in trends: magnitude & pvalue
# - relative dropout rate if available
##
getMinimalSummary = function(modelSummary, calculateRelativeDropout=FALSE, endPoint=FALSE) {

	cols = colnames(modelSummary)
	ids = intersect(c("reagentId", "geneId", "symbol"), cols)

	if(endPoint==TRUE) {
		magnitudes = cols[grep("^.*_estimate$", cols)]
		pvalues = cols[grep("^.*_pvalue$", cols)]

		# Because of end-point results formatting, need to explicitly exclude intercept values from the summary. 
    magnitudes = setdiff(magnitudes, "intercept_estimate")
    pvalues = setdiff(pvalues, "intercept_pvalue")
				
		minimalSummary = modelSummary[,c(ids, magnitudes, pvalues)]

	} else {
		dropouts = cols[grep("^dropout_rate_diff_.*estimate$", cols)]
		pvalues = cols[grep("^dropout_rate_diff_.*pvalue$", cols)]

		if(calculateRelativeDropout) {
		  relatives = cols[grep("^relative_dropout_rate_.*$", cols)]
		  minimalSummary = modelSummary[,c(ids, dropouts, relatives, pvalues)]
		} else {
		  minimalSummary = modelSummary[,c(ids, dropouts, pvalues)]
		}
	}

	return(minimalSummary)
}


##
# Instead of absolute difference, examine the difference relative to dropout rate
# for the category used as comparison baseline
# Ensure that the denominator of the ratio has a minimal non-zero value to prevent
# huge or undefined ratios due to flat dropout rate for the baseline category.
## 
getRelativeDropoutRates = function(results, geneOrReagent="gene") {

#	quant = ifelse(geneOrReagent == "gene", 0.75, 0.5)
	# Use the median of the baseline dropout rates as the minimum
	quant = 0.5
	baselineColumn = "dropout_rate_estimate"
	comparisonColumns = colnames(results)[grep("^dropout_rate_diff_.*_estimate$", colnames(results))]
	relativeDropoutColumns = gsub("^dropout_rate_diff_(.*)_estimate$", "relative_dropout_rate_\\1", comparisonColumns)

	# Calculate relative dropout rates for each non-baseline category and add them
	# to the results matrix
	for(i in 1:length(comparisonColumns)) {
		results[,relativeDropoutColumns[i]] = getRelativeDropout(results[,baselineColumn], results[,comparisonColumns[i]], quant=quant)
	}

	return(results)
}


getRelativeDropout = function(baseline, delta, quant=0.5) {

	if(quant < 0.01 | quant > 0.99) {
		quant = 0.5
	}

	# A minimum value for the ratio denominator so that the ratio doesn't blow up
	dropoutFloor = abs(quantile(baseline, probs=quant, na.rm=TRUE))
	# More or less lethal?
	deltaSign = sign(delta)

	alternative = baseline + delta

	minMax = t(apply(cbind(baseline, alternative), 1, sort))
	# Without proper NA handling in the sort function, rows with NAs are removed
	# and the casting to matrix/DF fails, thus breaking the code.
	minMax = t(apply(cbind(baseline, alternative), 1, sort, na.last=T))
	minMax = as.data.frame(minMax)
	colnames(minMax) = c("numerator", "denominator")

	# If the the greater numeric trend is above 0, shift both trends down so
	# that the greater one is now at -dropoutFloor.
	shifts = ifelse(minMax$denominator > -dropoutFloor, -(minMax$denominator+dropoutFloor), 0)
	minMax$numerator = minMax$numerator + shifts
	minMax$denominator = minMax$denominator + shifts
	# Set minimal magnitude of ratio denominator to avoid extreme/undefined ratios
  #	minMax$denominator = ifelse(minMax$denominator > dropoutFloor, dropoutFloor, minMax$denominator)

	# Negative relative dropouts are more essential in the second condition
	# whereas positive relative dropouts are more essential in the baseline condition
	ratio = deltaSign * abs(minMax$numerator/minMax$denominator)
	return(ratio)
}


##
# Given the t-values, and calculated degrees of freedom, generate pvalues for
# model fixed effects.
#
# For the t-based p-value calculation, we need to specify the Denominator Degrees of Freedom
# which is specified for the siMEM model structure using the Between-Within heuristic
# (as implemented in the 'nlme' R package)
##
getModelPValues = function(fixedEffects, pvalueType="gaussian", ddf=NULL) {

  # Get the t-statistics output by the model
	summaryT = unlist(fixedEffects[,grep("_t$", colnames(fixedEffects)), drop=T])
	cols = names(summaryT)
	# Generate column names for p-values
	colsPvalues = gsub("_t$", "_pvalue", cols)

	# For the current siMEM model structure, all 3 parameter (intercept, baseline slope, slope difference) ddfs are identical
	# This will not be the case, for example, if we have a 4 parameter (+ COVARIATE) model
	pvalueCount = length(cols)

	if(pvalueType=="t") {
	  colsDF = gsub("_t$", "_df", cols)
	  commonDF = ddf[["intercept"]]
	  matrixDF = rep(commonDF, pvalueCount)
	  pvalues = 2*pt(-abs(summaryT), df=commonDF)
	  pvalues = data.frame(t(c(matrixDF, pvalues)))
	  colnames(pvalues) = c(colsDF, colsPvalues)
	} else {
	  # Here we Rely on the Gaussian approximation of the t distribution as ddf increases past 30-50
	  pvalues = 2*pnorm(-abs(summaryT))
	  pvalues = data.frame(t(pvalues))
	  colnames(pvalues) = colsPvalues
	}

	return(pvalues)
}

##
# Apply Benjamini-Hochberg adjustment to each P-value column in the results
##
applyFDR = function(results) {

	if(is.data.frame(results) && nrow(results) > 1) {
		pvals = results[,grep("pvalue", colnames(results)),drop=F]
		cols = colnames(pvals)
		colsFDR = gsub("pvalue", "fdr", cols)
		fdr = apply(pvals, 2, p.adjust, method="BH")
		colnames(fdr) = colsFDR
		results = cbind.data.frame(results, fdr, stringsAsFactors=F)
	}
	
	return(results)
}

##
# given a list of T values, and associated degrees of freedom, calculate the
# associated pvalues 
##
pvaluesT = function(tvalues, dfs) {

	# if for whatever reason, the nlme call failed, use the smallest degrees of
 # freedom (most "conservative") to get pvalue associated with each T value
	mindf = min(dfs)
	dfs = ifelse(dfs == -1, mindf, dfs)
	tvalues = -abs(tvalues)

	vect = cbind.data.frame(tvalues, dfs)
	pvalues = apply(vect, 1, function(v){2*pt(v[1], df=v[2])})
	return(pvalues)	
}

##
# For larger T degrees of freedom (Rule of Thumb, T DoF > 30-50)
# T distribution is well-approximated by Gaussian
# so T statistic is treated as Z statistic and p-value computed using Gaussian
# approximation to T distribution with large degrees of Freedom
#
# For a panel of 20+ screens
# Estimated degrees of freedom for panel analysis range in the hundreds to thousands
# so large-sample approximation is applicable to panel analysis, either at reagent or gene level
## 
pvaluesNorm = function(tvalues) {
	tvalues = -abs(tvalues)
	pvalues = 2*pnorm(tvalues)
	return(pvalues)
}



getVarianceWeights = function(dat, min_var = 0.01, max_var=2) {

	dat$grouping = paste(dat$reagentId, dat$cellLineId, dat$timeGroup, dat$COVARIATE, sep="-")
	vars = by(dat$expr, list(dat$grouping), var, na.rm=TRUE)

	vars = cbind.data.frame(names(vars), as.numeric(vars), stringsAsFactors=F)
	colnames(vars) = c("grouping", "var")
	dat = merge(dat, vars)
	dat$var_bounded = dat$var
	dat$var_bounded = ifelse(dat$var_bounded < min_var, min_var, dat$var_bounded)
	dat$var_bounded = ifelse(dat$var_bounded > max_var, max_var, dat$var_bounded)

	# Weights are inversely proportional to the variance, multiplied by a scaling
	# factors s.t. the total sum of weights is still the number of rows/data points
	# without this, the significance values change substantially.
	dat$W = (nrow(dat)/sum(1/dat$var_bounded, na.rm=TRUE)) * (1/dat$var_bounded)
	return(dat)
}


##
# For the siMEM model structures, calculate the Denominator Degrees of Freedom
# using the Between-Within heuristic
#
# DDF = m_i - m_(i-1) - p_i
# by convention, the DDF of the intercept is evaluated at the lowest level
# the time and time:COVARIATE covariates are also evaluated at this level, since there
# are multiple -different- timed measurements associated with each screen, hence
# time and time:COVARIATE are "WITHIN/INNER" to the screenId grouping level
#
##
withinBetweenDF = function(dat, modelType="panelGene") {

	ddf = list()

	if(modelType == "panelGene") {

		denDFi = nrow(dat) - length(unique(dat$cellLineIdUnique)) - length(unique(dat$COVARIATE))

		# DDF = m_i - m_(i-1) - p_i
		# For Panel analyses, the COVARIATE is a cell-line specific variable, hence changes at
		# the level of cellLineId, but not within replicate measures of a given cellLineId
		denDFdelta = length(unique(dat$cellLineIdUnique)) - length(unique(dat$reagentId)) - length(unique(dat$COVARIATE)) + 1

		ddf[["intercept"]] = denDFi
		ddf[["interceptdelta"]] = denDFdelta
		ddf[["timenum"]] = denDFi
	} else if(modelType == "panelRg") {
		denDFi = nrow(dat) - length(unique(dat$cellLineIdUnique)) - length(unique(dat$COVARIATE))
		# Here, the m_(i-1) = m_0 = 1 by convention, since the model includes an intercept term.
		denDFdelta = length(unique(dat$cellLineIdUnique)) - 1 - length(unique(dat$COVARIATE)) + 1
		ddf[["intercept"]] = denDFi
		ddf[["interceptdelta"]] = denDFdelta
		ddf[["timenum"]] = denDFi

	} else if(modelType == "isogenicGene") {
		# DDF = m_i - m_(i-1) - p_i
		# by convention, the DDF of the intercept is evaluated at the lowest level
		# the time and time:COVARIATE covariates are also evaluated at this level, since there
		# are multiple -different- array expression values associated with each screen, hence
		# time and time:COVARIATE are "WITHIN/INNER" to the screenId grouping level
		denDFi = nrow(dat) - length(unique(dat$replicateIdUnique)) - length(unique(dat$COVARIATE))

		# For isogenic screens, the COVARIATE variable itself changes at the level of replicateIdUnique, 
		# so the DDF is the same as for the intercept
		denDFdelta = denDFi

		ddf[["intercept"]] = denDFi
		ddf[["interceptdelta"]] = denDFdelta
		ddf[["timenum"]] = denDFi

	} else if(modelType == "isogenicRg") {
		# DDF = m_i - m_(i-1) - p_i
		# by convention, the DDF of the intercept is evaluated at the lowest level
		# the time and time:COVARIATE covariates are also evaluated at this level, since there
		# are multiple -different- array expression values associated with each screen, hence
		# time and time:COVARIATE are "WITHIN/INNER" to the screenId grouping level
		denDFi = nrow(dat) - length(unique(dat$replicateIdUnique)) - length(unique(dat$COVARIATE))

		# For isogenic screens, the COVARIATE variable itself changes at the level of replicateIdUnique, 
		# so the DDF is the same as for the intercept
		denDFdelta = denDFi

		ddf[["intercept"]] = denDFi
		ddf[["interceptdelta"]] = denDFdelta
		ddf[["timenum"]] = denDFi
	} else {
		msg = paste("ERROR in withinBetweenDF():",
					"\n\nUnrecognized model type Between-Within Degrees of Freedom calculation: ", modelType, sep="")
		stop(msg)
	}

	return(ddf)
}

##
# Get per-row mean-variance for a set of replicates (eg: same screen, timepoint)
##
getReplicateMeanVar = function(reps) {

	if(ncol(reps) > 1) {
		mu = rowMeans(reps, na.rm=TRUE)
		variance = rowVars(reps, na.rm=T)
	} else {
	  # Since variance is used to weight measurements, setting variance to 1 by default
	  # will lead to all measurements being equally weighted
		mu = reps[,1]
		variance = rep(1, length(mu))
	}

	final = cbind.data.frame(mu=mu, variance=variance)
	return(final)
}

##
# smoothing with locfit causes segfault in locfit C code for some data (Project Achilles Cowley 2014 in particular)
##
getMeanVar = function(exprs, pheno, minVar=0.1, smoothing = 0.5, method="loess", variableMap=getDefaultVariableMap(), printProgress=FALSE) {

	meanVar = list()
	meanVarFit = list()
	meanVarPredicted = list()

	message("Beginning calculation of replicate means, and smoothed mean-variance function for grouped/replicate samples...")

	# Identify the column that identifies biological replicates for the mean-variance
	# calculation, and calculate the mean-variance function empirically, using local regression
	for(grp in unique(pheno[,variableMap$replicateGroupId])) {

	  if(printProgress) {
	    message(grp)
	  }

		# Get the unique sample identifiers for all cell line/timepoint replicates
		samples = pheno[pheno[,variableMap$replicateGroupId] == grp, variableMap$sampleId]
		sampleDat = exprs[,samples,drop=F]

		meanVarTable = getReplicateMeanVar(sampleDat)

		if(length(samples) > 1) {

		  if(method == "locfit") {
		    # the locfit C code sometimes segfaults - possibly when too many mean-variance pairs are identical
		    meanVarFitted = locfit(variance ~ mu, data=meanVarTable, family="gamma", alpha=smoothing)
		    # To get the precision (inverse-variance) weights, first get the smoothed mean-variance fit at the given mean.
		    meanVarTable$variance_smoothed = fitted(meanVarFitted, data=meanVarTable)
      } else {
        meanVarFitted = loess(variance ~ mu, data=meanVarTable, span=smoothing)
        # To get the precision (inverse-variance) weights, first get the smoothed mean-variance fit at the given mean.
        meanVarTable$variance_smoothed = predict(meanVarFitted, data=meanVarTable)
      }

		} else {
			# In the case of a single replicate, set the mean as the values of that replicate, and the variance to 1
			# inverse-variance weights will then be equal to 1, ie: equal weights.
			meanVarFitted = NA
			meanVarTable$variance_smoothed = meanVarTable$variance
		}

		meanVarTable$mu = round(meanVarTable$mu, digits=2)
		meanVarTable$variance = round(meanVarTable$variance, digits=3)
		meanVarTable$variance_smoothed = round(meanVarTable$variance_smoothed, digits=3)
		meanVarTable$variance_smoothed = ifelse(meanVarTable$variance_smoothed < minVar, minVar, meanVarTable$variance_smoothed)

		meanVar[[grp]] = meanVarTable
		meanVarFit[[grp]] = meanVarFitted
	}

	final = list()
	final[["mean_var"]] = meanVar
	final[["mean_var_fit"]] = meanVarFit
	return(final)
}

##
# Calculate smoothed mean-variance curves for each set of grouped/replicate samples, and
# assign inverse-variance weights to each data data point in the matrix
##
addPrecisionWeights = function(eset, variableMap = getDefaultVariableMap(), minVar=0.01, smoothing=0.05, printProgress=FALSE) {

	shrna = exprs(eset)
	pheno = pData(eset)
	fdat = fData(eset)

	meanVar = getMeanVar(shrna, pheno, minVar=minVar, smoothing=smoothing, variableMap=variableMap, printProgress=printProgress)

	message("Finalizing mean-variance function calculations and assigning precision weights to data points...")

	groupMeans = lapply(meanVar[[1]], function(dataf) {return(dataf[,"mu"])} )
	groupMeans = t(do.call(rbind, groupMeans))
	groupMeans = groupMeans[,pheno[,variableMap$replicateGroupId]]
	colnames(groupMeans) = pheno[,variableMap$sampleId]
	rownames(groupMeans) = rownames(shrna)

	groupVars = lapply(meanVar[[1]], function(dataf) {return(dataf[,"variance"])} )
	groupVars = t(do.call(rbind, groupVars))
	groupVars = groupVars[,pheno[,variableMap$replicateGroupId]]
	colnames(groupVars) = pheno[,variableMap$sampleId]
	rownames(groupVars) = rownames(shrna)

	varWeights = lapply(meanVar[[1]], function(dataf) {return(dataf[,"variance_smoothed"])} )
	varWeights = t(do.call(rbind, varWeights))
	# duplicate the weights for columns containing replicates for the same timepoint
	varWeights = varWeights[,pheno[,variableMap$replicateGroupId]]
	colnames(varWeights) = pheno[,variableMap$sampleId]
	rownames(varWeights) = rownames(shrna)
	varWeights = 1/varWeights
	varWeights = round(varWeights, digits=2)

	# Add the means and smoothed precision weights to the ExpressionSet in additional assayData slots
	# These slots will be accessed if the models are fit with precision weights to mitigate heteroscedasticity.
	storageMode(eset) = "environment"
	assayData(eset)[["varWeights"]] = varWeights
	assayData(eset)[["groupMeans"]] = groupMeans
	assayData(eset)[["groupVars"]] = groupVars
	storageMode(eset) = "lockedEnvironment"

	message("Done.")
	return(eset)
}

##
# Fit a local regression to estimate the posterior probability of "signal"
# as a function of measurement intensity
##
getSignalProb = function(signalProb) {
	# Use the predict method to interpolate posterior probability for values not covered in
	# input signalProb table, which has a 
	fit = splinefun(signalProb[,1], signalProb[,2], method="natural")
	return(fit)
}


addSignalProbWeights = function(eset, signalProbTable, variableMap=getDefaultVariableMap()) {

	posterior = getSignalProb(signalProbTable)

	groupMeans = assayData(eset)[["groupMeans"]]
	pheno = pData(eset)

	# This function calculates both group means and the variable map
	if(is.null(groupMeans)) {
		eset = addMeanVarianceWeights(eset, variableMap=variableMap)
	}

	if(length(grep(variableMap$screenId, colnames(pheno))) == 0) {
		pheno[,variableMap$screenId] = paste(pheno[,variableMap$cellLineId], pheno[,variableMap$condition], sep="_")
		pheno[,variableMap$screenId] = gsub("_$", "", pheno[,variableMap$screenId])
		print(table(pheno[,variableMap$screenId]))
	}

	# Order the samples by sample id, then timepoint
	meta = pheno[order(pheno[,variableMap$screenId], pheno[,variableMap$timeGroup]),]

	tempMx = matrix(1, nrow=nrow(groupMeans), ncol=ncol(groupMeans))
	rownames(tempMx) = rownames(groupMeans)
	colnames(tempMx) = colnames(groupMeans)

	message("Begin calculating signal probability weights...")

	for(screenId in unique(meta[,variableMap$screenId])) {

		message(screenId)

		# Reset screen-specific intermediate variables
		posteriorProd = rep(1, nrow(groupMeans))
		previousProb = rep(1, nrow(groupMeans))
	#	previousProb = 1:nrow(groupMeans)

		timePoints = unique(meta[ meta[,variableMap$screenId] == screenId, variableMap$timeGroup])
	
		for(timept in timePoints) {

			samples = meta[ meta[,variableMap$screenId] == screenId & meta[,variableMap$timeGroup] == timept, variableMap$sampleId]
	
			repeated = matrix(rep(previousProb, times=length(samples)),
								nrow=length(previousProb), ncol=length(samples))
			colnames(repeated) = samples

			# Store signal weights for this timepoint
			tempMx[,samples] = repeated
	
			# Get the posterior signal using the mean of all replicates
			# The mean is repeated for all replicate samples, so just take 1st one
			posteriorProb = posterior(groupMeans[,samples[1],drop=T])

			# Keep a running product of probabilities from all previous timepoints
			previousProb = previousProb * posteriorProb
			previousProb = ifelse(previousProb < 0.01, 0.01, previousProb)
		}

	#	stop("...")
	}

	# Ensure that the rows and columns match, so we don't incorrectly map weights
#	table(rownames(tempMx) == rownames(groupMeans))
#	table(colnames(tempMx) == colnames(groupMeans))

	message("Finalizing signal probability weights calculations ...")

	tempMx = round(tempMx, digits=3)

	storageMode(eset) = "environment"
	assayData(eset)[["signalProbWeights"]] = tempMx
	storageMode(eset) = "lockedEnvironment"
	return(eset)
}

###
# Store/access model formulas for later retrieval.
# The list returned by this function can be passed as a parameter to simem(),
# making it easier to run models with alternate structures.
###
getModelFormulas = function(type="") {

	# Gene-level models

	# For end-point analyses
	# f1 = as.formula("expr ~ 1 + COVARIATE + (1|reagentId/cellLineId)")
  # CROSSED random effects
  # f1 = as.formula("expr ~ 1 + COVARIATE + (1|reagentId) + (1|cellLineId)")
  # Reagent nested in cell line
	# f1 = as.formula("expr ~ 1 + COVARIATE + (1|cellLineId/reagentId)")

	# For short time-course analyses
	# f1 = as.formula("expr ~ 1 + timeNum + timeNum:COVARIATE + (timeNum|reagentId/cellLineId)")
  # CROSSED random effects
	# f1 = as.formula("expr ~ 1 + timeNum + timeNum:COVARIATE + (timeNum|reagentId) + (timeNum|cellLineId)")
  # Reagent nested in cell line
	# f1 = as.formula("expr ~ 1 + timeNum + timeNum:COVARIATE + (timeNum|cellLineId/reagentId)")

	modelMap = list(
		panelGene = as.formula("expr ~ 1 + timeNum + timeNum:COVARIATE + (timeNum|reagentId/cellLineId)"),
		panelRg = as.formula("expr ~ 1 + timeNum + timeNum:COVARIATE + (timeNum|cellLineId)"),
		isogenicGene = as.formula("expr ~ 1 + timeNum + timeNum:COVARIATE + (timeNum|reagentId) + (0+timeNum|reagentId/replicateId)"),
		isogenicRg = as.formula("expr ~ 1 + timeNum + timeNum:COVARIATE + (0+timeNum|replicateId)")
	)

	if(type == "endpoint") {
	  modelMap$panelGene = as.formula("expr ~ 1 + COVARIATE + (1|reagentId/cellLineId)")
	  modelMap$panelRg = as.formula("expr ~ 1 + COVARIATE + (1|cellLineId)")
	}

	return(modelMap)
}
