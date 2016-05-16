# Function to calculate correlation and plot
calcCompleteCorAndPlot <- function(COMPARE_data, COVAR_data, correlationType, title, 
                                   WEIGHTS = NULL, PLOT_ALL_COVARS=FALSE, EXCLUDE_VARS_FROM_FDR=NULL, MAX_FDR = 0.1) {
  
  # require(plyr)
  
  # Get factor and continuous covariates
  FactorCovariates <- colnames(COVAR_data)[sapply(COVAR_data,is.factor)]
  ContCovariates <- setdiff(colnames(COVAR_data),FactorCovariates)
  
  # Calculate correlation between compare_data and factor covariates
  if (length(FactorCovariates) > 0){
    comb <- expand.grid(colnames(COMPARE_data),FactorCovariates)
    factCont_cor <- plyr::ddply(comb,
                                .(Var1, Var2),
                                .fun = getFactorContAssociationStatistics,
                                cbind(COMPARE_data,COVAR_data[rownames(COMPARE_data),FactorCovariates,drop=F]),
                                alpha=MAX_FDR)
    factCont_cor_vals <- factCont_cor %>%
      dplyr::select(Var1, Var2, Estimate) %>%
      tidyr::spread(Var2, Estimate)
    rownames(factCont_cor_vals) <- factCont_cor_vals$Var1
    factCont_cor_vals$Var1 <- NULL
    factCont_cor_vals = factCont_cor_vals %>% data.matrix()
    
    factCont_cor_p <- factCont_cor %>%
      dplyr::select(Var1, Var2, Pval) %>%
      tidyr::spread(Var2, Pval)
    rownames(factCont_cor_p) <- factCont_cor_p$Var1
    factCont_cor_p$Var1 <- NULL
    factCont_cor_p = factCont_cor_p %>% data.matrix()
    
  } else {
    factCont_cor_vals <- NULL
    factCont_cor_p <- NULL
  }
  
  # Calculate correlation between compare_data and factor covariates
  if (length(ContCovariates) > 0){
    cont_cor <- corr.test(COMPARE_data,
                          COVAR_data[,ContCovariates, drop=F],
                          use='pairwise.complete.obs',
                          method=correlationType, 
                          adjust="none")
    cont_cor_vals <- cont_cor$r
    cont_cor_p <- cont_cor$p
    
    rownames(cont_cor_vals) <- colnames(COMPARE_data)
    colnames(cont_cor_vals) <- ContCovariates
    
    rownames(cont_cor_p) <- colnames(COMPARE_data)
    colnames(cont_cor_p) <- ContCovariates
  } else {
    cont_cor_vals <- NULL
    cont_cor_p <- NULL
  }
  
  all_cor_vals = cbind(factCont_cor_vals,cont_cor_vals)
  all_cor_p = cbind(factCont_cor_p,cont_cor_p)
  
  Effects.significantCovars = all_cor_vals
  Effects.significantCovars[all_cor_p>MAX_FDR] = 0
  Effects.significantCovars = colSums(abs(Effects.significantCovars)*replicate(dim(Effects.significantCovars)[2],WEIGHTS/sum(WEIGHTS)))
  Effects.significantCovars = Effects.significantCovars[order(abs(Effects.significantCovars),decreasing=T)]
    
  cor_mat = melt(all_cor_p, varnames=c("COMPARE", "COVAR"))
  colnames(cor_mat)[colnames(cor_mat) == "value"] = "pvalue"
  
  cor_mat$COMPARE = factor(cor_mat$COMPARE, levels=rownames(all_cor_p))
  cor_mat$COVAR = factor(cor_mat$COVAR, levels=colnames(all_cor_p))
  
  cor_mat$r = melt(all_cor_vals)$value
  
  calcFDRrows = rep(TRUE, nrow(cor_mat))
  markColumnsAsMissing = NULL
  if (!is.null(EXCLUDE_VARS_FROM_FDR)) {
    calcFDRrows = !(cor_mat$COVAR %in% EXCLUDE_VARS_FROM_FDR)
    markColumnsAsMissing = intersect(colnames(COVAR_data), EXCLUDE_VARS_FROM_FDR)
  }
  
  # Entries that pass the threshold of "significance":  
  markSignificantCorrelations = corMatFDRthreshFunc(cor_mat, indicesMask=calcFDRrows, MAX_FDR = 0.1)
  significantCorrelatedCovars = sort(unique(cor_mat$COVAR[markSignificantCorrelations]))
  
  markPotentialSignificantCorrelations = corMatFDRthreshFunc(cor_mat)
  # Specially mark only those incomplete covariates that would be significant in the context of all covariates:
  markPotentialSignificantCorrelations = markPotentialSignificantCorrelations & !calcFDRrows
  
  plotRows = 1:nrow(cor_mat)
  if (!PLOT_ALL_COVARS) {
    # Plot all correlations for:
    # 1) Covariates with at least one significant correlation
    # 2) Excluded covariates
    plotRows = (cor_mat$COVAR %in% significantCorrelatedCovars) | !calcFDRrows
  }
  plotCor = na.omit(cor_mat[plotRows, ])
  
  for (markCor in c("markSignificantCorrelations", "markPotentialSignificantCorrelations")) {
    useMarkCor = get(markCor)[plotRows]
    if (length(which(useMarkCor)) > 0) {
      plotCor[, markCor] = useMarkCor[ setdiff(1:length(useMarkCor), as.numeric(attr(plotCor, "na.action"))) ]
    }
  }
  
  if (!plyr::empty(plotCor)){
    plot = plotCorWithCompare(plotCor, title, paste("FDR <= ", MAX_FDR, sep=""), markColumnsAsMissing)
  } else{
    plot = NULL
  }
  
  return(list(plot=plot, significantCovars=as.character(significantCorrelatedCovars), Effects.significantCovars = Effects.significantCovars))
}