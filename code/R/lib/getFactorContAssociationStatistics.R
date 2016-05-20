# Find Inter Class Correlation between factor and continuous covariates
# Inspired from http://stats.stackexchange.com/questions/108007/correlations-with-categorical-variables
getFactorContAssociationStatistics <- function(factorContNames,COVARIATES, na.action='remove', 
                                               alpha = 0.05){
  require(psych)
  factorContNames = factorContNames %>% unlist %>% as.character()
  
  if (na.action == "remove")
    COVARIATES = na.omit(COVARIATES[,factorContNames])
  
  COVARIATES = sapply(COVARIATES, as.numeric)
  stats = psych::ICC(COVARIATES, missing = T, alpha = alpha)
  
  return(data.frame(Estimate = stats$results['Single_fixed_raters','ICC'],
                    Pval = stats$results['Single_fixed_raters','p']))
}
