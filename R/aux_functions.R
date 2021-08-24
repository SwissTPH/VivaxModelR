#' @title Converts annual incidence to daily
#'
#' @param incidence annual incidence (cases per 1000 person year)
#'
#' @return Scalar representing daily incidence (cases per 1 person day).
#' h=incidence/1000/365
#' @export
#'
incidence_year2day=function(incidence){
  return(incidence/365/1000)
}


#' @title Converts daily incidence to annual
#'
#' @param h daily incidence (cases per 1 person day)
#'
#' @return Scalar representing annual incidence (cases per 1000 person year)
#' incidence=h*365*1000
#' @export
#'
incidence_day2year=function(h){
  return(h*365*1000)
}




#' @title Samples uncertainty in incidence and importation data
#'
#' @param df dataframe containing the variables "cases" (reported cases),
#' "cases_local" (number of reported cases) and "population" (population size)
#' @param ndraw number of samples
#'
#' @return A dataframe with ndraw replicates of the original dataframe df, containing 4 additional variables.\cr
#' h is the sampled daily incidence (cases per person day), \cr
#' prop_import is the sampled proportion of imported cases
#' @export
#'
sample_uncertainty_incidence_import=function(df, ndraw=100){

  if(!"cases" %in% names(df)){ stop("no cases variable in df")}
  if(!"cases_local" %in% names(df)){ stop("no cases_local variable in df")}
  if(!"population" %in% names(df)){ stop("no population variable in df")}

  if(any(df$cases<df$cases_local)){ stop("error: cases_local > cases")}

  this.data=do.call("rbind", replicate(ndraw, df,
                                       simplify = FALSE))
  this.data$h=stats::rgamma(nrow(df)*ndraw,shape = (1/2)+this.data$cases, rate=this.data$population)/365
  this.data$prop_import=stats::rbeta(nrow(df)*ndraw, (1/2)+this.data$cases-this.data$cases_local, (1/2)+this.data$cases_local)

  return(this.data)
}
