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




#' @title Creates linear interpolation
#'
#' @param x time points
#' @param y values to interpolate
#' @param start rescales the x to x-start (useful to chain interventions properly)
#'
#' @return An interpolation function, like the one created by approxfun (stats package)
#' @export
create_approximation=function(x, y, start=0){
  return(stats::approxfun(x=x-start, y=y, method = "linear"))
}



#' @title Creates vector control exponential decay
#'
#' @param initial_omega the reduction in vectorial capacity at time=0
#' @param half_life half life of the exponential decay (in years)
#' @param every_x_years frequency at which the intervention is renewed (default is every 3 years)
#' @param maxtime maximum time point to be included in the database (from 0 to maxtime days)
#'
#' @return An interpolation function, like the one created by approxfun (stats package)
#' @export
vector_control_exponential_decay=function(initial_omega, half_life, every_x_years=3, maxtime){
  my_omega=data.frame(time=seq(0, maxtime))
  my_omega$value=1-(1-initial_omega)*exp(-log(2)*((my_omega$time/365)%%(every_x_years))/half_life)
  return(my_omega)
}




#' @title Function for time varying tau
#' @description From Chitnis et al. 2019 (supplementary information)
#' @param nu time
#' @param pr prevalence
#' @param N population
#'
#' @return A scalar value for tau
#' @export
varying_tau=function(nu, pr, N=10000){
  a1=0.23
  a2=-1.4
  a3=2.87
  return( exp((-a1*log(pr)+a2/nu-a3*log(pr)/nu)*(N-nu)/N))
}
