#' @title Check on the database input for simulating the model on a database of values
#'
#' @description Checks the presence of required variables and adds default values for missing variables.
#'
#' @param df a dataframe containing the data, with one column called h containing the daily incidence
#' and one variable called prop_import containing the proportion of imported cases among new infections.
#' additional optional variables are alpha (effective care), beta (proportion of liver stage cure),
#' rho (reporting rate) and omega (intensity of vector control)
#' @param delay a boolean indicating if the model including delays in treatment should be used.
#' @param rcd a boolean indicating if the model including reactive case detection should be used.
#' @param mda a boolean indicating if the model including mass drug administration (MDA) prophylaxis should be used.
#' @param sto a boolean indicating if the stochastic version of the model should be used.
#' @param rcd_at_baseline a boolean indicating if the model was calibrated using the RCD model (i.e. there is some RCD at baseline already). Default (FALSE) is the model without RCD at baseline
#'
#' @return A dataframe containing all required variables
#'
sanity_checks_inputs_simulate=function(df, delay, rcd, mda, sto, rcd_at_baseline){

  ####################################
  # sanity checks on the inputs
  if(!"id" %in% names(df)){ stop("no id variable in df")}
  if(!"lambda" %in% names(df)){ stop("no lambda variable in df")}
  if(!"rho.old" %in% names(df)){
    df$rho.old=1
    warning("no rho.old in df, assumed rho.old=1")
  }
  if(!"omega.old" %in% names(df)){
    df$omega.old=1
    warning("no omega.old in df, assumed omega.old=1")
  }
  if(!"alpha.old" %in% names(df)){
    df$alpha.old=0
    warning("no alpha.old in df, assumed alpha.old=0")
  }
  if(!"beta.old" %in% names(df)){
    df$beta.old=1
    warning("no beta.old in df, assumed beta.old=1")
  }
  if(!"sigma.old" %in% names(df) & delay==TRUE){
    stop("no sigma (sigma.old) variable in df, maybe use the model without delay (calculate_r0_rc_fromdata)")
  }

  if(!"rho.new" %in% names(df)){
    df$rho.new=df$rho.old
    warning("no rho.new in df, assumed rho.new=rho.old")
  }
  if(!"omega.new" %in% names(df)){
    df$omega.new=df$omega.old
    warning("no omega.new in df, assumed omega.new=omega.old")
  }
  if(!"alpha.new" %in% names(df)){
    df$alpha.new=df$alpha.old
    warning("no alpha.new in df, assumed alpha.new=alpha.old")
  }
  if(!"beta.new" %in% names(df)){
    df$beta.new=df$beta.old
    warning("no beta.new in df, assumed beta.new=beta.old")
  }
  if(!"sigma.new" %in% names(df) & delay==TRUE){
    df$sigma.new=df$sigma.old
    warning("no sigma.new in df, assumed sigma.new=sigma.old")
  }
  if(!"delta" %in% names(df)){
    df$delta=0
    warning("no delta in df, assumed delta=0")
  }
  if(!"delta.new" %in% names(df)){
    df$delta.new=df$delta
    warning("no delta.new in df, assumed delta.new=delta")
  }

  if((!"iota.new" %in% names(df) | !"tau.new" %in% names(df) | !"nu.new" %in% names(df) | !"eta.new" %in% names(df)| !"kappa.new" %in% names(df)) & rcd==TRUE){
    stop("RCD parameters are missing, please add them or use model without RCD")
  }

  if((!"iota_star" %in% names(df) | !"tau.old" %in% names(df) | !"nu.old" %in% names(df) | !"eta.old" %in% names(df)| !"kappa.old" %in% names(df)) & rcd_at_baseline==TRUE){
    stop("RCD parameters are missing at baseline, please add them or use model without RCD (rcd_at_baseline=F)")
  }

  if((!"MDAcov.new" %in% names(df) | !"MDAp_length.new" %in% names(df) | !"MDArad_cure.new" %in% names(df)  ) & mda==TRUE){
    stop("MDA parameters are missing, please add them or use model without MDA")
  }

  if((!"N" %in% names(df) ) & sto==TRUE){
    stop("population size N is missing for stochastic model. Please add it or use deterministic model")
  }

  return(df)
}


#' @title Check on the database input for calculating lambda at baseline
#'
#' @description Checks the presence of required variables and adds default values for missing variables.
#'
#' @param df a dataframe containing the data, with one column called h containing the daily incidence
#' and one variable called prop_import containing the proportion of imported cases among new infections.
#' additional optional variables are alpha (effective care), beta (proportion of liver stage cure),
#' rho (reporting rate) and omega (intensity of vector control)
#' @param delay a boolean indicating if the model including delays in treatment should be used.
#' @param rcd_at_baseline a boolean indicating if the model was calibrated using the RCD model (i.e. there is some RCD at baseline already). Default (FALSE) is the model without RCD at baseline
#'
#' @return A dataframe containing all required variables
#'
sanity_checks_inputs_calculate=function(df, delay, rcd_at_baseline){
  if(!"h" %in% names(df)){ stop("no h variable in df")}
  if(!"prop_import" %in% names(df)){ stop("no prop_import variable in df")}


  if(!"sigma" %in% names(df) & delay==TRUE){ stop("no sigma variable in df, maybe use the model without delay (calculate_r0_rc_fromdata)")}
  if(!"rho" %in% names(df)){
    df$rho=1
    warning("no rho in df, assumed rho=1")
  }
  if(!"omega" %in% names(df)){
    df$omega=1
    warning("no omega in df, assumed omega=1")
  }
  if(!"alpha" %in% names(df)){
    df$alpha=0
    warning("no alpha in df, assumed alpha=0")
  }
  if(!"beta" %in% names(df)){
    df$beta=1
    warning("no beta in df, assumed beta=1")
  }

  if((!"iota" %in% names(df) | !"tau" %in% names(df) | !"nu" %in% names(df) | !"eta" %in% names(df)| !"kappa" %in% names(df)) & rcd_at_baseline==TRUE ){
    stop("RCD parameters are missing, please add them or use model without RCD")
  }
  return(df)
}
