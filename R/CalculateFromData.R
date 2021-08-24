#' @title Calculate R0 and Rc on an incidence dataset
#'
#' @description Uses the compartmental model to calculate R0 and Rc using data on
#' incidence and the proportion of imported cases
#'
#' @param df a dataframe containing the data, with one column called h containing the daily incidence
#' and one variable called prop_import containing the proportion of imported cases among new infections.
#' additional optional variables are alpha (effective care), beta (proportion of liver stage cure),
#' rho (reporting rate) and omega (intensity of vector control)
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param return.all if TRUE, also returns delta and I estimates
#' @param h.cutoff R0/Rc are not calculated for h values which are strictly inferior to h.cutoff
#'
#' @return A dataframe with 3 additional columns: lambda is the transmission rate,
#' R0 is the basic reproduction number, Rc is the controlled reproduction number. \cr
#' If lambda = -2, it means that the parameter combination does not correspond to a positive lambda solution with the model. \cr
#' If lambda = -3, it means the incidence data contains missing values or infinite values or values below h.cutoff.  \cr
#' If return.all=T, delta and I are also included.
#'
#'
#' @details If alpha is not provided in df, alpha=0. If beta is not provided in df, beta=1.
#' If rho is not provided in df, rho=1. If omega is not provided in df, omega=1.
#'
#' @example
#' mydata=data.frame(incidence=c(23,112),prop_import=c(0,0.1))
#' mydata$h=incidence_year2day(mydata$incidence)
#' mydata$beta=c(0.43,0.42)
#' mydata$alpha=c(0.17, 0.12)
#' mydata$rho=c(0.17, 0.12)
#' mydata$omega=c(1,1)
#' calculate_r0_rc_fromdata(df=mydata, f=1/72, gamma=1/223, r=1/60, return.all = TRUE)
#'
#' @export
#'
calculate_r0_rc_fromdata=function(df, f=1/72, gamma=1/223, r=1/60,
                              return.all=F, h.cutoff=5e-08){

  if(!"h" %in% names(df)){ stop("no h variable in df")}
  if(!"prop_import" %in% names(df)){ stop("no prop_import variable in df")}
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
  # convert incidence and initialisation
  dataTransform = df %>%
    dplyr::mutate( lambda=-1, Rc=-1,R0=-1, delta=NA, I=NA)

  pb <- utils::txtProgressBar(min = 0, max = nrow(dataTransform), style = 3)
  # calculate R0 and Rc with formula
  for(i in 1:nrow(dataTransform)){
    my.h=dataTransform$h[i]
    my.p=dataTransform$prop_import[i]
    my.rho=dataTransform$rho[i]
    my.alpha=dataTransform$alpha[i]
    my.beta=dataTransform$beta[i]
    my.omega=dataTransform$omega[i]
    if(is.na(my.h)==FALSE & is.infinite(my.h)==FALSE  &  (my.h)>=h.cutoff){
      my.lambda.solve=solve_lambda_vivax(h=my.h, r=r,  gamma=gamma, f=f,
                                                        alpha=my.alpha, beta=my.beta,
                                                        rho=my.rho,
                                                        p=my.p, omega=my.omega,
                                                        return.all = T)
      my.lambda=my.lambda.solve$lambda_real
      my.lambda=ifelse( (length(my.lambda)==0), -2, my.lambda)

    } else{my.lambda =-3}
    dataTransform$delta[i]=ifelse( my.lambda>0, my.lambda.solve$delta , NA)
    dataTransform$I[i]=ifelse( my.lambda>0, my.lambda.solve$I , NA)

    dataTransform$lambda[i]=my.lambda
    dataTransform$R0[i]=ifelse(my.lambda>0, get_r0_vivax(my.lambda,  f=f,r=r, gamma=gamma), NA)
    dataTransform$Rc[i]=ifelse(my.lambda>0, get_rc_vivax(my.lambda,  f=f,r=r, gamma=gamma, alpha=my.alpha, beta=my.beta, omega=my.omega), NA)
    utils::setTxtProgressBar(pb, i)
  }

  output_small=dataTransform[,! names(dataTransform) %in% c("delta", "I")]

  if(return.all){
    return(dataTransform)
  } else{
    return(output_small)
  }

}





#' @title Simulate from equilibrium
#'
#' @description Uses the compartmental model to simulate a trajectory, starting from equilibrium variable
#'
#' @param df a dataframe containing the data, with one column called I containing the proportion
#' of infectious individuals (I0+Il) at equilibrium, one column called lambda containing the transmission rate,
#' and one variable called id which identifies uniquely each row in the dataset.
#' Additional optional variables are: \cr
#' rho (reporting rate), delta (importation rate) \cr
#' intervention levels in the past (when lambda was calculated):  alpha.old (effective care), beta.old (proportion of liver stage cure), omega.old (intensity of vector control)
#' intervention levels in the future (in the simulation):  alpha.new (effective care), beta.new (proportion of liver stage cure), omega.new (intensity of vector control)
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param maxtime number of time steps for simulation
#' @param year if TRUE, aggregates the outputs per year (h would be in cases per person year). if FALSE, returns daily outputs (h would be in cases per person day).
#'
#' @return A dataframe with the simulated state variables for each parameter combination in df
#'
#' @details If alpha is not provided in df, alpha=0. If beta is not provided in df, beta=1.
#' If rho is not provided in df, rho=1. If omega is not provided in df, omega=1.
#' If delta is not provided in df, delta=0.
#'
#' @examples
#' mydata=data.frame(incidence=c(23,112),lambda=c(0.0063,0.0071),I=c(0.017,0.12),id=c(1,2))
#' mydata$rho=c(0.18,0.13)
#' mydata$beta.old=c(0.43,0.42)
#' mydata$alpha.old=c(0.17, 0.12)
#' mydata$delta=c(0,0)
#' mydata$omega.old=c(1,1)
#' simulate_from_equilibrium_fromdata(df=mydata, f=1/69, gamma=1/383, r=1/60,maxtime=2000,year=TRUE)
#'
#' @export
#'
simulate_from_equilibrium_fromdata=function(df, f=1/72, gamma=1/223, r=1/60,
                                        maxtime,year){

  if(!"I" %in% names(df)){ stop("no I variable in df")}
  if(!"id" %in% names(df)){ stop("no id variable in df")}
  if(!"lambda" %in% names(df)){ stop("no lambda variable in df")}
  if(!"rho" %in% names(df)){
    df$rho=1
    warning("no rho in df, assumed rho=1")
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
  if(!"delta" %in% names(df)){
    df$delta=0
    warning("no delta in df, assumed delta=0")
  }

  db_all=data.frame()
  for(i in 1:nrow(df)){
    # get equilibrium values for Sl, S0, Il and I0
    equ_states=get_equilibrium_states_vivax(I=df[i,]$I, lambda=as.numeric(df[i,]$lambda), r=r, gamma=gamma, f=f,
                                      alpha=df$alpha.old[i], beta=df$beta.old[i], rho=df$rho[i], delta=df$delta[i], omega=df$omega.old[i])

    # simul approx model delay
    this.simul=simulate_vivax_ode(parameters=list("r"=r,"gamma"=gamma, "f"=f, "lambda"=as.numeric(df[i,]$lambda), "delta"=df$delta[i],
                                                  "alpha"=df$alpha.new[i],
                                                 "beta"=df$beta.new[i],"rho"=df$rho[i],"omega"=df$omega.new[i],
                                                 "I0"=equ_states$I0, "S0"=equ_states$S0, "Sl"=equ_states$Sl, "Il"=equ_states$Il,
                                                 "h"=equ_states$h, "hr"=equ_states$hr),
                                 ODEmodel =ode_vivax_cm_importation_vc , maxtime = maxtime, year=year)
    if(year) {
      this.simul$incidence= this.simul$h*1000
      } else { this.simul$incidence = incidence_day2year(this.simul$h)}
    this.simul$id=df[i,]$id
    # combine the files
    db_all=rbind(db_all, this.simul)

  } #end for

  return(db_all)

}



#' @title Formating intervention parameters within a database
#'
#' @description Helper function to format dataframes
#'
#' @param df a dataframe containing at least 3 variables : alpha, beta and omega
#' @param intervention_object an named list containing the intervention description. It should have the follwing structure:
#' list(intervention_name="string", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA )
#' where NA can be replaced by scalars or kepts as such.
#'
#' @return A dataframe with the required input parameters for future simulation
#'
#' @details If alpha.new is not provided in intervention_object it is equal to alpha.
#' If beta.new is not provided in intervention_object it is equal to beta.
#' If omega.new is not provided in intervention_object it is equal to omega.
#'
format_data_simulation=function(df, intervention_object){

  if(!"alpha" %in% names(df)){
    df$alpha=0
    warning("no alpha in df, assumed alpha=0")
  }
  if(!"beta" %in% names(df)){
    df$beta=1
    warning("no beta in df, assumed beta=1")
  }
  if(!"omega" %in% names(df)){
    df$omega=1
    warning("no omega in df, assumed omega=1")
  }

  new_df=df
  new_df$alpha.old=new_df$alpha
  new_df$beta.old=new_df$beta
  new_df$omega.old=new_df$omega

  if(is.na(intervention_object$alpha.new)) new_df$alpha.new = new_df$alpha else new_df$alpha.new = intervention_object$alpha.new
  if(is.na(intervention_object$beta.new)) new_df$beta.new = new_df$beta else new_df$beta.new = intervention_object$beta.new
  if(is.na(intervention_object$omega.new)) new_df$omega.new = new_df$omega else new_df$omega.new = intervention_object$omega.new

  new_df$intervention=intervention_object$intervention_name
  return(new_df)
}


#' @title Simulate the effect of several interventions
#'
#' @description Uses the compartmental model to simulate the effect of interventions
#'
#' @param df a dataframe containing the data, with one column called I containing the proportion
#' of infectious individuals (I0+Il) at equilibrium, one column called lambda containing the transmission rate,
#' and one variable called id which identifies uniquely each row in the dataset.
#' Additional optional variables are: \cr
#' rho (reporting rate), delta (importation rate) \cr
#' intervention levels in the past (when lambda was calculated):  alpha.old (effective care), beta.old (proportion of liver stage cure), omega.old (intensity of vector control)
#' intervention levels in the future (in the simulation):  alpha.new (effective care), beta.new (proportion of liver stage cure), omega.new (intensity of vector control)
#' In practice, this dataframe can be the output of the function calculate_r0_rc_fromdata
#' @param intervention_list a list of intervention objects.
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param maxtime number of time steps for simulation
#' @param year if TRUE, aggregates the outputs per year (h would be in cases per person year). if FALSE, returns daily outputs (h would be in cases per person day).
#'
#' @return A dataframe with the simulated state variables for each parameter combination in df
#'
#' @examples
#' mydata=data.frame(incidence=c(23,112),lambda=c(0.0063,0.0071),I=c(0.017,0.12),id=c(1,2))
#' mydata$rho=c(0.18,0.13)
#' mydata$beta.old=c(0.43,0.42)
#' mydata$alpha.old=c(0.17, 0.12)
#' mydata$delta=c(0,0)
#' mydata$omega.old=c(1,1)
#' int_0=list(intervention_name="0", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA )
#' int_A=list(intervention_name="A", "alpha.new"=NA, "beta.new"=0.6, "omega.new"=NA )
#' my_intervention_list=list(int_0, int_A)
#' simulate_vivax_interventions(df=mydata, my_intervention_list)
#'
#' @export
#'
#' @details
#' An intervention object is named list containing the intervention description. It should have the follwing structure:
#' list(intervention_name="string", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA )
#' where NA can be replaced by scalars or kepts as such.
#' If alpha.new is not provided in intervention_object it is equal to alpha.
#' If beta.new is not provided in intervention_object it is equal to beta.
#' If omega.new is not provided in intervention_object it is equal to omega.
#'
simulate_vivax_interventions=function(df, intervention_list, f=1/72, gamma=1/223, r=1/60, year=T,maxtime=2000){

  df_full=data.frame()
  for (interv in intervention_list){
    df_full=rbind(df_full, format_data_simulation(df, interv ))
  }
  df_full$id0=df_full$id
  df_full$id=paste0(df_full$id0, df_full$intervention)
  merging_table=df_full[c("id", "id0", "intervention")]

  simulation_model= simulate_from_equilibrium_fromdata(df=df_full, f=f, gamma=gamma, r=r,maxtime=maxtime,year=year)
  simulation_model=merge(simulation_model,merging_table )
  simulation_model$id=simulation_model$id0
  simulation_model$id0=NULL

  return(simulation_model)

}
