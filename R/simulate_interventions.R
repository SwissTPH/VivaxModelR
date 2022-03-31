#' @title Formating intervention parameters within a database
#'
#' @description Helper function to format dataframes
#'
#' @param df a dataframe containing at least 3 variables : alpha, beta and omega
#' @param intervention_object an named list containing the intervention description. It should have the follwing structure:
#' list(intervention_name="string", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA )
#' where NA can be replaced by scalars or kepts as such.
#' @param delay a boolean indicating if the model including delays in treatment should be used. Default (FALSE) is the model without delay in treatment
#' @param rcd a boolean indicating if the model including reactive case detection should be used. Default (FALSE) is the model without RCD
#' @param mda a boolean indicating if the model including mass drug administration (MDA) prophylaxis should be used. Default (FALSE) is the model without MDA
#' @param mystart an integer indicating the time index of the start of the simulation (this is used to chain interventions appropriately)
#' @param rcd_at_baseline a boolean indicating if the model was calibrated using the RCD model (i.e. there is some RCD at baseline already). Default (FALSE) is the model without RCD at baseline
#'
#' @return A dataframe with the required input parameters for future simulation
#'
#' @details If alpha.new is not provided in intervention_object it is equal to alpha.
#' If beta.new is not provided in intervention_object it is equal to beta.
#' If omega.new is not provided in intervention_object it is equal to omega.
#' If rho.new is not provided in intervention_object it is equal to rho
#'
format_data_simulation=function(df, intervention_object, delay=FALSE, rcd=FALSE, mda=FALSE, mystart, rcd_at_baseline=FALSE){

  if(!"rho" %in% names(df)){
    df$rho=1
    warning("no rho in df, assumed rho=0")
  }
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


  if(rcd==FALSE & rcd_at_baseline==TRUE){
    stop("there is RCD at baseline but no RCD in the future: this is not supported by the model")
  }

  new_df=df
  new_df$rho.old=new_df$rho
  new_df$alpha.old=new_df$alpha
  new_df$beta.old=new_df$beta
  new_df$omega.old=new_df$omega

  if(is.na(intervention_object$rho.new)) new_df$rho.new = new_df$rho else new_df$rho.new = intervention_object$rho.new
  if(is.na(intervention_object$alpha.new)) new_df$alpha.new = new_df$alpha else new_df$alpha.new = intervention_object$alpha.new
  if(is.na(intervention_object$beta.new)) new_df$beta.new = new_df$beta else new_df$beta.new = intervention_object$beta.new


  if(is.data.frame(intervention_object$omega.new)){
    if((!"time" %in% names(intervention_object$omega.new) ) | (!"value" %in% names(intervention_object$omega.new) )){
      stop("omega is a dataframe but does not have the correct variable names, i.e. time and value")
    } else{
      # create the linear interpolation, and ensure that the time start matches if interventions are chained
      omega_t=create_approximation(x=intervention_object$omega.new$time, y=intervention_object$omega.new$value,
                                   start=mystart)
      new_df$omega.new = rep(list(omega_t),nrow(new_df))
    }

  } else {
    if(is.na(intervention_object$omega.new)){
      new_df$omega.new =new_df$omega
    } else{
      new_df$omega.new =  intervention_object$omega.new
    }
  }


  if(!"delta.new" %in% names(intervention_object)) {
    new_df$delta.new = new_df$delta
  } else if(is.numeric(intervention_object$delta.new)){
    if(is.na(intervention_object$delta.new)){
      new_df$delta.new =new_df$delta
    } else{
      new_df$delta.new =  intervention_object$delta.new
    }
  } else { # accomodate for potential time varying delta
    if(is.data.frame(intervention_object$delta.new)){
      if((!"time" %in% names(intervention_object$delta) ) | (!"value" %in% names(intervention_object$delta) )){
        stop("delta is a dataframe but does not have the correct variable names, i.e. time and value")
      } else{
        if(!"id" %in% names(intervention_object$delta) ){
          # create the linear interpolation, and ensure that the time start matches if interventions are chained
          delta_t=create_approximation(x=intervention_object$delta$time, y=intervention_object$delta$value,
                                       start=mystart)
          new_df$delta.new=rep(list(delta_t),nrow(new_df))
        } else{
          new_df$delta.new=NA
          if(length(unique(new_df$id)) != length(unique(intervention_object$delta.new$id))){
            stop("some id values are missing in the delta.new dataset")
          }
          for (my_id in unique(new_df$id)){
            delta_t=create_approximation(x=intervention_object$delta.new$time[intervention_object$delta.new$id==my_id],
                                         y=intervention_object$delta.new$value[intervention_object$delta.new$id==my_id],
                                         start=mystart)
            new_df$delta.new[new_df$id==my_id]=list(delta_t)
          }
        }

      }

    }
  }



  if(!"sigma" %in% names(df) & delay==TRUE){
    stop("no sigma in df, please use model without delay")
  }
  if(delay==TRUE) {
    new_df$sigma.old=new_df$sigma
    if(is.na(intervention_object$sigma.new)) new_df$sigma.new = new_df$sigma else new_df$sigma.new = intervention_object$sigma.new
  }

  if(rcd==TRUE){
    if(!"iota.new" %in% names(intervention_object) | !"tau.new" %in% names(intervention_object)|
        !"nu.new" %in% names(intervention_object)| !"eta.new" %in% names(intervention_object)| !"kappa.new" %in% names(intervention_object) ){
      stop("some rcd parameters (iota.new, nu.new, tau.new, eta.new or kappa.new) are missing in df, please use model without rcd")
    }

    if(is.na(intervention_object$iota.new)){
      new_df$iota.new = 0
      warning("no iota parameter, assumed iota=0")
    } else new_df$iota.new = intervention_object$iota.new
    if(is.na(intervention_object$nu.new)){
      new_df$nu.new = 0
      warning("no nu parameter, assumed nu=0")
    }  else new_df$nu.new = intervention_object$nu.new
    if(is.na(intervention_object$tau.new)){
      new_df$tau.new = 1
      warning("no tau parameter, assumed tau=1")
    } else new_df$tau.new = ifelse(is.numeric(intervention_object$tau.new), intervention_object$tau.new, list(intervention_object$tau.new))
    if(is.na(intervention_object$eta.new)){
      new_df$eta.new = 0
      warning("no eta parameter, assumed eta=0")
    } else new_df$eta.new = intervention_object$eta.new
    if(is.na(intervention_object$kappa.new)){
      new_df$kappa.new = 1
      warning("no kappa parameter, assumed kappa=1")
    } else new_df$kappa.new = intervention_object$kappa.new


    if(rcd_at_baseline==TRUE){

      if(!"iota_star" %in% names(intervention_object) | !"tau.old" %in% names(intervention_object)|
         !"nu.old" %in% names(intervention_object)| !"eta.old" %in% names(intervention_object)| !"kappa.old" %in% names(intervention_object) ){
        stop("some rcd parameters (iota_star, nu.old, tau.old, eta.old or kappa.old) are missing in df, please use model without rcd")
      }

      new_df$nu.old=new_df$nu
      new_df$tau.old=new_df$tau
      new_df$eta.old=new_df$eta
      new_df$kappa.old=new_df$kappa
    }

  }

  if(mda==TRUE){
    if(!"MDAcov.new" %in% names(intervention_object) | !"MDAp_length.new" %in% names(intervention_object)|
       !"MDArad_cure.new" %in% names(intervention_object)){
      stop("some mda parameters (MDAcov.new, MDArad_cure.new or MDAp_length.new) are missing in df, please use model without MDA")
    }

    if(is.na(intervention_object$MDAcov.new)){
      new_df$MDAcov.new = 0
      warning("no MDAcov.new parameter, assumed MDAcov=0")
    } else new_df$MDAcov.new = intervention_object$MDAcov.new
    if(is.na(intervention_object$MDAp_length.new)){
      new_df$MDAp_length.new = 1
      warning("no MDAp_length parameter, assumed MDAp_length=1")
    }  else new_df$MDAp_length.new = intervention_object$MDAp_length.new
    if(is.na(intervention_object$MDArad_cure.new)){
      new_df$MDArad_cure.new = 0
      warning("no MDArad_cure parameter, assumed MDArad_cure=1")
    } else new_df$MDArad_cure.new = intervention_object$MDArad_cure.new

  }

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
#' @param previous_simulation the result of a previous simulation (e.g. the outcome of simulate_vivax_interventions or chain_vivax_interventions).
#' If NULL (default), the model is simulated from equilibrium
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param maxtime number of time steps for simulation
#' @param year if TRUE, aggregates the outputs per year (h would be in cases per person year). if FALSE, returns daily outputs (h would be in cases per person day).
#' @param delay a boolean indicating if the model including delays in treatment should be used. Default (FALSE) is the model without delay in treatment
#' @param rcd a boolean indicating if the model including reactive case detection should be used. Default (FALSE) is the model without RCD
#' @param referral a boolean indicating if the rcd model includes referral. Default (FALSE) is the model with referral for RCD. This parameter is used only if rcd==TRUE.
#' @param mda a boolean indicating if the model including mass drug administration (MDA) prophylaxis should be used. Default (FALSE) is the model without MDA
#' @param rcd_at_baseline a boolean indicating if the model was calibrated using the RCD model (i.e. there is some RCD at baseline already). Default (FALSE) is the model without RCD at baseline
#' @param sto a boolean indicating if the stochastic model is used. Default (FALSE) is the deterministic (ODE) model
#' @param sto_method a scalarindicating which simulation method is used.
#' Default ("exact") is Gillespie's direct method. Other options are "approximate" (tau-leap) or "mixed".
#' cf. the documentation of the TiPS package for more information.
#' @param runs number of draws of the stochastic model
#' @param seeds a vector of the length of runs containing the seeds for each simulation (don't use "0" which has another use in TiPS)
#'
#' @return A dataframe with the simulated state variables for each parameter combination in df
#'
#' @examples
#' mydata=data.frame(incidence=c(23,112),lambda=c(0.0063,0.0071),I=c(0.017,0.12),id=c(1,2))
#' mydata$rho.old=c(0.18,0.13)
#' mydata$beta.old=c(0.43,0.42)
#' mydata$alpha.old=c(0.17, 0.12)
#' mydata$delta=c(0,0)
#' mydata$omega.old=c(1,1)
#' int_0=list(intervention_name="0", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "rho.new"=NA )
#' int_A=list(intervention_name="A", "alpha.new"=NA, "beta.new"=0.6, "omega.new"=NA, "rho.new"=NA )
#' my_intervention_list=list(int_0, int_A)
#' simulate_vivax_interventions(df=mydata, my_intervention_list)
#'
#' @export
#'
#' @details
#' An intervention object is named list containing the intervention description. It should have the follwing structure:
#' list(intervention_name="string", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "rho.new"=NA )
#' where NA can be replaced by scalars or kepts as such.
#' If alpha.new is not provided in intervention_object it is equal to alpha.
#' If beta.new is not provided in intervention_object it is equal to beta.
#' If omega.new is not provided in intervention_object it is equal to omega.
#' If rho.new is not provided in intervention_object it is equal to rho.
#'
simulate_vivax_interventions=function(df, intervention_list, previous_simulation=NULL, f=1/72, gamma=1/223, r=1/60, year=T,maxtime=2000, delay=FALSE, rcd=FALSE ,referral=FALSE, mda=FALSE, rcd_at_baseline=FALSE,
                                      sto=FALSE,sto_method="exact", runs=1, seeds=NULL){

  #############################
  # sanity check
  if(delay==FALSE & "sigma" %in% names(df)){
    warning("There is a sigma parameter in df: are you sure you don't want to use the model with delay in treatment?")
  }

  # sanity check
  if(delay==FALSE & sto==TRUE){
    stop("The stochastic version of the model without delays is not available")
  }
  #############################
  # prepare intervention list and make them part of the id column
  df_full=data.frame()
  for (interv in intervention_list){

    df_full=rbind(df_full, format_data_simulation(df, interv, delay=delay, rcd=rcd, mda=mda, rcd_at_baseline=rcd_at_baseline,
                                                  mystart=ifelse(is.null(previous_simulation), 0, max(previous_simulation$time)) ))
  }
  df_full$id0=df_full$id
  df_full$id=paste0(df_full$id0, df_full$intervention)
  merging_table=df_full[c("id", "id0", "intervention")]

  #############################
  # extract initial condition from previous simulation if needed
  if(!is.null(previous_simulation)){
    if(((!"Tl" %in% names(previous_simulation) )| (!"T0" %in% names(previous_simulation)))& delay==TRUE){
      stop("Tl or T0 are missing in previous_simulation, please add them or use the model without delay in treatment")
    }
    if((("Tl" %in% names(previous_simulation) )| ("T0" %in% names(previous_simulation)))& delay==FALSE){
      stop("Tl and/or T0 are present in previous_simulation, maybe you should use the model with delay in treatment")
    }
    # extract the last time step of the previous simulation
    chain_time=max(previous_simulation$time)
    myinitialstate=previous_simulation[previous_simulation$time==chain_time,]

    if(sto){myvars= c("id", "intervention", "step", "run")} else{myvars=c("id", "intervention", "step")}

    names(myinitialstate)[!names(myinitialstate) %in% myvars]=paste0(names(previous_simulation)[!names(myinitialstate) %in% myvars], "_init")
    myinitialstate$id0=myinitialstate$id_init
    myinitialstate$id0=myinitialstate$id
    myinitialstate$id=paste0(myinitialstate$id0, myinitialstate$intervention)

    if(any(!myinitialstate$id %in% df_full$id) | any(!df_full$id %in% myinitialstate$id)){
      stop("The intervention names are not the same as in the previous simulation")
    }

    this_from_equilibrium=FALSE

    if(sto){  #if stochastic, match the number of runs from previous simulation
      runs_list=rep(1:runs, times=nrow(df_full))
      df_full <- df_full[rep(seq_len(nrow(df_full)), each = runs), ]
      df_full$run=runs_list
      runs=1
    }

  } else{
    message("Starting from equilibrium")
    this_from_equilibrium=TRUE
    myinitialstate=NULL

    if(rcd){
      message("We start from the equilibrium without RCD")
    }
  }


  if(delay==TRUE){
    simulation_model= simulate_from_data_delay(df=df_full, from_equilibrium=this_from_equilibrium, initial_states=myinitialstate,
                                         f=f, gamma=gamma, r=r,maxtime=maxtime,year=year, rcd=rcd, referral = referral, mda=mda, rcd_at_baseline=rcd_at_baseline,
                                         sto=sto,sto_method=sto_method, runs=runs)
 } else {
    simulation_model= simulate_from_data(df=df_full, from_equilibrium=this_from_equilibrium, initial_states=myinitialstate,
                                         f=f, gamma=gamma, r=r,maxtime=maxtime,year=year, rcd=rcd, mda=mda, rcd_at_baseline=rcd_at_baseline)
  }

  simulation_model=merge(simulation_model,merging_table )
  simulation_model$id=simulation_model$id0
  simulation_model$id0=NULL

  if(!is.null(previous_simulation)){
    simulation_model$time=simulation_model$time+chain_time

    simulation_model$step=max(previous_simulation$step)+1
    simulation_model=rbind(previous_simulation, simulation_model[simulation_model$time !=chain_time,])
  } else{
    simulation_model$step=1
  }

  return(simulation_model)

}




#' @title Calibrate the vivax model
#'
#' @description Uses the compartmental model to simulate the effect of interventions
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
#' @param delay a boolean indicating if the model including delays in treatment should be used. Default (FALSE) is the model without delay in treatment
#' @param rcd a boolean indicating if the model including reactive case detection should be used. Default (FALSE) is the model without RCD
#' @param referral a boolean indicating if the rcd model includes referral. Default (FALSE) is the model with referral for RCD. This parameter is used only if rcd==TRUE.
#'
#' @return A dataframe with 3 additional columns: lambda is the transmission rate,
#' R0 is the basic reproduction number, Rc is the controlled reproduction number. \cr
#' If lambda = -2, it means that the parameter combination does not correspond to a positive lambda solution with the model. \cr
#' If lambda = -3, it means the incidence data contains missing values or infinite values or values below h.cutoff.  \cr
#' If return.all=T, delta and I are also included.
#'
#' @details If alpha is not provided in df, alpha=0. If beta is not provided in df, beta=1.
#' If rho is not provided in df, rho=1. If omega is not provided in df, omega=1.
#' @export
#'
#' @example
#' mydata=data.frame(incidence=c(23,112),prop_import=c(0,0.1))
#' mydata$h=incidence_year2day(mydata$incidence)
#' mydata$beta=c(0.43,0.42)
#' mydata$alpha=c(0.17, 0.12)
#' mydata$rho=c(0.17, 0.12)
#' mydata$omega=c(1,1)
#' calibrate_vivax_equilibrium(df=mydata, f=1/72, gamma=1/223, r=1/60, return.all = TRUE)
#'
calibrate_vivax_equilibrium=function(df, f=1/72, gamma=1/223, r=1/60, delay=FALSE, rcd=FALSE ,referral=FALSE,
                                     return.all=F, h.cutoff=5e-08){

  #############################
  # sanity check
  if(delay==FALSE & "sigma" %in% names(df)){
    warning("There is a sigma parameter in df: are you sure you don't want to use the model with delay in treatment?")
  }

  if((!"iota" %in% names(df) | !"tau" %in% names(df) | !"nu" %in% names(df) | !"eta" %in% names(df)| !"kappa" %in% names(df)) & rcd==TRUE ){
    stop("RCD parameters are missing, please add them or use model without RCD")
  }

  if(delay){
    if(rcd){
      df_calibrated=calculate_r0_rc_fromdata_delay_rcd(df, f=f, gamma=gamma, r=r,return.all=return.all, h.cutoff=h.cutoff, referral = referral)
    } else{
      df_calibrated=calculate_r0_rc_fromdata_delay(df, f=f, gamma=gamma, r=r,return.all=return.all, h.cutoff=h.cutoff)
    }
  }else {
    if(rcd){
      df_calibrated=calculate_r0_rc_fromdata_rcd(df, f=f, gamma=gamma, r=r,return.all=return.all, h.cutoff=h.cutoff)
    } else{
      df_calibrated=calculate_r0_rc_fromdata(df, f=f, gamma=gamma, r=r,return.all=return.all, h.cutoff=h.cutoff)
    }
  }

  return(df_calibrated)

}
