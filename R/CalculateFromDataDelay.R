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
#' mydata$sigma=c(1/15, 1/15)
#' mydata$omega=c(1,1)
#' calculate_r0_rc_fromdata_delay(df=mydata, f=1/72, gamma=1/223, r=1/60, return.all = TRUE)
#'
#' @export
#'
calculate_r0_rc_fromdata_delay=function(df, f=1/72, gamma=1/223, r=1/60,
                                  return.all=F, h.cutoff=5e-08){

  if(!"h" %in% names(df)){ stop("no h variable in df")}
  if(!"prop_import" %in% names(df)){ stop("no prop_import variable in df")}
  if(!"sigma" %in% names(df)){ stop("no sigma variable in df, maybe use the model without delay (calculate_r0_rc_fromdata)")}
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
    my.sigma=dataTransform$sigma[i]
    my.omega=dataTransform$omega[i]
    if(is.na(my.h)==FALSE & is.infinite(my.h)==FALSE  &  (my.h)>=h.cutoff){
      my.lambda.solve=solve_lambda_vivax_delay(h=my.h, r=r,  gamma=gamma, f=f,
                                         alpha=my.alpha, beta=my.beta, sigma=my.sigma,
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
    dataTransform$Rc[i]=ifelse(my.lambda>0, get_rc_vivax_delay(my.lambda,  f=f,r=r, gamma=gamma, alpha=my.alpha, beta=my.beta, sigma=my.sigma, omega=my.omega), NA)
    utils::setTxtProgressBar(pb, i)
  }

  output_small=dataTransform[,! names(dataTransform) %in% c("delta", "I")]

  if(return.all){
    return(dataTransform)
  } else{
    return(output_small)
  }

}


#' @title Simulate from equilibrium, with delays in treatment
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
#' @param from_equilibrium boolean indicating if the model is run from equilibrium (TRUE, default) or from a pre-specified initial condition, which should be specified in initial_states
#' @param initial_states given initial condition, as a dataframe containing the variables
#' "Il_init", "I0_init", "Sl_init", "S0_init", "Tl_init", "T0_init", "h_init",
#' with one variable called id which identifies uniquely each row in the dataset. This input is not used when from_equilibrium=TRUE (default).
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param maxtime number of time steps for simulation
#' @param year if TRUE, aggregates the outputs per year (h would be in cases per person year). if FALSE, returns daily outputs (h would be in cases per person day).
#' @param rcd a boolean indicating if the model including reactive case detection should be used. Default (FALSE) is the model without RCD
#' @param referral a boolean indicating if the rcd model includes referral. Default (FALSE) is the model without referral for RCD. This parameter is used only if rcd==TRUE.
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
#' mydata$sigma.old=c(1/15,1/15)
#' simulate_from_data_delay(df=mydata, f=1/69, gamma=1/383, r=1/60,maxtime=2000,year=TRUE)
#'
#' @export
#'
simulate_from_data_delay=function(df, from_equilibrium=TRUE, initial_states=NULL, f=1/72, gamma=1/223, r=1/60,   maxtime, year, rcd=FALSE, referral=FALSE, mda=FALSE, rcd_at_baseline=FALSE, sto=FALSE,
                                  sto_method="exact", runs=1, seeds=NULL){

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
  if(!"sigma.old" %in% names(df)){
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
  if(!"sigma.new" %in% names(df)){
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

  if((!"iota.old" %in% names(df) | !"tau.old" %in% names(df) | !"nu.old" %in% names(df) | !"eta.old" %in% names(df)| !"kappa.old" %in% names(df)) & rcd_at_baseline==TRUE){
    stop("RCD parameters are missing at baseline, please add them or use model without RCD (rcd_at_baseline=F)")
  }

  if((!"MDAcov.new" %in% names(df) | !"MDAp_length.new" %in% names(df)) & mda==TRUE){
    stop("MDA parameters are missing, please add them or use model without MDA")
  }

  if((!"N" %in% names(df) ) & sto==TRUE){
    stop("population size N is missing for stochastic model. Please add it or use deterministic model")
  }

  # further sanity checks, depending if equilibrium or non equilibrium is used
  if(from_equilibrium==FALSE){
    if(is.null(initial_states)){ stop("initial_states is missing")}
    if(!"id" %in% names(initial_states)){ stop("no id variable in initial_states")}
    if((!"Il_init" %in% names(initial_states) | !"I0_init" %in% names(initial_states) | !"Sl_init" %in% names(initial_states) |
        !"S0_init" %in% names(initial_states)|!"Tl_init" %in% names(initial_states)|!"T0_init" %in% names(initial_states)|
        !"h_init" %in% names(initial_states) |!"hl_init" %in% names(initial_states)|
       !"hh_init" %in% names(initial_states) |!"hhl_init" %in% names(initial_states))) {
      stop("missing state variables or wrong names in initial_states")}
    if(nrow(df) != nrow(initial_states) ){ stop("mydata and myinitialstate don't have the same number of rows")}

    ######################
    # processing the initial condition, if from_equilibrium=FALSE
    df=merge(df, initial_states)

    if(sto){
      df[,c("Il_init", "I0_init","Sl_init", "S0_init","Tl_init", "T0_init")]=df[,c("Il_init", "I0_init","Sl_init", "S0_init","Tl_init", "T0_init")]/df$N
    }
    if( any(abs((df$Il_init+df$I0_init+df$Sl_init+df$S0_init+df$Tl_init+df$T0_init)-1)>1e-10)){ stop("initial states do not sum to 1")}


  } else{
    message("simulating from equilibrium")
    if(!"I" %in% names(df)){ stop("no I variable in df")}
    if(rcd){ message("We start from the equilibrium without RCD")}
  }
  #####################################
  # indicate the model to simulate from (rcd or not, mda or not)

  if(sto==TRUE & rcd==TRUE & any(is.numeric(df$tau))==FALSE){stop("stochastic version of RCD model with time varying tau is not implemented yet")}

  my_ode_model=ifelse(rcd, ifelse(referral,ode_vivax_delay_rcd_referral,ode_vivax_delay_rcd_no_referral),
                     ode_vivax_delay)

  if(sto){
    if(rcd){
      if(referral){
        my_sto_model=model_sto_vivax_delay_rcd_referral()
      } else{
        my_sto_model=model_sto_vivax_delay_rcd_no_referral()
      }
    } else {
      my_sto_model=  model_sto_vivax_delay()
    }
  }

  if(mda){
    my_ode_model_mda=ifelse(rcd, ifelse(referral,ode_vivax_delay_rcd_referral_mda, ode_vivax_delay_rcd_no_referral_mda),
                            ode_vivax_delay_mda)

    if(sto){
      if(rcd){
        if(referral){
          my_sto_model_mda=model_sto_vivax_delay_rcd_referral_mda()
        } else{
          my_sto_model_mda=model_sto_vivax_delay_rcd_no_referral_mda()
        }
      } else {
        my_sto_model_mda=  model_sto_vivax_delay_mda()
      }
    }
  }

  #####################################
  # Large loop over all areas/scenarios
  db_all=data.frame()
  for(i in 1:nrow(df)){

    # create parameter object
    myparameters=list("r"=r,"gamma"=gamma, "f"=f, "lambda"=as.numeric(df[i,]$lambda), "delta"=df$delta.new[i][[1]],
                      "alpha"=df$alpha.new[i],
                      "beta"=df$beta.new[i],"sigma"=df$sigma.new[i],"rho"=df$rho.new[i],"omega"=df$omega.new[i][[1]])

    if(rcd){
      myparameters$iota=df$iota.new[i]
      myparameters$tau=ifelse(is.list(df$tau.new[i]),df$tau.new[i][[1]],df$tau.new[i])
      myparameters$nu=df$nu.new[i]
      myparameters$eta=df$eta.new[i]
      myparameters$kappa=df$kappa.new[i]
    }

    if(mda){
      myparameters$MDAcov=df$MDAcov.new[i]
      myparameters$MDAp_length=df$MDAp_length.new[i]
      myparameters$MDArad_cure=df$MDArad_cure.new[i]
    }

    if(sto){
      myparameters$N=df$N[i]
    }

    if(from_equilibrium==TRUE){
      # get equilibrium values for Sl, S0, Il and I0
      equ_states=get_equilibrium_states_vivax_delay(I=df[i,]$I, lambda=as.numeric(df[i,]$lambda), r=r, gamma=gamma, f=f,
                                                    alpha=df$alpha.old[i], beta=df$beta.old[i],sigma=df$sigma.old[i], rho=df$rho.old[i],
                                                    delta=ifelse(is.numeric(df$delta[i]), df$delta[i], df$delta[i](0)),
                                                    omega=ifelse(is.numeric(df$omega.old[i]), df$omega.old[i], df$omega.old[i](0)))

      if(rcd_at_baseline){
        if(referral){
          equ_states=get_equilibrium_states_vivax_rcd_referral(I=df[i,]$I, lambda=as.numeric(df[i,]$lambda), r=r, gamma=gamma, f=f,
                                                                  alpha=df$alpha.old[i], beta=df$beta.old[i],sigma=df$sigma.old[i], rho=df$rho.old[i],
                                                                  delta=ifelse(is.numeric(df$delta[i]), df$delta[i], df$delta[i](0)),
                                                                  omega=ifelse(is.numeric(df$omega.old[i]), df$omega.old[i], df$omega.old[i](0)),
                                                                  iota_star =df$iota_star[i], nu=df$nu.old[i], tau=df$tau.old[i],eta=df$eta.old[i],kappa=df$kappa.old[i])
        }else{
          equ_states=get_equilibrium_states_vivax_rcd_no_referral(I=df[i,]$I, lambda=as.numeric(df[i,]$lambda), r=r, gamma=gamma, f=f,
                                                                  alpha=df$alpha.old[i], beta=df$beta.old[i],sigma=df$sigma.old[i], rho=df$rho.old[i],
                                                                  delta=ifelse(is.numeric(df$delta[i]), df$delta[i], df$delta[i](0)),
                                                                  omega=ifelse(is.numeric(df$omega.old[i]), df$omega.old[i], df$omega.old[i](0)),
                                                                  iota_star =df$iota_star[i], nu=df$nu.old[i], tau=df$tau.old[i],eta=df$eta.old[i],kappa=df$kappa.old[i])
        }
      }

      myparameters$I0=equ_states$I0
      myparameters$S0=equ_states$S0
      myparameters$Il=equ_states$Il
      myparameters$Sl=equ_states$Sl
      myparameters$Tl=equ_states$Tl
      myparameters$T0=equ_states$T0
      myparameters$h=equ_states$h
      myparameters$hl=equ_states$hl
      myparameters$hh=equ_states$hh
      myparameters$hhl=equ_states$hhl
    } else{
      # extract from file
      myparameters$I0=df$I0_init[i]
      myparameters$S0=df$S0_init[i]
      myparameters$Il=df$Il_init[i]
      myparameters$Sl=df$Sl_init[i]
      myparameters$Tl=df$Tl_init[i]
      myparameters$T0=df$T0_init[i]
      myparameters$h=df$h_init[i]
      myparameters$hl=df$hl_init[i]
      myparameters$hh=df$hh_init[i]
      myparameters$hhl=df$hhl_init[i]

      seeds=df$run[i]
    }

    # simul approx model delay

    if(sto){

      if(mda){
        this.simul=simulate_vivax_delay_mda_sto(parameters=myparameters, STOmodel=my_sto_model, STOmodel_mda = my_sto_model_mda , maxtime = maxtime, year=year, sto_method=sto_method, runs=runs, seeds=seeds)
      } else{
        this.simul=simulate_vivax_delay_sto(parameters=myparameters, STOmodel=my_sto_model, maxtime=maxtime, year=year, sto_method=sto_method, runs=runs, seeds=seeds)
      }
      if("run" %in% names(df)){this.simul$run=df$run[i]}
    } else {
      if(mda){
        this.simul=simulate_vivax_delay_mda_ode(parameters=myparameters, ODEmodel =my_ode_model, ODEmodel_mda = my_ode_model_mda , maxtime = maxtime, year=year)
      } else{
        this.simul=simulate_vivax_delay_ode(parameters=myparameters,   ODEmodel =my_ode_model , maxtime = maxtime, year=year)
      }
    }


    # post-process
    if(year) {
      this.simul$incidence= this.simul$h*1000
    } else { this.simul$incidence = incidence_day2year(this.simul$h)}

    if(sto){
      this.simul$incidence=this.simul$incidence/myparameters$N
    }

    this.simul$id=df[i,]$id

    # combine the files
    db_all=rbind(db_all, this.simul)

  } #end for

  return(db_all)

}
