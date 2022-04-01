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

  df=sanity_checks_inputs_calculate(df=df, delay=FALSE, rcd_at_baseline=FALSE)

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






#' @title Simulate the model for several values in a database
#'
#' @description Uses the compartmental model to simulate a trajectory, starting from equilibrium equilibrium state or a given initial condition
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
#' "Il_init", "I0_init", "Sl_init", "S0_init", "h_init", "hr_init",
#' with one variable called id which identifies uniquely each row in the dataset. This input is not used when from_equilibrium=TRUE (default).
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param maxtime number of time steps for simulation
#' @param year if TRUE, aggregates the outputs per year (h would be in cases per person year). if FALSE, returns daily outputs (h would be in cases per person day).
#' @param rcd if TRUE, includes RCD in the model. Default is FALSE (no RCD)
#' @param mda a boolean indicating if the model including mass drug administration (MDA) prophylaxis should be used. Default (FALSE) is the model without MDA
#' @param rcd_at_baseline a boolean indicating if the model was calibrated using the RCD model (i.e. there is some RCD at baseline already). Default (FALSE) is the model without RCD at baseline
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
#' simulate_from_data(df=mydata, f=1/69, gamma=1/383, r=1/60,maxtime=2000,year=TRUE)
#'
#' @export
#'
simulate_from_data=function(df, from_equilibrium=TRUE, initial_states=NULL, f=1/72, gamma=1/223, r=1/60,  maxtime, year, rcd=FALSE, mda=FALSE, rcd_at_baseline=FALSE){

  ####################################
  # sanity checks on the inputs
  df=sanity_checks_inputs_simulate(df=df, delay=FALSE, rcd=rcd, mda=mda, sto=FALSE, rcd_at_baseline=rcd_at_baseline)

  # further sanity checks, depending if equilibrium or non equilibrium is used
  if(from_equilibrium==FALSE){
    if(is.null(initial_states)){ stop("initial_states is missing")}
    if(!"id" %in% names(initial_states)){ stop("no id variable in initial_states")}
    if((!"Il_init" %in% names(initial_states) | !"I0_init" %in% names(initial_states) | !"Sl_init" %in% names(initial_states) |
        !"S0_init" %in% names(initial_states)|
        !"h_init" %in% names(initial_states)| !"hr_init" %in% names(initial_states)|
        !"hl_init" %in% names(initial_states)| !"hh_init" %in% names(initial_states)| !"hhl_init" %in% names(initial_states))){
      stop("missing state variables or wrong names in initial_states")}
    if(nrow(df) != nrow(initial_states) ){ stop("mydata and myinitialstate don't have the same number of rows")}
    if( any((initial_states$Il_init+initial_states$I0_init+initial_states$Sl_init+initial_states$S0_init)-1>1e-10)){ stop("initial states do not sum to 1")}

    ######################
    # processing the initial condition, if from_equilibrium=FALSE
    df=merge(df, initial_states)

  } else{
    message("simulating from equilibrium")
    if(!"I" %in% names(df)){ stop("no I variable in df")}
    if(rcd){message("We start from the equilibrium without RCD")}
  }


  #####################################
  # indicate the model to simulate from (rcd or not, mda or not)
  my_ode_model=ifelse(rcd, ode_vivax_rcd, ode_vivax_cm_importation_vc)

  if(mda){
    my_ode_model_mda=ifelse(rcd, ode_vivax_rcd_mda, ode_vivax_mda)
  }

  #####################################
  # Large loop over all areas/scenarios
  db_all=data.frame()
  for(i in 1:nrow(df)){

    # create parameter object
    myparameters=list("r"=r,"gamma"=gamma, "f"=f, "lambda"=as.numeric(df[i,]$lambda), "delta"=df$delta.new[i][[1]],
                      "alpha"=df$alpha.new[i],
                      "beta"=df$beta.new[i],"rho"=df$rho.new[i],"omega"=df$omega.new[i][[1]])

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

    if(from_equilibrium==TRUE){
      # get equilibrium values for Sl, S0, Il and I0
      equ_states=get_equilibrium_states_vivax(I=df[i,]$I, lambda=as.numeric(df[i,]$lambda), r=r, gamma=gamma, f=f,
                                              alpha=df$alpha.old[i], beta=df$beta.old[i], rho=df$rho.old[i],
                                              delta=ifelse(is.numeric(df$delta[i]), df$delta[i], df$delta[i](0)),
                                              omega=ifelse(is.numeric(df$omega.old[i]), df$omega.old[i], df$omega.old[i](0)))

      if(rcd_at_baseline){
        equ_states=get_equilibrium_states_vivax_rcd(I=df[i,]$I, lambda=as.numeric(df[i,]$lambda), r=r, gamma=gamma, f=f,
                                                alpha=df$alpha.old[i], beta=df$beta.old[i], rho=df$rho.old[i],
                                                delta=ifelse(is.numeric(df$delta[i]), df$delta[i], df$delta[i](0)),
                                                omega=ifelse(is.numeric(df$omega.old[i]), df$omega.old[i], df$omega.old[i](0)),
                                                iota_star =df$iota_star[i], nu=df$nu.old[i], tau=df$tau.old[i],eta=df$eta.old[i],kappa=df$kappa.old[i])
      }

      myparameters$I0=equ_states$I0
      myparameters$S0=equ_states$S0
      myparameters$Il=equ_states$Il
      myparameters$Sl=equ_states$Sl
      myparameters$h=equ_states$h
      myparameters$hr=equ_states$hr
      myparameters$hl=equ_states$hl
      myparameters$hh=equ_states$hh
      myparameters$hhl=equ_states$hhl
    } else{
      # extract from file
      myparameters$I0=df$I0_init[i]
      myparameters$S0=df$S0_init[i]
      myparameters$Il=df$Il_init[i]
      myparameters$Sl=df$Sl_init[i]
      myparameters$h=df$h_init[i]
      myparameters$hr=df$hr_init[i]
      myparameters$hl=df$hl_init[i]
      myparameters$hh=df$hh_init[i]
      myparameters$hhl=df$hhl_init[i]
    }

    # simulate the model

    if(mda){
      this.simul=simulate_vivax_mda_ode(parameters=myparameters, ODEmodel =my_ode_model, ODEmodel_mda = my_ode_model_mda , maxtime = maxtime, year=year)
    } else{
      this.simul=simulate_vivax_ode(parameters=myparameters, ODEmodel =my_ode_model , maxtime = maxtime, year=year)
    }

    # post-process
    if(year) {
      this.simul$incidence= this.simul$h*1000
    } else { this.simul$incidence = incidence_day2year(this.simul$h)}
    this.simul$id=df[i,]$id

    # combine the files
    db_all=rbind(db_all, this.simul)

  } #end for

  return(db_all)

}

