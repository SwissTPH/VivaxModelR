#' @title Reaction associated with the vivax model with delays, in stochastic version
#' @description This function creates the stochastic model for vivax dynamics with delays to treatment.
#' @return A list with 2 objects:
#' model: the stochastic model compiled (as in the TiPS package)
#' reactions_incidence: a list of vectors indicating which reactions correspond to the incidence variables hh and hhl
#' @export
model_sto_vivax_delay=function(){
  reactions <- c("S0 [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*S0/N] -> Il", # untreated infection
                 "Sl [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*Sl/N] -> Il", # untreated infection
                 "I0 [(omega*lambda*(Il+I0+Tl+T0))*I0/N] -> Il", # untreated reinfection
                 "S0 [alpha*(omega*lambda*(Il+I0+Tl+T0))*S0/N] -> Tl", # treated infection
                 "Sl [alpha*(omega*lambda*(Il+I0+Tl+T0))*Sl/N] -> Tl", # treated infection
                 "T0 [(omega*lambda*(Il+I0+Tl+T0))*T0/N] -> Tl", # treated reinfection

                 "S0 [(1-alpha)*(delta)*S0] -> Il", # untreated infection imported
                 "Sl [(1-alpha)*(delta)*Sl] -> Il", # untreated infection imported
                 "I0 [(delta)*I0] -> Il", # untreated reinfection imported
                 "S0 [alpha*(delta)*S0] -> Tl", # treated infection imported
                 "Sl [alpha*(delta)*Sl] -> Tl", # treated infection imported
                 "T0 [(delta)*T0] -> Tl", # treated reinfection imported

                 "Sl [(1-alpha)*f*Sl] -> Il", # untreated relapse
                 "Sl [alpha*f*Sl] -> Tl", # treated relapse
                 "Il [gamma*Il] -> I0",  # liver clearance
                 "Tl [gamma*Tl] -> T0",  # liver clearance
                 "Sl [gamma*Sl] -> S0",  # liver clearance
                 "Il [r*Il] -> Sl", # recovery
                 "Tl [r*Tl] -> Sl", # recovery
                 "I0 [r*I0] -> S0", # recovery
                 "T0 [r*T0] -> S0", # recovery
                 "Tl [(1-beta)*sigma*Tl] -> Sl", # treatment no rad cure
                 "Tl [beta*sigma*Tl] -> S0", # treatment rad cure
                 "T0 [sigma*T0 ] -> S0" # treatment
  )

  vivax_simu <- TiPS::build_simulator(reactions)

  safe_run <- function(f, ...) {
    out <- list()
    while(! length(out)) {out <- f(...)}
    out
  }

  safe_vivax_simu <- function(...) safe_run(vivax_simu, ...)

  reactions_hhl=reactions[c(1,2,4,5,13,14)]
  reactions_importation=reactions[c(7,8,10,11)]
  reactions_hh=c(reactions_hhl, reactions_importation)
  reactions_rcd=c()

  return(list(model=safe_vivax_simu, reactions_incidence=list(reactions_hh=reactions_hh, reactions_hhl=reactions_hhl, reactions_rcd=reactions_rcd)))
}



#' @title Simulate vivax model with delay to treatment, with stochastic model
#' @description Simulates a draw of the chosen vivax stochastic model
#' @param parameters model parameters
#' @param STOmodel  stochastic model to be simulated (reactions)
#' @param maxtime maximal time step
#' @param year if TRUE, aggregates outputs per year. if FALSE, returns daily output
#' @param sto_method a scalarindicating which simulation method is used.
#' Default ("exact") is Gillespie's direct method. Other options are "approximate" (tau-leap) or "mixed".
#' cf. the documentation of the TiPS package for more information.
#' @param runs number of draws of the stochastic model
#' @param seeds a vector of the length of runs containing the seeds for each simulation (don't use "0" which has another use in TiPS)
#'
#' @return A dataframe containing the simulated model
#' @export
#' @importFrom rlang .data
simulate_vivax_delay_sto=function(parameters, STOmodel=model_sto_vivax_delay(), maxtime=1465, year=FALSE, sto_method="exact", runs=1, seeds=NULL){

  #simulation
  state      = round(c(Il= parameters$Il,
                 I0= parameters$I0,
                 Sl= parameters$Sl,
                 S0 = parameters$S0,
                 Tl= parameters$Tl,
                 T0 = parameters$T0)*parameters$N)
  if(sum(state)!=parameters$N){
    state[which.max(state)]=state[which.max(state)] - (sum(state) - parameters$N)
  }

  times <- c(0, maxtime)

  # if omega is a function, extract values
  if(is.function(parameters$omega)){
    parameters$omega=parameters$omega(seq(0, maxtime-1))
    times <- seq(0, maxtime)
  }

  if(is.function(parameters$delta)){
    parameters$delta=parameters$delta(seq(0, maxtime-1))
    times <- seq(0, maxtime)
  }

  safe_vivax_simu=STOmodel$model

  if(is.null(seeds)){seeds=1:runs}
  if(length(seeds)!=runs)(stop("seeds length does not match runs"))


  solutionVivax=data.frame()
  for (rr in 1:runs){
    traj_dm <- safe_vivax_simu(
      paramValues = parameters,
      initialStates = state,
      times = times, method=sto_method, seed=seeds[rr])

    this_solutionVivax=traj_dm$traj
    this_solutionVivax$run=seeds[rr]
    solutionVivax=rbind(solutionVivax, this_solutionVivax)
  }

  solutionVivax$is_hh=solutionVivax$Reaction %in% STOmodel$reactions_incidence$reactions_hh
  solutionVivax$is_hhl=solutionVivax$Reaction %in% STOmodel$reactions_incidence$reactions_hhl
  solutionVivax$is_rcd=solutionVivax$Reaction %in% STOmodel$reactions_incidence$reactions_rcd


  if(year){
    solutionVivax$time=ceiling(solutionVivax$Time /365)*365
    all_time_steps=unique(ceiling(seq(0, maxtime)/365))*365
  } else {
    solutionVivax$time=ceiling(solutionVivax$Time)
    all_time_steps=seq(0, maxtime)
  }
  solutionVivax=solutionVivax[solutionVivax$time <= maxtime,]
  all_time_steps=all_time_steps[all_time_steps <= maxtime]

  solutionVivax=solutionVivax  %>% dplyr::group_by(.data$time,.data$run) %>%
    dplyr::summarise(Il=.data$Il[.data$Time==max(.data$Time)],I0=.data$I0[.data$Time==max(.data$Time)],
                     Sl=.data$Sl[.data$Time==max(.data$Time)],S0=.data$S0[.data$Time==max(.data$Time)],
                     Tl=.data$Tl[.data$Time==max(.data$Time)],T0=.data$T0[.data$Time==max(.data$Time)],
                     hh=sum(.data$is_hh), hhl=sum(.data$is_hhl), rcd_reac=sum(.data$is_rcd))

  corr_factor_rcd=ifelse(any(solutionVivax$rcd_reac), (1-parameters$alpha-parameters$rho+parameters$rho*parameters$alpha/parameters$kappa)/(1-parameters$alpha),0)

  solutionVivax$h=solutionVivax$hh*parameters$rho + solutionVivax$rcd_reac*corr_factor_rcd
  solutionVivax$hl=solutionVivax$hhl*parameters$rho + solutionVivax$rcd_reac*corr_factor_rcd
  solutionVivax$I=solutionVivax$Il+solutionVivax$I0+solutionVivax$Tl+solutionVivax$T0

  solutionVivax$hh[solutionVivax$time==0]=ifelse(year, parameters$hh*365, parameters$hh)*parameters$N
  solutionVivax$hhl[solutionVivax$time==0]=ifelse(year, parameters$hhl*365, parameters$hhl)*parameters$N
  solutionVivax$h[solutionVivax$time==0]=ifelse(year, parameters$h*365, parameters$h)*parameters$N
  solutionVivax$hl[solutionVivax$time==0]=ifelse(year, parameters$hl*365, parameters$hl)*parameters$N

  solutionVivax$rcd_reac=NULL
  # add the lines for the time steps without reactions
  all_times=list()
  all_times$time=all_time_steps
  all_times$run=as.numeric(seeds)
  all_times_full = base::expand.grid( all_times )

  solutionVivax_full=base::merge(solutionVivax, all_times_full, all=T)

  solutionVivax_full=solutionVivax_full[  with(solutionVivax_full, order(run, time)),]
  index_list=base::which(is.na(solutionVivax_full$Il))
  for(i in index_list){
    solutionVivax_full[i,c("Il", "I0", "Sl", "S0", "Tl", "T0", "I")]=solutionVivax_full[(i-1),c("Il", "I0", "Sl", "S0", "Tl", "T0", "I")]
    solutionVivax_full[i,c("hh", "hhl", "h", "hl")]=c(0,0,0,0)
  }


  # calculate proportion of imported cases
  solutionVivax_full$p=(solutionVivax_full$h-solutionVivax_full$hl)/solutionVivax_full$h

  return(as.data.frame(solutionVivax_full))
}


#' @title Reaction associated with the vivax model with delays and MDA during prophylaxis, in stochastic version
#' @description This function creates the stochastic model for vivax dynamics with delays to treatment and MDA.
#' @return A list with 2 objects:
#' model: the stochastic model compiled (as in the TiPS package)
#' reactions_incidence: a list of vectors indicating which reactions correspond to the incidence variables hh and hhl
#' @export
model_sto_vivax_delay_mda=function(){
  reactions <- c("SS0 [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*SS0/N] -> Il", # untreated infection
                 "SSl [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*SSl/N] -> Il", # untreated infection
                 "I0 [(omega*lambda*(Il+I0+Tl+T0))*I0/N] -> Il", # untreated reinfection
                 "SS0 [alpha*(omega*lambda*(Il+I0+Tl+T0))*SS0/N] -> Tl", # treated infection
                 "SSl [alpha*(omega*lambda*(Il+I0+Tl+T0))*SSl/N] -> Tl", # treated infection
                 "T0 [(omega*lambda*(Il+I0+Tl+T0))*T0/N] -> Tl", # treated reinfection

                 "SS0 [(1-alpha)*(delta)*SS0] -> Il", # untreated infection imported
                 "SSl [(1-alpha)*(delta)*SSl] -> Il", # untreated infection imported
                 "I0 [(delta)*I0] -> Il", # untreated reinfection imported
                 "SS0 [alpha*(delta)*SS0] -> Tl", # treated infection imported
                 "SSl [alpha*(delta)*SSl] -> Tl", # treated infection imported
                 "T0 [(delta)*T0] -> Tl", # treated reinfection imported

                 "SSl [(1-alpha)*f*SSl] -> Il", # untreated relapse
                 "SSl [alpha*f*SSl] -> Tl", # treated relapse
                 "Il [gamma*Il] -> I0",  # liver clearance
                 "Tl [gamma*Tl] -> T0",  # liver clearance
                 "SSl [gamma*SSl] -> SS0",  # liver clearance
                 "Il [r*Il] -> SSl", # recovery
                 "Tl [r*Tl] -> SSl", # recovery
                 "I0 [r*I0] -> SS0", # recovery
                 "T0 [r*T0] -> SS0", # recovery
                 "Tl [(1-beta)*sigma*Tl] -> SSl", # treatment no rad cure
                 "Tl [beta*sigma*Tl] -> SS0", # treatment rad cure
                 "T0 [sigma*T0 ] -> SS0", # treatment
                 "Pl [gamma*Pl ] -> P0" # liver clearance in prophylaxis compartment
  )

  vivax_simu <- TiPS::build_simulator(reactions)

  safe_run <- function(f, ...) {
    out <- list()
    while(! length(out)) {out <- f(...)}
    out
  }

  safe_vivax_simu <- function(...) safe_run(vivax_simu, ...)

  reactions_hhl=reactions[c(1,2,4,5,13,14)]
  reactions_importation=reactions[c(7,8,10,11)]
  reactions_hh=c(reactions_hhl, reactions_importation)

  return(list(model=safe_vivax_simu, reactions_incidence=list(reactions_hh=reactions_hh, reactions_hhl=reactions_hhl)))
}



#' @title Simulate vivax model with delay to treatment and MDA, with stochastic model
#' @description Simulates a draw of the chosen vivax stochastic model
#' @param parameters model parameters
#' @param STOmodel  stochastic model to be simulated (reactions)
#' @param STOmodel_mda  stochastic model with MDA to be simulated (during prophylaxis)
#' @param maxtime maximal time step
#' @param year if TRUE, aggregates outputs per year. if FALSE, returns daily output
#' @param sto_method a scalarindicating which simulation method is used.
#' Default ("exact") is Gillespie's direct method. Other options are "approximate" (tau-leap) or "mixed".
#' cf. the documentation of the TiPS package for more information.
#' @param runs number of draws of the stochastic model
#' @param seeds a vector of the length of runs containing the seeds for each simulation (don't use "0" which has another use in TiPS)
#'
#' @return A dataframe containing the simulated model
#' @export
#' @importFrom rlang .data
simulate_vivax_delay_mda_sto=function(parameters, STOmodel=model_sto_vivax_delay(), STOmodel_mda=model_sto_vivax_delay_mda(), maxtime=1465, year=FALSE, sto_method="exact", runs=1, seeds=NULL){


  MDAcov=parameters["MDAcov"][[1]]
  MDArad_cure=parameters["MDArad_cure"][[1]]
  MDAp_length=parameters["MDAp_length"][[1]]


  # apply the MDA on initial conditions
  state_mda      = round(c(Il= parameters$Il*(1-MDAcov),
                       I0= parameters$I0*(1-MDAcov),
                       SSl= parameters$Sl*(1-MDAcov),
                       SS0 = parameters$S0*(1-MDAcov),
                       Tl= parameters$Tl*(1-MDAcov),
                       T0 = parameters$T0*(1-MDAcov),
                       Pl=MDAcov*(1-MDArad_cure)*(parameters$Il+parameters$Tl+parameters$Sl),
                       P0=MDAcov*(parameters$I0+ parameters$T0+parameters$S0+MDArad_cure*(parameters$Il+parameters$Tl+parameters$Sl))
  )*parameters$N)

  if(sum(state_mda)!=parameters$N){
    state_mda[which.max(state_mda)]=state_mda[which.max(state_mda)] - (sum(state_mda) - parameters$N)
  }


  times_mda <- c(0, MDAp_length)
  times_post_mda = c(MDAp_length, maxtime)

  # if omega is a function, extract values
  if(is.function(parameters$omega)){
    parameters$omega=parameters$omega(seq(0, maxtime-1))
    times_mda <- seq(0, MDAp_length)
    times_post_mda = seq(MDAp_length, maxtime, by = 1)
  }

  # if delta is a function, extract values
  if(is.function(parameters$delta)){
    parameters$delta=parameters$delta(seq(0, maxtime-1))
    times_mda <- seq(0, MDAp_length)
    times_post_mda = seq(MDAp_length, maxtime, by = 1)
  }

  safe_vivax_simu_mda=STOmodel_mda$model
  safe_vivax_simu_post_mda=STOmodel$model

  if(is.null(seeds)){seeds=1:runs}
  if(length(seeds)!=runs)(stop("seeds length does not match runs"))


  solutionVivax=data.frame()
  for (rr in 1:runs){

    #simulation during MDA prophylaxis
    traj_dm_mda <- safe_vivax_simu_mda(
      paramValues = parameters,
      initialStates = state_mda,
      times = times_mda, method=sto_method, seed=seeds[rr])

    this_solutionVivax_mda=traj_dm_mda$traj
    this_solutionVivax_mda$Sl=this_solutionVivax_mda$SSl+ this_solutionVivax_mda$Pl
    this_solutionVivax_mda$S0=this_solutionVivax_mda$SS0+ this_solutionVivax_mda$P0

    this_solutionVivax_mda=this_solutionVivax_mda[c("Time", "Reaction", "Nrep","Il","I0","Sl","S0","Tl","T0")]

    # simulation after prophylaxis
    state_post_mda= c(Il= this_solutionVivax_mda$Il[nrow(this_solutionVivax_mda)],
                      I0= this_solutionVivax_mda$I0[nrow(this_solutionVivax_mda)],
                      Sl= this_solutionVivax_mda$Sl[nrow(this_solutionVivax_mda)],
                      S0= this_solutionVivax_mda$S0[nrow(this_solutionVivax_mda)],
                      Tl= this_solutionVivax_mda$Tl[nrow(this_solutionVivax_mda)],
                      T0= this_solutionVivax_mda$T0[nrow(this_solutionVivax_mda)])


    traj_dm_post_mda <- safe_vivax_simu_post_mda(
      paramValues = parameters,
      initialStates = state_post_mda,
      times = times_post_mda, method=sto_method, seed=seeds[rr])

    this_solutionVivax_post_mda=traj_dm_post_mda$traj

    this_solutionVivax=rbind(this_solutionVivax_mda, this_solutionVivax_post_mda)
    this_solutionVivax$run=seeds[rr]
    solutionVivax=rbind(solutionVivax, this_solutionVivax)
  }

  solutionVivax$is_hh=solutionVivax$Reaction %in% c(STOmodel_mda$reactions_incidence$reactions_hh,STOmodel$reactions_incidence$reactions_hh)
  solutionVivax$is_hhl=solutionVivax$Reaction %in% c(STOmodel_mda$reactions_incidence$reactions_hhl,STOmodel$reactions_incidence$reactions_hhl)
  solutionVivax$is_rcd=solutionVivax$Reaction %in% c(STOmodel_mda$reactions_incidence$reactions_rcd,STOmodel$reactions_incidence$reactions_rcd)

  if(year){
    solutionVivax$time=ceiling(solutionVivax$Time /365)*365
    all_time_steps=unique(ceiling(seq(0, maxtime)/365))*365
  } else {
    solutionVivax$time=ceiling(solutionVivax$Time)
    all_time_steps=seq(0, maxtime)
  }
  solutionVivax=solutionVivax[solutionVivax$time <= maxtime,]
  all_time_steps=all_time_steps[all_time_steps <= maxtime]


  solutionVivax=solutionVivax  %>% dplyr::group_by(.data$time,.data$run) %>%
    dplyr::summarise(Il=.data$Il[.data$Time==max(.data$Time)],I0=.data$I0[.data$Time==max(.data$Time)],
                     Sl=.data$Sl[.data$Time==max(.data$Time)],S0=.data$S0[.data$Time==max(.data$Time)],
                     Tl=.data$Tl[.data$Time==max(.data$Time)],T0=.data$T0[.data$Time==max(.data$Time)],
                     hh=sum(.data$is_hh), hhl=sum(.data$is_hhl), rcd_reac=sum(.data$is_rcd))

  corr_factor_rcd=ifelse(any(solutionVivax$rcd_reac), (1-parameters$alpha-parameters$rho+parameters$rho*parameters$alpha/parameters$kappa)/(1-parameters$alpha),0)

  solutionVivax$h=solutionVivax$hh*parameters$rho + solutionVivax$rcd_reac*corr_factor_rcd
  solutionVivax$hl=solutionVivax$hhl*parameters$rho + solutionVivax$rcd_reac*corr_factor_rcd
  solutionVivax$I=solutionVivax$Il+solutionVivax$I0+solutionVivax$Tl+solutionVivax$T0

  solutionVivax$hh[solutionVivax$time==0]=ifelse(year, parameters$hh*365, parameters$hh)*parameters$N
  solutionVivax$hhl[solutionVivax$time==0]=ifelse(year, parameters$hhl*365, parameters$hhl)*parameters$N
  solutionVivax$h[solutionVivax$time==0]=ifelse(year, parameters$h*365, parameters$h)*parameters$N
  solutionVivax$hl[solutionVivax$time==0]=ifelse(year, parameters$hl*365, parameters$hl)*parameters$N

  solutionVivax$rcd_reac=NULL
  # add the lines for the time steps without reactions
  all_times=list()
  all_times$time=all_time_steps
  all_times$run=seeds
  all_times_full = base::expand.grid( all_times )

  solutionVivax_full=base::merge(solutionVivax, all_times_full, all=T)

  solutionVivax_full=solutionVivax_full[  with(solutionVivax_full, order(run, time)),]
  index_list=base::which(is.na(solutionVivax_full$Il))
  for(i in index_list){
    solutionVivax_full[i,c("Il", "I0", "Sl", "S0", "Tl", "T0", "I")]=solutionVivax_full[(i-1),c("Il", "I0", "Sl", "S0", "Tl", "T0", "I")]
    solutionVivax_full[i,c("hh", "hhl", "h", "hl")]=c(0,0,0,0)
  }


  # calculate proportion of imported cases
  solutionVivax_full$p=(solutionVivax_full$h-solutionVivax_full$hl)/solutionVivax_full$h

  return(as.data.frame(solutionVivax_full))
}



#' @title Reaction associated with the vivax model with delays and RCD, in stochastic version
#' @description This function creates the stochastic model for vivax dynamics with delays to treatment and reactive case detection (referral to health facilities).
#' @return A list with 2 objects:
#' model: the stochastic model compiled (as in the TiPS package)
#' reactions_incidence: a list of vectors indicating which reactions correspond to the incidence variables hh and hhl
#' @export
model_sto_vivax_delay_rcd_referral =function(){
  reactions <- c("S0 [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*S0/N] -> Il", # untreated infection
                 "Sl [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*Sl/N] -> Il", # untreated infection
                 "I0 [(omega*lambda*(Il+I0+Tl+T0))*I0/N] -> Il", # untreated reinfection
                 "S0 [alpha*(omega*lambda*(Il+I0+Tl+T0))*S0/N] -> Tl", # treated infection
                 "Sl [alpha*(omega*lambda*(Il+I0+Tl+T0))*Sl/N] -> Tl", # treated infection
                 "T0 [(omega*lambda*(Il+I0+Tl+T0))*T0/N] -> Tl", # treated reinfection

                 "S0 [(1-alpha)*(delta)*S0] -> Il", # untreated infection imported
                 "Sl [(1-alpha)*(delta)*Sl] -> Il", # untreated infection imported
                 "I0 [(delta)*I0] -> Il", # untreated reinfection imported
                 "S0 [alpha*(delta)*S0] -> Tl", # treated infection imported
                 "Sl [alpha*(delta)*Sl] -> Tl", # treated infection imported
                 "T0 [(delta)*T0] -> Tl", # treated reinfection imported

                 "Sl [(1-alpha)*f*Sl] -> Il", # untreated relapse
                 "Sl [alpha*f*Sl] -> Tl", # treated relapse
                 "Il [gamma*Il] -> I0",  # liver clearance
                 "Tl [gamma*Tl] -> T0",  # liver clearance
                 "Sl [gamma*Sl] -> S0",  # liver clearance
                 "Il [r*Il] -> Sl", # recovery
                 "Tl [r*Tl] -> Sl", # recovery
                 "I0 [r*I0] -> S0", # recovery
                 "T0 [r*T0] -> S0", # recovery
                 "Tl [(1-beta)*sigma*Tl] -> Sl", # treatment no rad cure
                 "Tl [beta*sigma*Tl] -> S0", # treatment rad cure
                 "T0 [sigma*T0 ] -> S0", # treatment

                 "Il [Il*nu*tau*eta*(!(iota<rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl)/N +f*Sl)/N)?(rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl)/N +f*Sl)/N):iota)] -> Tl", # RCD referral
                 "I0 [I0*nu*tau*eta*(!(iota<rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl)/N +f*Sl)/N)?(rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl)/N +f*Sl)/N):iota)] -> T0" # RCD referral



  )

  vivax_simu <- TiPS::build_simulator(reactions)

  safe_run <- function(f, ...) {
    out <- list()
    while(! length(out)) {out <- f(...)}
    out
  }

  safe_vivax_simu <- function(...) safe_run(vivax_simu, ...)

  reactions_hhl=reactions[c(1,2,4,5,13,14)]
  reactions_importation=reactions[c(7,8,10,11)]
  reactions_hh=c(reactions_hhl, reactions_importation)
  reactions_rcd=reactions[c(25,26)]

  return(list(model=safe_vivax_simu, reactions_incidence=list(reactions_hh=reactions_hh, reactions_hhl=reactions_hhl, reactions_rcd=reactions_rcd)))
}

#' @title Reaction associated with the vivax model with delays and RCD, in stochastic version
#' @description This function creates the stochastic model for vivax dynamics with delays to treatment and reactive case detection (without referral to health facilities).
#' @return A list with 2 objects:
#' model: the stochastic model compiled (as in the TiPS package)
#' reactions_incidence: a list of vectors indicating which reactions correspond to the incidence variables hh and hhl
#' @export
model_sto_vivax_delay_rcd_no_referral =function(){
  reactions <- c("S0 [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*S0/N] -> Il", # untreated infection
                 "Sl [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*Sl/N] -> Il", # untreated infection
                 "I0 [(omega*lambda*(Il+I0+Tl+T0))*I0/N] -> Il", # untreated reinfection
                 "S0 [alpha*(omega*lambda*(Il+I0+Tl+T0))*S0/N] -> Tl", # treated infection
                 "Sl [alpha*(omega*lambda*(Il+I0+Tl+T0))*Sl/N] -> Tl", # treated infection
                 "T0 [(omega*lambda*(Il+I0+Tl+T0))*T0/N] -> Tl", # treated reinfection

                 "S0 [(1-alpha)*(delta)*S0] -> Il", # untreated infection imported
                 "Sl [(1-alpha)*(delta)*Sl] -> Il", # untreated infection imported
                 "I0 [(delta)*I0] -> Il", # untreated reinfection imported
                 "S0 [alpha*(delta)*S0] -> Tl", # treated infection imported
                 "Sl [alpha*(delta)*Sl] -> Tl", # treated infection imported
                 "T0 [(delta)*T0] -> Tl", # treated reinfection imported

                 "Sl [(1-alpha)*f*Sl] -> Il", # untreated relapse
                 "Sl [alpha*f*Sl] -> Tl", # treated relapse
                 "Il [gamma*Il] -> I0",  # liver clearance
                 "Tl [gamma*Tl] -> T0",  # liver clearance
                 "Sl [gamma*Sl] -> S0",  # liver clearance
                 "Il [r*Il] -> Sl", # recovery
                 "Tl [r*Tl] -> Sl", # recovery
                 "I0 [r*I0] -> S0", # recovery
                 "T0 [r*T0] -> S0", # recovery
                 "Tl [(1-beta)*sigma*Tl] -> Sl", # treatment no rad cure
                 "Tl [beta*sigma*Tl] -> S0", # treatment rad cure
                 "T0 [sigma*T0 ] -> S0", # treatment

                 "Il [(1-beta)*Il*nu*tau*eta*(!(iota<rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl)/N +f*Sl)/N)?(rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl)/N +f*Sl)/N):iota)] -> Sl", # RCD referral
                 "Il [beta*Il*nu*tau*eta*(!(iota<rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl)/N +f*Sl)/N)?(rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl)/N +f*Sl)/N):iota)] -> S0", # RCD referral
                 "I0 [I0*nu*tau*eta*(!(iota<rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl)/N +f*Sl)/N)?(rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl)/N +f*Sl)/N):iota)] -> S0"

  )

  vivax_simu <- TiPS::build_simulator(reactions)

  safe_run <- function(f, ...) {
    out <- list()
    while(! length(out)) {out <- f(...)}
    out
  }

  safe_vivax_simu <- function(...) safe_run(vivax_simu, ...)

  reactions_hhl=reactions[c(1,2,4,5,13,14)]
  reactions_importation=reactions[c(7,8,10,11)]
  reactions_hh=c(reactions_hhl, reactions_importation)
  reactions_rcd=reactions[c(25,26, 27)]

  return(list(model=safe_vivax_simu, reactions_incidence=list(reactions_hh=reactions_hh, reactions_hhl=reactions_hhl, reactions_rcd=reactions_rcd)))
}



#' @title Reaction associated with the vivax model with delays, RCD and MDA during prophylaxis, in stochastic version
#' @description This function creates the stochastic model for vivax dynamics with delays to treatment, reactive case detection (referral to health facilities) and MDA during prophylaxis.
#' @return A list with 2 objects:
#' model: the stochastic model compiled (as in the TiPS package)
#' reactions_incidence: a list of vectors indicating which reactions correspond to the incidence variables hh and hhl
#' @export
model_sto_vivax_delay_rcd_referral_mda =function(){
  reactions <- c("SS0 [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*SS0/N] -> Il", # untreated infection
                 "SSl [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*SSl/N] -> Il", # untreated infection
                 "I0 [(omega*lambda*(Il+I0+Tl+T0))*I0/N] -> Il", # untreated reinfection
                 "SS0 [alpha*(omega*lambda*(Il+I0+Tl+T0))*SS0/N] -> Tl", # treated infection
                 "SSl [alpha*(omega*lambda*(Il+I0+Tl+T0))*SSl/N] -> Tl", # treated infection
                 "T0 [(omega*lambda*(Il+I0+Tl+T0))*T0/N] -> Tl", # treated reinfection

                 "SS0 [(1-alpha)*(delta)*SS0] -> Il", # untreated infection imported
                 "SSl [(1-alpha)*(delta)*SSl] -> Il", # untreated infection imported
                 "I0 [(delta)*I0] -> Il", # untreated reinfection imported
                 "SS0 [alpha*(delta)*SS0] -> Tl", # treated infection imported
                 "SSl [alpha*(delta)*SSl] -> Tl", # treated infection imported
                 "T0 [(delta)*T0] -> Tl", # treated reinfection imported

                 "SSl [(1-alpha)*f*SSl] -> Il", # untreated relapse
                 "SSl [alpha*f*SSl] -> Tl", # treated relapse
                 "Il [gamma*Il] -> I0",  # liver clearance
                 "Tl [gamma*Tl] -> T0",  # liver clearance
                 "SSl [gamma*SSl] -> SS0",  # liver clearance
                 "Il [r*Il] -> SSl", # recovery
                 "Tl [r*Tl] -> SSl", # recovery
                 "I0 [r*I0] -> SS0", # recovery
                 "T0 [r*T0] -> SS0", # recovery
                 "Tl [(1-beta)*sigma*Tl] -> SSl", # treatment no rad cure
                 "Tl [beta*sigma*Tl] -> SS0", # treatment rad cure
                 "T0 [sigma*T0 ] -> SS0", # treatment

                 "Il [Il*nu*tau*eta*(!(iota<rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(SS0+SSl)/N +f*SSl)/N)?(rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(SS0+SSl)/N +f*SSl)/N):iota)] -> Tl", # RCD referral
                 "I0 [I0*nu*tau*eta*(!(iota<rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(SS0+SSl)/N +f*SSl)/N)?(rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(SS0+SSl)/N +f*SSl)/N):iota)] -> T0", # RCD referral

                 "Pl [gamma*Pl ] -> P0" # liver clearance in prophylaxis compartment

  )

  vivax_simu <- TiPS::build_simulator(reactions)

  safe_run <- function(f, ...) {
    out <- list()
    while(! length(out)) {out <- f(...)}
    out
  }

  safe_vivax_simu <- function(...) safe_run(vivax_simu, ...)

  reactions_hhl=reactions[c(1,2,4,5,13,14)]
  reactions_importation=reactions[c(7,8,10,11)]
  reactions_hh=c(reactions_hhl, reactions_importation)
  reactions_rcd=reactions[c(25,26)]

  return(list(model=safe_vivax_simu, reactions_incidence=list(reactions_hh=reactions_hh, reactions_hhl=reactions_hhl, reactions_rcd=reactions_rcd)))
}

#' @title Reaction associated with the vivax model with delays, RCD and MDA during prophylaxis, in stochastic version
#' @description This function creates the stochastic model for vivax dynamics with delays to treatment, reactive case detection (referral to health facilities) and MDA during prophylaxis.
#' @return A list with 2 objects:
#' model: the stochastic model compiled (as in the TiPS package)
#' reactions_incidence: a list of vectors indicating which reactions correspond to the incidence variables hh and hhl
#' @export
model_sto_vivax_delay_rcd_no_referral_mda =function(){
  reactions <- c("SS0 [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*SS0/N] -> Il", # untreated infection
                 "SSl [(1-alpha)*(omega*lambda*(Il+I0+Tl+T0))*SSl/N] -> Il", # untreated infection
                 "I0 [(omega*lambda*(Il+I0+Tl+T0))*I0/N] -> Il", # untreated reinfection
                 "SS0 [alpha*(omega*lambda*(Il+I0+Tl+T0))*SS0/N] -> Tl", # treated infection
                 "SSl [alpha*(omega*lambda*(Il+I0+Tl+T0))*SSl/N] -> Tl", # treated infection
                 "T0 [(omega*lambda*(Il+I0+Tl+T0))*T0/N] -> Tl", # treated reinfection

                 "SS0 [(1-alpha)*(delta)*SS0] -> Il", # untreated infection imported
                 "SSl [(1-alpha)*(delta)*SSl] -> Il", # untreated infection imported
                 "I0 [(delta)*I0] -> Il", # untreated reinfection imported
                 "SS0 [alpha*(delta)*SS0] -> Tl", # treated infection imported
                 "SSl [alpha*(delta)*SSl] -> Tl", # treated infection imported
                 "T0 [(delta)*T0] -> Tl", # treated reinfection imported

                 "SSl [(1-alpha)*f*SSl] -> Il", # untreated relapse
                 "SSl [alpha*f*SSl] -> Tl", # treated relapse
                 "Il [gamma*Il] -> I0",  # liver clearance
                 "Tl [gamma*Tl] -> T0",  # liver clearance
                 "SSl [gamma*SSl] -> SS0",  # liver clearance
                 "Il [r*Il] -> SSl", # recovery
                 "Tl [r*Tl] -> SSl", # recovery
                 "I0 [r*I0] -> SS0", # recovery
                 "T0 [r*T0] -> SS0", # recovery
                 "Tl [(1-beta)*sigma*Tl] -> SSl", # treatment no rad cure
                 "Tl [beta*sigma*Tl] -> SS0", # treatment rad cure
                 "T0 [sigma*T0 ] -> SS0", # treatment

                 "Il [(1-beta)*Il*nu*tau*eta*(!(iota<rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(SS0+SSl)/N +f*SSl)/N)?(rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(SS0+SSl)/N +f*SSl)/N):iota)] -> SSl", # RCD referral
                 "Il [beta*Il*nu*tau*eta*(!(iota<rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(SS0+SSl)/N +f*SSl)/N)?(rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(SS0+SSl)/N +f*SSl)/N):iota)] -> SS0", # RCD referral
                 "I0 [I0*nu*tau*eta*(!(iota<rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(SS0+SSl)/N +f*SSl)/N)?(rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(SS0+SSl)/N +f*SSl)/N):iota)] -> SS0",

                 "Pl [gamma*Pl ] -> P0" # liver clearance in prophylaxis compartment

  )

  vivax_simu <- TiPS::build_simulator(reactions)

  safe_run <- function(f, ...) {
    out <- list()
    while(! length(out)) {out <- f(...)}
    out
  }

  safe_vivax_simu <- function(...) safe_run(vivax_simu, ...)

  reactions_hhl=reactions[c(1,2,4,5,13,14)]
  reactions_importation=reactions[c(7,8,10,11)]
  reactions_hh=c(reactions_hhl, reactions_importation)
  reactions_rcd=reactions[c(25,26,27)]

  return(list(model=safe_vivax_simu, reactions_incidence=list(reactions_hh=reactions_hh, reactions_hhl=reactions_hhl, reactions_rcd=reactions_rcd)))
}
