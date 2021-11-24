#' @title ODE vivax model with case management, importation and vector control, during MDA prophylaxis
#' @description Defines the equations of the model
#' @param t time
#' @param y state variables
#' @param parameters model parameters
#'
#' @return a list of the state variables for all time steps
#' @export
ode_vivax_mda <- function(t, y, parameters) {
  lambda=parameters["lambda"][[1]]
  f=parameters["f"][[1]]
  gamma=parameters["gamma"][[1]]
  r=parameters["r"][[1]]
  alpha=parameters["alpha"][[1]]
  beta=parameters["beta"][[1]]
  rho=parameters["rho"][[1]]
  delta=parameters["delta"][[1]](t)
  omega=parameters["omega"][[1]](t)

  Il=y[1]
  I0=y[2]
  Sl=y[3]
  S0=y[4]
  Pl=y[5]
  P0=y[6]
  h=y[7]
  hr=y[8]
  hl=y[9]
  hh=y[10]
  hhl=y[11]


  dIl= (1-alpha)*(omega*lambda*(Il+I0)+delta)*(S0+Sl) + (omega*lambda*(Il+I0)+delta)*I0 + (1-alpha)*f*Sl - gamma*Il - r*Il
  dI0= -(omega*lambda*(Il+I0)+delta)*(I0) +  gamma*Il - r*I0
  dSl= -(1-alpha*(1-beta))*(omega*lambda*(Il+I0)+delta)*(Sl) - (1-alpha*(1-beta))*f*Sl +(omega*lambda*(Il+I0)+delta)*alpha*(1-beta)*S0 -  gamma*Sl + r*Il
  dS0= -(1-alpha*beta)*(omega*lambda*(Il+I0)+delta)*(S0) + (omega*lambda*(Il+I0)+delta)*alpha*beta*Sl +alpha*beta*f*Sl +  gamma*Sl + r*I0
  dPl= -gamma*Pl
  dP0= gamma*Pl
  dh= rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl)
  dhr= rho*f*Sl
  dhl= rho*((omega*lambda*(Il+I0))*(S0+Sl) + f*Sl)
  dhh= ((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl)
  dhhl=((omega*lambda*(Il+I0))*(S0+Sl) + f*Sl)

  res=c(dIl, dI0, dSl, dS0, dPl, dP0, dh, dhr, dhl, dhh, dhhl)

  return(list(res))
}


#' @title ODE vivax model with case management, RCD, importation and vector control, during MDA prophylaxis
#' @description Defines the equations of the model
#' @param t time
#' @param y state variables
#' @param parameters model parameters
#'
#' @return a list of the state variables for all time steps
#' @export
ode_vivax_rcd_mda <- function(t, y, parameters) {
  lambda=parameters["lambda"][[1]]
  f=parameters["f"][[1]]
  gamma=parameters["gamma"][[1]]
  r=parameters["r"][[1]]
  alpha=parameters["alpha"][[1]]
  beta=parameters["beta"][[1]]
  rho=parameters["rho"][[1]]
  kappa=parameters["kappa"][[1]]
  delta=parameters["delta"][[1]](t)
  omega=parameters["omega"][[1]](t)
  iota=parameters["iota"][[1]]
  nu=parameters["nu"][[1]]
  eta=parameters["eta"][[1]]
  tau0=parameters["tau"][[1]]

  if(is.numeric(tau0)){
    tau=function(pr){return(tau0)}
  } else{ tau=tau0}

  corr_factor_rcd=(1-alpha-rho+rho*alpha/kappa)/(1-alpha)

  Il=y[1]
  I0=y[2]
  Sl=y[3]
  S0=y[4]
  Pl=y[5]
  P0=y[6]
  h=y[7]
  hr=y[8]
  hl=y[9]
  hh=y[10]
  hhl=y[11]

  dIl= (1-alpha)*(omega*lambda*(Il+I0)+delta)*(S0+Sl) + (omega*lambda*(Il+I0)+delta)*I0 + (1-alpha)*f*Sl - gamma*Il - r*Il -Il*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl) )
  dI0= -(omega*lambda*(Il+I0)+delta)*(I0) +  gamma*Il - r*I0 -I0*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl) )
  dSl= -(1-alpha*(1-beta))*(omega*lambda*(Il+I0)+delta)*(Sl) - (1-alpha*(1-beta))*f*Sl +(omega*lambda*(Il+I0)+delta)*alpha*(1-beta)*S0 -  gamma*Sl + r*Il +(1-beta)*Il*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl) )
  dS0= -(1-alpha*beta)*(omega*lambda*(Il+I0)+delta)*(S0) + (omega*lambda*(Il+I0)+delta)*alpha*beta*Sl +alpha*beta*f*Sl +  gamma*Sl + r*I0 + (beta*Il+I0)*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl) )
  dPl= -gamma*Pl
  dP0= gamma*Pl
  dh= rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl)+(Il+I0)*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl) )*corr_factor_rcd
  dhr= rho*f*Sl
  dhl= rho*((omega*lambda*(Il+I0))*(S0+Sl) + f*Sl)+ (Il+I0)*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl))*corr_factor_rcd
  dhh= ((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl)
  dhhl= ((omega*lambda*(Il+I0))*(S0+Sl) + f*Sl)

  res=c(dIl, dI0, dSl, dS0, dPl, dP0, dh, dhr, dhl, dhh, dhhl)

  return(list(res))
}

#' @title Simulate vivax model with MDA
#' @description Simulates a draw of the chosen vivax ODE model
#' @param parameters model parameters
#' @param ODEmodel  ODE model to be simulated
#' @param ODEmodel_mda  ODE model with MDA to be simulated (during prophylaxis)
#' @param maxtime maximal time step
#' @param year if TRUE, aggregates outputs per year. if FALSE, returns daily output
#'
#' @return A dataframe containing the simulated model
#' @export
#'
simulate_vivax_mda_ode <- function(parameters, ODEmodel=ode_vivax_cm_importation_vc, ODEmodel_mda=ode_vivax_mda, maxtime=1465, year=FALSE){

  # if omega is a scalar, make it a constant function
  if(is.numeric(parameters$omega)){
    omega_t <- stats::approxfun(data.frame(t=seq(0, maxtime+1)) %>% dplyr::mutate(vc=parameters$omega))
    parameters$omega=omega_t
  }
  # if delta is a scalar, make it a constant function
  if(is.numeric(parameters$delta)){
    delta_t <- stats::approxfun(data.frame(t=seq(0, maxtime+1)) %>% dplyr::mutate(val=parameters$delta))
    parameters$delta=delta_t
  }

  MDAcov=parameters["MDAcov"][[1]]
  MDArad_cure=parameters["MDArad_cure"][[1]]
  MDAp_length=parameters["MDAp_length"][[1]]

  #simulation during MDA prophylaxis
  state_mda   = c(parameters["Il"][[1]]*(1-MDAcov),
                  parameters["I0"][[1]]*(1-MDAcov),
                  parameters["Sl"][[1]]*(1-MDAcov),
                  parameters["S0"][[1]]*(1-MDAcov),
                  MDAcov*(1-MDArad_cure)*(parameters["Il"][[1]]+parameters["Sl"][[1]]),
                  MDAcov*(parameters["I0"][[1]]+parameters["S0"][[1]]+MDArad_cure*(parameters["Il"][[1]]+parameters["Sl"][[1]])),
                  parameters["h"][[1]],
                  parameters["hr"][[1]],
                  parameters["hl"][[1]],
                  parameters["hh"][[1]],
                  parameters["hhl"][[1]])

  times_mda = seq(0, MDAp_length, by = 1)

  solveSIS_mda = deSolve::ode(y = state_mda, times = times_mda, func = ODEmodel_mda, parms =parameters, method = "lsoda")
  solutionVivax_mda=as.data.frame(solveSIS_mda)
  names(solutionVivax_mda)=c("time", "Il", "I0", "SSl", "SS0", "Pl", "P0", "h", "hr", "hl", "hh", "hhl")

  solutionVivax_mda$Sl=solutionVivax_mda$SSl+ solutionVivax_mda$Pl
  solutionVivax_mda$S0=solutionVivax_mda$SS0+ solutionVivax_mda$P0

  #simulation after MDA has vanished
  state_post_mda= c(solutionVivax_mda$Il[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$I0[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$Sl[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$S0[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$h[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$hr[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$hl[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$hh[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$hhl[solutionVivax_mda$time==MDAp_length])

  times_post_mda = seq(MDAp_length, maxtime, by = 1)

  solveSIS_post_mda = deSolve::ode(y = state_post_mda, times = times_post_mda, func = ODEmodel, parms =parameters, method = "lsoda")
  solutionVivax_post_mda=as.data.frame(solveSIS_post_mda)
  names(solutionVivax_post_mda)=c("time", "Il", "I0", "Sl", "S0", "h", "hr", "hl", "hh", "hhl")

  # Combine all results into one single file
  solutionVivax=rbind(solutionVivax_mda[c("time", "Il", "I0", "Sl", "S0", "h", "hr", "hl", "hh", "hhl")], solutionVivax_post_mda[solutionVivax_post_mda$time>MDAp_length,])

  if(year){solutionVivax=solutionVivax[(solutionVivax$time %%365 ==0),]}
  h0=ifelse(year, parameters["h"][[1]]*365, parameters["h"][[1]])
  hr0=ifelse(year, parameters["hr"][[1]]*365, parameters["hr"][[1]])
  solutionVivax$h=c(h0,diff(solutionVivax$h))
  solutionVivax$hr=c(hr0,diff(solutionVivax$hr))
  hl0=ifelse(year, parameters["hl"][[1]]*365, parameters["hl"][[1]])
  solutionVivax$hl=c(hl0,diff(solutionVivax$hl))
  hh0=ifelse(year, parameters["hh"][[1]]*365, parameters["hh"][[1]])
  solutionVivax$hh=c(hh0,diff(solutionVivax$hh))
  hhl0=ifelse(year, parameters["hhl"][[1]]*365, parameters["hhl"][[1]])
  solutionVivax$hhl=c(hhl0,diff(solutionVivax$hhl))

  solutionVivax$I=solutionVivax$Il+solutionVivax$I0
  solutionVivax$p=(solutionVivax$h-solutionVivax$hl)/solutionVivax$h

  return(solutionVivax)
}





#' @title ODE vivax model with case management, importation and vector control, including delay to treatment, during MDA prophylaxis
#' @description Defines the equations of the model
#' @param t time
#' @param y state variables
#' @param parameters model parameters
#'
#' @return a list of the state variables for all time steps
#' @export

ode_vivax_delay_mda <- function(t, y, parameters) {
  lambda=parameters["lambda"][[1]]
  f=parameters["f"][[1]]
  gamma=parameters["gamma"][[1]]
  r=parameters["r"][[1]]
  alpha=parameters["alpha"][[1]]
  beta=parameters["beta"][[1]]
  rho=parameters["rho"][[1]]
  sigma=parameters["sigma"][[1]]
  omega=parameters["omega"][[1]](t)
  delta=parameters["delta"][[1]](t)

  Il=y[1]
  I0=y[2]
  Sl=y[3]
  S0=y[4]
  Tl=y[5]
  T0=y[6]
  Pl=y[7]
  P0=y[8]
  h=y[9]
  hl=y[10]
  hh=y[11]
  hhl=y[12]

  dIl= (1-alpha)*(omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) + (1-alpha)*f*Sl + (omega*lambda*(Il+I0+Tl+T0)+delta)*I0 - gamma*Il - r*Il
  dI0= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(I0) +  gamma*Il - r*I0
  dSl= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(Sl) - f*Sl +(1-beta)*sigma*Tl -  gamma*Sl + r*Il+ r*Tl
  dS0= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(S0) +beta*sigma*Tl+ sigma *T0 +  gamma*Sl + r*I0+ r*T0
  dTl= alpha*(omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) + alpha*f*Sl - sigma*Tl -r*Tl -gamma*Tl+ (omega*lambda*(Il+I0+Tl+T0)+delta)*T0
  dT0= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(T0) - sigma*T0+gamma*Tl - r*T0
  dPl= -gamma*Pl
  dP0= gamma*Pl
  dh= rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl)
  dhl= rho*((omega*lambda*(Il+I0+Tl+T0))*(S0+Sl) +f*Sl)
  dhh= ((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl)
  dhhl= ((omega*lambda*(Il+I0+Tl+T0))*(S0+Sl) +f*Sl)

  res=c(dIl, dI0, dSl, dS0, dTl, dT0, dPl, dP0, dh, dhl, dhh, dhhl)

  return(list(res))
}




#' @title ODE vivax model with case management, RCD with referral, importation and vector control, including delay to treatment, during PDA prophylaxis
#' @description Defines the equations of the model
#' @param t time
#' @param y state variables
#' @param parameters model parameters
#'
#' @return a list of the state variables for all time steps
#' @export
ode_vivax_delay_rcd_referral_mda <- function(t, y, parameters) {
  lambda=parameters["lambda"][[1]]
  f=parameters["f"][[1]]
  gamma=parameters["gamma"][[1]]
  r=parameters["r"][[1]]
  alpha=parameters["alpha"][[1]]
  beta=parameters["beta"][[1]]
  rho=parameters["rho"][[1]]
  sigma=parameters["sigma"][[1]]
  omega=parameters["omega"][[1]](t)
  delta=parameters["delta"][[1]](t)
  iota=parameters["iota"][[1]]
  nu=parameters["nu"][[1]]
  eta=parameters["eta"][[1]]
  kappa=parameters["kappa"][[1]]
  tau0=parameters["tau"][[1]]

  if(is.numeric(tau0)){
    tau=function(pr){return(tau0)}
  } else{ tau=tau0}

  corr_factor_rcd=(1-alpha-rho+rho*alpha/kappa)/(1-alpha)

  Il=y[1]
  I0=y[2]
  Sl=y[3]
  S0=y[4]
  Tl=y[5]
  T0=y[6]
  Pl=y[7]
  P0=y[8]
  h=y[9]
  hl=y[10]
  hh=y[11]
  hhl=y[12]

  dIl= (1-alpha)*(omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) + (1-alpha)*f*Sl + (omega*lambda*(Il+I0+Tl+T0)+delta)*I0 - gamma*Il - r*Il -Il*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dI0= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(I0) +  gamma*Il - r*I0 -I0*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dSl= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(Sl) - f*Sl +(1-beta)*sigma*Tl -  gamma*Sl + r*Il+ r*Tl
  dS0= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(S0) +beta*sigma*Tl+ sigma *T0 +  gamma*Sl + r*I0+ r*T0
  dTl= alpha*(omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) + alpha*f*Sl - sigma*Tl -r*Tl -gamma*Tl+ (omega*lambda*(Il+I0+Tl+T0)+delta)*T0 +Il*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dT0= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(T0) - sigma*T0+gamma*Tl - r*T0 +I0*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dPl= -gamma*Pl
  dP0= gamma*Pl
  dh= rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl) + (Il+I0)*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))*corr_factor_rcd
  dhl= rho*((omega*lambda*(Il+I0+Tl+T0))*(S0+Sl) +f*Sl)+ (Il+I0)*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))*corr_factor_rcd
  dhh= ((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl)
  dhhl= ((omega*lambda*(Il+I0+Tl+T0))*(S0+Sl) +f*Sl)

  res=c(dIl, dI0, dSl, dS0, dTl, dT0, dPl, dP0, dh, dhl, dhh, dhhl)
  return(list(res))
}



#' @title ODE vivax model with case management, RCD without referral, importation and vector control, including delay to treatment, during MDA prophylaxis
#' @description Defines the equations of the model
#' @param t time
#' @param y state variables
#' @param parameters model parameters
#'
#' @return a list of the state variables for all time steps
#' @export
ode_vivax_delay_rcd_no_referral_mda <- function(t, y, parameters) {
  lambda=parameters["lambda"][[1]]
  f=parameters["f"][[1]]
  gamma=parameters["gamma"][[1]]
  r=parameters["r"][[1]]
  alpha=parameters["alpha"][[1]]
  beta=parameters["beta"][[1]]
  rho=parameters["rho"][[1]]
  sigma=parameters["sigma"][[1]]
  omega=parameters["omega"][[1]](t)
  delta=parameters["delta"][[1]](t)
  iota=parameters["iota"][[1]]
  nu=parameters["nu"][[1]]
  eta=parameters["eta"][[1]]
  kappa=parameters["kappa"][[1]]
  tau0=parameters["tau"][[1]]

  if(is.numeric(tau0)){
    tau=function(pr){return(tau0)}
  } else{ tau=tau0}

  corr_factor_rcd=(1-alpha-rho+rho*alpha/kappa)/(1-alpha)

  Il=y[1]
  I0=y[2]
  Sl=y[3]
  S0=y[4]
  Tl=y[5]
  T0=y[6]
  Pl=y[7]
  P0=y[8]
  h=y[9]
  hl=y[10]
  hh=y[11]
  hhl=y[12]

  dIl= (1-alpha)*(omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) + (1-alpha)*f*Sl + (omega*lambda*(Il+I0+Tl+T0)+delta)*I0 - gamma*Il - r*Il -Il*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dI0= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(I0) +  gamma*Il - r*I0 -I0*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dSl= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(Sl) - f*Sl +(1-beta)*sigma*Tl -  gamma*Sl + r*Il+ r*Tl +(1-beta)*Il*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dS0= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(S0) +beta*sigma*Tl+ sigma *T0 +  gamma*Sl + r*I0+ r*T0 +(beta*Il+I0)*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dTl= alpha*(omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) + alpha*f*Sl - sigma*Tl -r*Tl -gamma*Tl+ (omega*lambda*(Il+I0+Tl+T0)+delta)*T0
  dT0= -(omega*lambda*(Il+I0+Tl+T0)+delta)*(T0) - sigma*T0+gamma*Tl - r*T0
  dPl= -gamma*Pl
  dP0= gamma*Pl
  dh= rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl) + (Il+I0)*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))*corr_factor_rcd
  dhl= rho*((omega*lambda*(Il+I0+Tl+T0))*(S0+Sl) +f*Sl)+ (Il+I0)*nu*tau(pr=(I0+Il+T0+Tl))*eta*min(iota,rho*((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl))*corr_factor_rcd
  dhh= ((omega*lambda*(Il+I0+Tl+T0)+delta)*(S0+Sl) +f*Sl)
  dhhl= ((omega*lambda*(Il+I0+Tl+T0))*(S0+Sl) +f*Sl)

  res=c(dIl, dI0, dSl, dS0, dTl, dT0, dPl, dP0, dh, dhl, dhh, dhhl)
  return(list(res))
}

#' @title Simulate vivax model  with delay to treatment and MDA
#' @description Simulates a draw of the chosen vivax ODE model
#' @param parameters model parameters
#' @param ODEmodel  ODE model to be simulated
#' @param ODEmodel_mda  ODE model with MDA to be simulated (during prophylaxis)
#' @param maxtime maximal time step
#' @param year if TRUE, aggregates outputs per year. if FALSE, returns daily output
#'
#' @return A dataframe containing the simulated model
#' @export
simulate_vivax_delay_mda_ode=function(parameters, ODEmodel=ode_vivax_delay, ODEmodel_mda=ode_vivax_delay_mda, maxtime=1465, year=FALSE){

  # if omega is a scalar, make it a constant function
  if(is.numeric(parameters$omega)){
    omega_t <- stats::approxfun(data.frame(t=seq(0, maxtime+1)) %>% dplyr::mutate(vc=parameters$omega))
    parameters$omega=omega_t
  }
  # if delta is a scalar, make it a constant function
  if(is.numeric(parameters$delta)){
    delta_t <- stats::approxfun(data.frame(t=seq(0, maxtime+1)) %>% dplyr::mutate(val=parameters$delta))
    parameters$delta=delta_t
  }

  MDAcov=parameters["MDAcov"][[1]]
  MDArad_cure=parameters["MDArad_cure"][[1]]
  MDAp_length=parameters["MDAp_length"][[1]]

  #simulation during MDA prophylaxis
  state_mda   = c(parameters["Il"][[1]]*(1-MDAcov),
                  parameters["I0"][[1]]*(1-MDAcov),
                  parameters["Sl"][[1]]*(1-MDAcov),
                  parameters["S0"][[1]]*(1-MDAcov),
                  parameters["Tl"][[1]]*(1-MDAcov),
                  parameters["T0"][[1]]*(1-MDAcov),
                  MDAcov*(1-MDArad_cure)*(parameters["Il"][[1]]+parameters["Tl"][[1]]+parameters["Sl"][[1]]),
                  MDAcov*(parameters["I0"][[1]]+parameters["T0"][[1]]+parameters["S0"][[1]]+MDArad_cure*(parameters["Il"][[1]]+parameters["Tl"][[1]]+parameters["Sl"][[1]])),
                  parameters["h"][[1]],
                  parameters["hl"][[1]],
                  parameters["hh"][[1]],
                  parameters["hhl"][[1]])

  times_mda = seq(0, MDAp_length, by = 1)

  solveSIS_mda = deSolve::ode(y = state_mda, times = times_mda, func = ODEmodel_mda, parms =parameters, method = "lsoda")
  solutionVivax_mda=as.data.frame(solveSIS_mda)
  names(solutionVivax_mda)=c("time", "Il", "I0", "SSl", "SS0", "Tl", "T0", "Pl", "P0", "h", "hl", "hh", "hhl")

  solutionVivax_mda$Sl=solutionVivax_mda$SSl+ solutionVivax_mda$Pl
  solutionVivax_mda$S0=solutionVivax_mda$SS0+ solutionVivax_mda$P0

  #simulation after MDA has vanished
  state_post_mda= c(solutionVivax_mda$Il[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$I0[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$Sl[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$S0[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$Tl[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$T0[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$h[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$hl[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$hh[solutionVivax_mda$time==MDAp_length],
                    solutionVivax_mda$hhl[solutionVivax_mda$time==MDAp_length])

  times_post_mda = seq(MDAp_length, maxtime, by = 1)

  solveSIS_post_mda = deSolve::ode(y = state_post_mda, times = times_post_mda, func = ODEmodel, parms =parameters, method = "lsoda")
  solutionVivax_post_mda=as.data.frame(solveSIS_post_mda)
  names(solutionVivax_post_mda)=c("time", "Il", "I0", "Sl", "S0", "Tl", "T0", "h", "hl", "hh", "hhl")

  # Combine all results into one single file
  solutionVivax=rbind(solutionVivax_mda[c("time", "Il", "I0", "Sl", "S0", "Tl", "T0", "h", "hl", "hh", "hhl")], solutionVivax_post_mda[solutionVivax_post_mda$time>MDAp_length,])

  if(year){solutionVivax=solutionVivax[(solutionVivax$time %%365 ==0),]}
  h0=ifelse(year, parameters["h"][[1]]*365, parameters["h"][[1]])
  solutionVivax$h=c(h0,diff(solutionVivax$h))
  hl0=ifelse(year, parameters["hl"][[1]]*365, parameters["hl"][[1]])
  solutionVivax$hl=c(hl0,diff(solutionVivax$hl))
  hh0=ifelse(year, parameters["hh"][[1]]*365, parameters["hh"][[1]])
  solutionVivax$hh=c(hh0,diff(solutionVivax$hh))
  hhl0=ifelse(year, parameters["hhl"][[1]]*365, parameters["hhl"][[1]])
  solutionVivax$hhl=c(hhl0,diff(solutionVivax$hhl))

  solutionVivax$I=solutionVivax$Il+solutionVivax$I0+solutionVivax$Tl+solutionVivax$T0
  solutionVivax$p=(solutionVivax$h-solutionVivax$hl)/solutionVivax$h

  return(solutionVivax)
}
