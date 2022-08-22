#' @title ODE vivax model with case management, importation and vector control, including delay to treatment
#' @description Defines the equations of the model
#' @param t time
#' @param y state variables
#' @param parameters model parameters
#'
#' @return a list of the state variables for all time steps
#' @export

ode_vivax_delay <- function(t, y, parameters) {
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

  Ul=y[1]
  U0=y[2]
  Sl=y[3]
  S0=y[4]
  Tl=y[5]
  T0=y[6]
  h=y[7]
  hl=y[8]
  hh=y[9]
  hhl=y[10]

  dUl= (1-alpha)*(omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) + (1-alpha)*f*Sl + (omega*lambda*(Ul+U0+Tl+T0)+delta)*U0 - gamma*Ul - r*Ul
  dU0= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(U0) +  gamma*Ul - r*U0
  dSl= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(Sl) - f*Sl +(1-beta)*sigma*Tl -  gamma*Sl + r*Ul+ r*Tl
  dS0= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0) +beta*sigma*Tl+ sigma *T0 +  gamma*Sl + r*U0+ r*T0
  dTl= alpha*(omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) + alpha*f*Sl - sigma*Tl -r*Tl -gamma*Tl+ (omega*lambda*(Ul+U0+Tl+T0)+delta)*T0
  dT0= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(T0) - sigma*T0+gamma*Tl - r*T0
  dh= rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl)
  dhl= rho*((omega*lambda*(Ul+U0+Tl+T0))*(S0+Sl) +f*Sl)
  dhh= ((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl)
  dhhl= ((omega*lambda*(Ul+U0+Tl+T0))*(S0+Sl) +f*Sl)
  res=c(dUl, dU0, dSl, dS0, dTl, dT0, dh, dhl, dhh, dhhl)
  return(list(res))
}


#' @title Simulate vivax model with delay to treatment
#' @description Simulates a draw of the chosen vivax ODE model
#' @param parameters model parameters
#' @param ODEmodel  ODE model to be simulated
#' @param maxtime maximal time step
#' @param year if TRUE, aggregates outputs per year. if FALSE, returns daily output
#'
#' @return A dataframe containing the simulated model
#' @export
simulate_vivax_delay_ode=function(parameters, ODEmodel=ode_vivax_delay, maxtime=1465, year=FALSE){

  # if omega is a scalar, make it a constant function
  if(is.numeric(parameters$omega)){
    omega_t <- stats::approxfun(data.frame(t=seq(0, maxtime+1)) %>% dplyr::mutate(vc=parameters$omega))
    parameters$omega=omega_t
  }
  # if delta is a scalar, make it a constant function
  if(is.numeric(parameters$delta)){
    delta_t <- stats::approxfun(data.frame(t=seq(0, maxtime+1)) %>% dplyr::mutate(vc=parameters$delta))
    parameters$delta=delta_t
  }


  #simulation
  state      <- c(parameters["Ul"][[1]],
                  parameters["U0"][[1]],
                  parameters["Sl"][[1]],
                  parameters["S0"][[1]],
                  parameters["Tl"][[1]],
                  parameters["T0"][[1]],
                  parameters["h"][[1]],
                  parameters["hl"][[1]],
                  parameters["hh"][[1]],
                  parameters["hhl"][[1]])
  times <- seq(0, maxtime, by = 1)
  solveSIS <- deSolve::ode(y = state, times = times, func = ODEmodel, parms =parameters, method = "lsoda")
  solutionVivax=as.data.frame(solveSIS)
  names(solutionVivax)=c("time", "Ul", "U0", "Sl", "S0", "Tl", "T0", "h", "hl", "hh", "hhl")
  if(year){solutionVivax=solutionVivax[(solutionVivax$time %%365 ==0),]}
  h0=ifelse(year, parameters["h"][[1]]*365, parameters["h"][[1]])
  solutionVivax$h=c(h0,diff(solutionVivax$h))
  hl0=ifelse(year, parameters["hl"][[1]]*365, parameters["hl"][[1]])
  solutionVivax$hl=c(hl0,diff(solutionVivax$hl))
  hh0=ifelse(year, parameters["hh"][[1]]*365, parameters["hh"][[1]])
  solutionVivax$hh=c(hh0,diff(solutionVivax$hh))
  hhl0=ifelse(year, parameters["hhl"][[1]]*365, parameters["hhl"][[1]])
  solutionVivax$hhl=c(hhl0,diff(solutionVivax$hhl))

  solutionVivax$I=solutionVivax$Ul+solutionVivax$U0+solutionVivax$Tl+solutionVivax$T0
  solutionVivax$p=(solutionVivax$h-solutionVivax$hl)/solutionVivax$h

  return(solutionVivax)
}




#' @title Rc calculation, with case management and vector control and delay to treatment
#'
#' @description Calculation of Rc in vivax model
#'
#' @param lambda transmission parameter
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param sigma delay to treatment (1/duration of the delay)
#' @param omega intensity of vector control. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#'
#' @return A scalar Rc value
#' @export
#'
get_rc_vivax_delay <- function(lambda, f, r, gamma, alpha, beta, sigma, omega=1){
  mynum=omega*lambda*(gamma+r)*(f+gamma)*(gamma+r+sigma)*(r+sigma*(1-alpha))
  mydenom=r*(r+sigma)*(gamma*(f+gamma+r)*(gamma+r+sigma)+alpha*f*sigma*(beta*(r+gamma)-gamma))
  return(mynum/mydenom)
}


get_lambda_from_rc_vivax_delay <- function(Rc, f, r, gamma, alpha, beta, sigma, omega=1){
  mynum=omega*(gamma+r)*(f+gamma)*(gamma+r+sigma)*(r+sigma*(1-alpha))
  mydenom=r*(r+sigma)*(gamma*(f+gamma+r)*(gamma+r+sigma)+alpha*f*sigma*(beta*(r+gamma)-gamma))
  return(Rc*mydenom/mynum)
}


#' @title Calculate equilibrium state variables, in model with delay to treatment
#'
#' @description Calculates equilibrium state variables U0, Ul, T0, Tl, S0 and Sl based on I and model parameters
#'
#' @param I proportion of infectious individuals at equilibrium (U0+Ul)
#' @param lambda transmission rate
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param sigma delay to treatment (1/duration of the delay)
#' @param rho observation rate
#' @param omega vector control intensity. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#' @param delta importation rate
#'
#' @return A list with the equilibrium states (U0, Ul, S0 and Sl)
#' @export
#'
get_equilibrium_states_vivax_delay <- function(I, lambda, r, gamma, f, alpha, beta, sigma, rho, delta, omega){
  lambda=lambda*omega
  TT=I*alpha*r/(r+(1-alpha)*sigma)  # T0+Tl
  UU=I*(1-alpha)*(r+sigma)/(r+(1-alpha)*sigma) #U0 + ul
  T0=TT*gamma/(lambda*I+delta+gamma+r+sigma)
  Tl=TT-T0
  U0=UU*gamma/(lambda*I+delta+gamma+r)
  Ul=UU-U0
  Sl=(r*Ul+(r+(1-beta)*sigma)*Tl)/(lambda*I+delta+gamma+f)
  S0=1-I-Sl
  h=rho*((lambda*I+delta)*(1-I)+f*Sl)
  hl=rho*((lambda*I)*(1-I)+f*Sl)
  hh=(lambda*I+delta)*(1-I)+f*Sl
  hhl=(lambda*I)*(1-I)+f*Sl
  return(list("Ul"=Ul, "U0"=U0, "Tl"=Tl, "T0"=T0, "Sl"=Sl, "S0"=S0, "h"=h, "hl"=hl, "hh"=hh, "hhl"=hhl))
}


##################################
# Back-calculate lambda
##################################

#' @title Transmission rate calculation, in model with delay to treatment
#'
#' @description Solve the equation for lambda, vivax model with CM, vector control and importation
#'
#' @param h daily incidence
#' @param p proportion of imported cases
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param sigma delay to treatment (1/duration of the delay)
#' @param rho observation rate
#' @param omega vector control intensity. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#' @param return.all if TRUE, returns lambda, I and delta. If FALSE, returns only lambda
#'
#' @return If return.all=TRUE, returns a list with lambda, I and delta. If return.all=FALSE, returns only a lambda value
#'
#' @export
#'
solve_lambda_vivax_delay <- function(h, r,  f, gamma, alpha, beta, rho, sigma, p, omega, return.all=FALSE){

  a_1 <- function(I, r,  f, gamma, alpha, beta, rho,sigma,  delta){
    return(I^4*(1-I))
  }

  a_2 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta){
    return((I^3)*(1-I)*(3*gamma+ 4*delta + 2*r +sigma+ f) - r*(I^4)*(r+sigma)/(r+(1-alpha)*sigma))
  }

  a_3 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta){
    return((I^2)*(1-I)*((gamma +delta + r)*(3*gamma +3*delta +r+ 2*f)+delta*(3*gamma+2*r+f+3*delta)+sigma*(3*delta+2*gamma+r+f) )
           - r*(I^3)*((3*gamma +2*r+sigma +3*delta)*(r+sigma)+alpha*beta*f*sigma)/(r+(1-alpha)*sigma))
  }

  a_4 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta){
    return(I*(1-I)*((delta+gamma+f)*(delta+gamma+r)*(delta+gamma+r+sigma)+
                      delta*((delta+gamma+f)*(delta+gamma+r)+(delta+gamma+r)*(delta+gamma+r+sigma)+(delta+gamma+f)*(delta+gamma+r+sigma)))+
             I*I*r*(f*(r+(1-alpha*beta)*sigma)*(2*delta+2*r+sigma+gamma)-(r+sigma)*((delta+gamma+f)*(delta+gamma+r)+(delta+gamma+r)*(delta+gamma+r+sigma)+(delta+gamma+f)*(delta+gamma+r+sigma)))/(r+(1-alpha)*sigma))

  }

  a_5 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta){
    return((1-I)*delta*(delta+gamma+f)*(delta+gamma+r)*(delta+gamma+r+sigma)-
             I*r*(delta+gamma+r)*((delta+gamma)*(delta+gamma+r)*(r+sigma)+f*gamma*(r+sigma)+alpha*beta*f*sigma*(delta+r))/(r+(1-alpha)*sigma)-
             I*r*sigma*((delta+gamma)*(delta+gamma+r)*(r+sigma)+(1-alpha)*f*gamma*(r+sigma)+alpha*beta*f*sigma*(delta+gamma+r))/(r+(1-alpha)*sigma))
  }


  I=h*(r+(1-alpha)*sigma)/r/rho/(r+sigma)
  delta=p * h/(1-I)/rho
  # solving the equation for lambda
  lambda_complex =polyroot(c(a_5(I, r,  f, gamma, alpha, beta, rho, sigma, delta),
                                       a_4(I, r,  f, gamma, alpha, beta, rho, sigma, delta),
                                       a_3(I, r,  f, gamma, alpha, beta, rho, sigma, delta),
                                       a_2(I, r,  f, gamma, alpha, beta, rho, sigma, delta),
                                       a_1(I, r,  f, gamma, alpha, beta, rho, sigma, delta)))


  lambda_real=as.double(lambda_complex[abs(Im(lambda_complex)) < 1e-10])
  lambda_real=lambda_real[lambda_real>=0]

  if(return.all){
    return(list("lambda_real"=lambda_real/omega, "I"=I, "delta"=delta))
  } else{
    return(lambda_real/omega)
  }

}

