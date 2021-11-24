#' @title ODE vivax model with case management, importation and vector control
#' @description Defines the equations of the model
#' @param t time
#' @param y state variables
#' @param parameters model parameters
#'
#' @return a list of the state variables for all time steps
#' @export
#'
ode_vivax_cm_importation_vc <- function(t, y, parameters) {
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
  h=y[5]
  hr=y[6]
  hl=y[7]
  hh=y[8]
  hhl=y[9]

  dIl= (1-alpha)*(omega*lambda*(Il+I0)+delta)*(S0+Sl) + (omega*lambda*(Il+I0)+delta)*I0 + (1-alpha)*f*Sl - gamma*Il - r*Il
  dI0= -(omega*lambda*(Il+I0)+delta)*(I0) +  gamma*Il - r*I0
  dSl= -(1-alpha*(1-beta))*(omega*lambda*(Il+I0)+delta)*(Sl) - (1-alpha*(1-beta))*f*Sl +(omega*lambda*(Il+I0)+delta)*alpha*(1-beta)*S0 -  gamma*Sl + r*Il
  dS0= -(1-alpha*beta)*(omega*lambda*(Il+I0)+delta)*(S0) + (omega*lambda*(Il+I0)+delta)*alpha*beta*Sl +alpha*beta*f*Sl +  gamma*Sl + r*I0
  dh= rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl)
  dhr= rho*f*Sl
  dhl= rho*((omega*lambda*(Il+I0))*(S0+Sl) + f*Sl)
  dhh= ((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl)
  dhhl=((omega*lambda*(Il+I0))*(S0+Sl) + f*Sl)
  res=c(dIl, dI0, dSl, dS0, dh, dhr, dhl, dhh, dhhl)
  return(list(res))
}





#' @title Simulate vivax model
#' @description Simulates a draw of the chosen vivax ODE model
#' @param parameters model parameters
#' @param ODEmodel  ODE model to be simulated
#' @param maxtime maximal time step
#' @param year if TRUE, aggregates outputs per year. if FALSE, returns daily output
#'
#' @return A dataframe containing the simulated model
#' @export
#'
simulate_vivax_ode <- function(parameters, ODEmodel=ode_vivax_cm_importation_vc, maxtime=1465, year=FALSE){

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
  #simulation
  state      = c(parameters["Il"][[1]],
                  parameters["I0"][[1]],
                  parameters["Sl"][[1]],
                  parameters["S0"][[1]],
                  parameters["h"][[1]],
                  parameters["hr"][[1]],
                 parameters["hl"][[1]],
                 parameters["hh"][[1]],
                 parameters["hhl"][[1]])

  times = seq(0, maxtime, by = 1)
  solveSIS = deSolve::ode(y = state, times = times, func = ODEmodel, parms =parameters, method = "lsoda")
  solutionVivax=as.data.frame(solveSIS)
  names(solutionVivax)=c("time", "Il", "I0", "Sl", "S0", "h", "hr", "hl", "hh", "hhl")
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



##################################
# Calculate R0
##################################

#' @title R0 calculation
#'
#' @description Calculation of R0 in vivax model
#'
#' @param lambda transmission parameter
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#'
#' @return A scalar R0 value
#' @export
#'
get_r0_vivax <- function(lambda, f, r, gamma){
  return(lambda*(gamma+f)*(gamma+r)/r/gamma/(gamma+f+r))
}





#' @title Rc calculation, with case management and vector control
#'
#' @description Calculation of Rc in vivax model
#'
#' @param lambda transmission parameter
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param omega intensity of vector control. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#'
#' @return A scalar Rc value
#' @export
#'
get_rc_vivax <- function(lambda, f, r, gamma, alpha, beta, omega=1){
  mynum=omega*lambda*(1-alpha)*(gamma+r)*(f+gamma)
  mydenom=r*(f*alpha*beta*r+gamma*f*(1-alpha*(1-beta))+gamma*(gamma+r))
  return(mynum/mydenom)
}





##################################
# Back-calculate lambda
##################################

#' @title Transmission rate calculation
#'
#' @description Solve the equation for lambda, vivax with CM, vector control and importation
#'
#' @param h daily incidence
#' @param p proportion of imported cases
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param rho observation rate
#' @param omega vector control intensity. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#' @param return.all if TRUE, returns lambda, I and delta. If FALSE, returns only lambda
#'
#' @return If return.all=TRUE, returns a list with lambda, I and delta. If return.all=FALSE, returns only a lambda value
#'
#' @export
#'
solve_lambda_vivax <- function(h, r,  f, gamma, alpha, beta, rho, p, omega, return.all=FALSE){

  a_1 <- function(I, r,  f, gamma, alpha, beta, rho, delta){
    return(I^3*(1-I))
  }

  a_2 <- function(I, r,  f, gamma, alpha, beta, rho, delta){
    return((I^2)*(1-I)*(2*gamma+ 3*delta + r + f) - r*(I^3)/(1-alpha))
  }

  a_3 <- function(I, r,  f, gamma, alpha, beta, rho, delta){
    return(I*(1-I)*((gamma +delta + r)*(gamma +delta + f)+delta*(2*gamma+r+f+2*delta) )
           - r*(I^2)*(2*gamma +r+f*alpha*beta +2*delta)/(1-alpha))
  }

  a_4 <- function(I, r,  f, gamma, alpha, beta, rho, delta){
    return((gamma+delta+r)*(delta*(1-I)*(delta+f+gamma)-r*I*(gamma+delta+(1-alpha*(1-beta))*f)/(1-alpha))+f*r*(delta+r)*I)
  }

  I=h*(1-alpha)/r/rho
  delta=p * h/(1-I)/rho
  # solving the equation for lambda
  lambda_complex =polyroot(c(a_4(I, r,  f, gamma, alpha, beta, rho, delta),
                             a_3(I, r,  f, gamma, alpha, beta, rho, delta),
                             a_2(I, r,  f, gamma, alpha, beta, rho, delta),
                             a_1(I, r,  f, gamma, alpha, beta, rho, delta)))



  lambda_real=as.double(lambda_complex[abs(Im(lambda_complex)) < 1e-10])
  #print(lambda_complex)
  lambda_real=lambda_real[lambda_real>=0]

  if(return.all){
    return(list("lambda_real"=lambda_real/omega, "I"=I, "delta"=delta))
  } else{
    return(lambda_real/omega)
  }

}



#' @title Calculate proportion of relapses
#'
#' @description Calculates the proportion of new infections due to relapse at equilibrium
#'
#' @param h daily incidence
#' @param p proportion of imported cases
#' @param lambda transmission rate
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param rho observation rate
#' @param omega vector control intensity. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#'
#' @return A scalar giving the proportion of relapses
#' @export
#'
get_prop_relapse <- function(lambda, h, r,  f, gamma, alpha, beta, rho, p, omega){
  lambda=omega*lambda
  I=h*(1-alpha)/r/rho
  delta=p * h/(1-I)/rho
  Sl=(r*(I*(lambda*I+delta+r)/(lambda*I+delta+gamma+r))+alpha*(1-beta)*(lambda*I+delta)*(1-I))/(lambda*I+delta+(1-alpha*(1-beta))*f+gamma)

  return(rho*f*Sl/h)
}


#' @title Calculate equilibrium state variables
#'
#' @description Calculates equilibrium state variables I0, Il, S0 and Sl based on I and model parameters
#'
#' @param I proportion of infectious individuals at equilibrium (I0+Il)
#' @param lambda transmission rate
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param rho observation rate
#' @param omega vector control intensity. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#' @param delta importation rate
#'
#' @return A list with the equilibrium states (I0, Il, S0 and Sl)
#' @export
#'
get_equilibrium_states_vivax <- function(I, lambda, r, gamma, f, alpha, beta, rho, delta, omega){
  lambda=lambda*omega
  Il=I*(lambda*I+delta+r)/(lambda*I+delta+gamma+r)
  I0=I-Il
  Sl=(r*Il+alpha*(1-beta)*(lambda*I+delta)*(1-I))/(lambda*I+delta+(1-alpha*(1-beta))*f+gamma)
  S0=1-I-Sl
  h=rho*((lambda*(I0+Il)+delta)*(S0+Sl)+f*Sl)
  hr=rho*f*Sl
  hl=rho*((lambda*(I0+Il))*(S0+Sl)+f*Sl)
  hh=((lambda*(I0+Il)+delta)*(S0+Sl)+f*Sl)
  hhl=((lambda*(I0+Il))*(S0+Sl)+f*Sl)
  return(list("Il"=Il, "I0"=I0, "Sl"=Sl, "S0"=S0,
              "h"=h, "hr"=hr, "hl"=hl, "hh"=hh, "hhl"=hhl))
}



#' @title Checks validity condition
#'
#' @description Check if validity condition for transmission rate calculation
#' based on equilibrium prevalence is verified.
#'
#' @param I proportion of infectious individuals at equilibrium (I0+Il)
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param delta importation rate
#'
#' @return A scalar giving the condition. If negative, the solution exists. If positive there is no positive solution for lambda.
#' @export
#'
get_validity_condition_vivax <- function(I, r,  f, gamma, alpha, beta, delta){
  return((gamma+delta+r)*(delta*(1-I)*(delta+f+gamma)-r*I*(gamma+delta+(1-alpha*(1-beta))*f)/(1-alpha))+f*r*(delta+r)*I)
}
