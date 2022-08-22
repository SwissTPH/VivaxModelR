#' @title ODE vivax model with case management, RCD, importation and vector control, including delay to treatment
#' @description Defines the equations of the model
#' @param t time
#' @param y state variables
#' @param parameters model parameters
#'
#' @return a list of the state variables for all time steps
#' @export

ode_vivax_rcd <- function(t, y, parameters) {
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
  h=y[5]
  hr=y[6]
  hl=y[7]
  hh=y[8]
  hhl=y[9]

  dIl= (1-alpha)*(omega*lambda*(Il+I0)+delta)*(S0+Sl) + (omega*lambda*(Il+I0)+delta)*I0 + (1-alpha)*f*Sl - gamma*Il - r*Il -Il*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl) )
  dI0= -(omega*lambda*(Il+I0)+delta)*(I0) +  gamma*Il - r*I0 -I0*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl) )
  dSl= -(1-alpha*(1-beta))*(omega*lambda*(Il+I0)+delta)*(Sl) - (1-alpha*(1-beta))*f*Sl +(omega*lambda*(Il+I0)+delta)*alpha*(1-beta)*S0 -  gamma*Sl + r*Il +(1-beta)*Il*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl) )
  dS0= -(1-alpha*beta)*(omega*lambda*(Il+I0)+delta)*(S0) + (omega*lambda*(Il+I0)+delta)*alpha*beta*Sl +alpha*beta*f*Sl +  gamma*Sl + r*I0 + (beta*Il+I0)*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl) )
  dh= rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl) + (Il+I0)*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl))*corr_factor_rcd
  dhr= rho*f*Sl
  dhl= rho*((omega*lambda*(Il+I0))*(S0+Sl) + f*Sl)+ (Il+I0)*nu*tau(pr=(I0+Il))*eta*min(iota,rho*((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl))*corr_factor_rcd
  dhh= ((omega*lambda*(Il+I0)+delta)*(S0+Sl) + f*Sl)
  dhhl= ((omega*lambda*(Il+I0))*(S0+Sl) + f*Sl)
  res=c(dIl, dI0, dSl, dS0, dh, dhr, dhl, dhh, dhhl)
  return(list(res))
}



#' @title ODE vivax model with case management, RCD with referral, importation and vector control, including delay to treatment
#' @description Defines the equations of the model
#' @param t time
#' @param y state variables
#' @param parameters model parameters
#'
#' @return a list of the state variables for all time steps
#' @export

ode_vivax_delay_rcd_referral <- function(t, y, parameters) {
  lambda=parameters["lambda"][[1]]
  f=parameters["f"][[1]]
  gamma=parameters["gamma"][[1]]
  r=parameters["r"][[1]]
  alpha=parameters["alpha"][[1]]
  beta=parameters["beta"][[1]]
  rho=parameters["rho"][[1]]
  kappa=parameters["kappa"][[1]]
  sigma=parameters["sigma"][[1]]
  omega=parameters["omega"][[1]](t)
  delta=parameters["delta"][[1]](t)
  iota=parameters["iota"][[1]]
  nu=parameters["nu"][[1]]
  eta=parameters["eta"][[1]]
  tau0=parameters["tau"][[1]]

  if(is.numeric(tau0)){
    tau=function(pr){return(tau0)}
  } else{ tau=tau0}

  corr_factor_rcd=(1-alpha-rho+rho*alpha/kappa)/(1-alpha)

  Ul=y[1]
  U0=y[2]
  Sl=y[3]
  S0=y[4]
  Tl=y[5]
  T0=y[6]
  h=y[7]
  hl=y[8]
  hhh=y[9]
  hhl=y[10]

  dUl= (1-alpha)*(omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) + (1-alpha)*f*Sl + (omega*lambda*(Ul+U0+Tl+T0)+delta)*U0 - gamma*Ul - r*Ul -Ul*nu*tau(pr=(Ul+U0+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dU0= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(U0) +  gamma*Ul - r*U0 -U0*nu*tau(pr=(Ul+U0+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dSl= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(Sl) - f*Sl +(1-beta)*sigma*Tl -  gamma*Sl + r*Ul+ r*Tl
  dS0= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0) +beta*sigma*Tl+ sigma *T0 +  gamma*Sl + r*U0+ r*T0
  dTl= alpha*(omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) + alpha*f*Sl - sigma*Tl -r*Tl -gamma*Tl+ (omega*lambda*(Ul+U0+Tl+T0)+delta)*T0 +Ul*nu*tau(pr=(Ul+U0+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dT0= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(T0) - sigma*T0+gamma*Tl - r*T0 +U0*nu*tau(pr=(Ul+U0+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dh= rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl) + (Ul+U0)*nu*tau(pr=(Ul+U0+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))*corr_factor_rcd
  dhl= rho*((omega*lambda*(Ul+U0+Tl+T0))*(S0+Sl) +f*Sl)+ (Ul+U0)*nu*tau(pr=(Ul+U0+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))*corr_factor_rcd
  dhh= ((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl)
  dhhl= ((omega*lambda*(Ul+U0+Tl+T0))*(S0+Sl) +f*Sl)
  res=c(dUl, dU0, dSl, dS0, dTl, dT0, dh, dhl, dhh, dhhl)
  return(list(res))
}



#' @title ODE vivax model with case management, RCD without referral, importation and vector control, including delay to treatment
#' @description Defines the equations of the model
#' @param t time
#' @param y state variables
#' @param parameters model parameters
#'
#' @return a list of the state variables for all time steps
#' @export

ode_vivax_delay_rcd_no_referral <- function(t, y, parameters) {
  lambda=parameters["lambda"][[1]]
  f=parameters["f"][[1]]
  gamma=parameters["gamma"][[1]]
  r=parameters["r"][[1]]
  alpha=parameters["alpha"][[1]]
  beta=parameters["beta"][[1]]
  rho=parameters["rho"][[1]]
  kappa=parameters["kappa"][[1]]
  sigma=parameters["sigma"][[1]]
  omega=parameters["omega"][[1]](t)
  delta=parameters["delta"][[1]](t)
  iota=parameters["iota"][[1]]
  nu=parameters["nu"][[1]]
  eta=parameters["eta"][[1]]
  tau0=parameters["tau"][[1]]

  if(is.numeric(tau0)){
    tau=function(pr){return(tau0)}
  } else{ tau=tau0}

  corr_factor_rcd=(1-alpha-rho+rho*alpha/kappa)/(1-alpha)

  Ul=y[1]
  U0=y[2]
  Sl=y[3]
  S0=y[4]
  Tl=y[5]
  T0=y[6]
  h=y[7]
  hl=y[8]
  hhh=y[9]
  hhl=y[10]

  dUl= (1-alpha)*(omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) + (1-alpha)*f*Sl + (omega*lambda*(Ul+U0+Tl+T0)+delta)*U0 - gamma*Ul - r*Ul -Ul*nu*tau(pr=(U0+Ul+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dU0= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(U0) +  gamma*Ul - r*U0 -U0*nu*tau(pr=(U0+Ul+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dSl= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(Sl) - f*Sl +(1-beta)*sigma*Tl -  gamma*Sl + r*Ul+ r*Tl +(1-beta)*Ul*nu*tau(pr=(U0+Ul+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dS0= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0) +beta*sigma*Tl+ sigma *T0 +  gamma*Sl + r*U0+ r*T0 +(beta*Ul+U0)*nu*tau(pr=(U0+Ul+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))
  dTl= alpha*(omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) + alpha*f*Sl - sigma*Tl -r*Tl -gamma*Tl+ (omega*lambda*(Ul+U0+Tl+T0)+delta)*T0
  dT0= -(omega*lambda*(Ul+U0+Tl+T0)+delta)*(T0) - sigma*T0+gamma*Tl - r*T0
  dh= rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl) + (Ul+U0)*nu*tau(pr=(U0+Ul+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))*corr_factor_rcd
  dhl= rho*((omega*lambda*(Ul+U0+Tl+T0))*(S0+Sl) +f*Sl)+ (Ul+U0)*nu*tau(pr=(U0+Ul+Tl+T0))*eta*min(iota,rho*((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl))*corr_factor_rcd
  dhh= ((omega*lambda*(Ul+U0+Tl+T0)+delta)*(S0+Sl) +f*Sl)
  dhhl= ((omega*lambda*(Ul+U0+Tl+T0))*(S0+Sl) +f*Sl)
  res=c(dUl, dU0, dSl, dS0, dTl, dT0, dh, dhl, dhh, dhhl)
  return(list(res))
}




##################################
# Back-calculate lambda
##################################

#' @title Calculate h1 from aggregated incidence data
#'
#' @description Calculate h1 in RCD model for back-calculation when
#'
#' @param h daily incidence
#' @param alpha effective treatment probability
#' @param rho observation rate
#' @param kappa common part in rho and alpha
#' @param iota RCD parameter: maximum number of index cases investigated per population per unit of time
#' @param nu RCD parameter: number of secondary cases investigated per index case
#' @param tau RCD parameter: targeting ratio of the RCD
#' @param eta RCD parameter: probability that an investigated cases is detected (symptoms, test sensitivity, etc.)
#' @param r blood clearance rate
#'
#' @return A scalar equal to the equilibrium value for h1
#'
#' @export
#'
calculate_h1_rcd=function(h, alpha, rho, kappa, iota, nu, tau, eta, r ){
  corr_factor_rcd_sub=1-alpha+rho*alpha/kappa

  rcd_sub=nu*tau*eta
  h1_capped=(rho*rcd_sub*h-rho*r+sqrt((rho*rcd_sub*h-rho*r)^2+4*rho*r*h*corr_factor_rcd_sub*rcd_sub))/(2*rcd_sub*corr_factor_rcd_sub)
  h1_fixed=h*rho*(r+iota*rcd_sub)/(rho*r+iota*rcd_sub*corr_factor_rcd_sub)
  h1=ifelse(h1_capped<=iota & rcd_sub>0, h1_capped, h1_fixed)

  return(h1)
}


#' @title Calculate the targeting ratio from incidence data with information on RCD
#'
#' @description Calculate tau in RCD model when h1 and h2 are available
#'
#' @param h1 daily incidence, directly detected cases
#' @param h2 daily incidence, reactively detected cases
#' @param alpha effective treatment probability
#' @param rho observation rate
#' @param kappa common part in rho and alpha
#' @param iota RCD parameter: maximum number of index cases investigated per population per unit of time
#' @param nu RCD parameter: number of secondary cases investigated per index case
#' @param eta RCD parameter: probability that an investigated cases is detected (symptoms, test sensitivity, etc.)
#' @param r blood clearance rate
#'
#' @return a scalar equal to the targeting ratio tau.
#'
#' @export
#'
calculate_tau_rcd=function(h1, h2, alpha, rho, kappa, iota, nu, eta, r ){

  if(h2< 1e-09){stop("h2=0, we cannot compute tau")}
  if(iota==0 | nu==0 | eta==0){stop("Some RCD parameters =0, we cannot compute tau")}
  iota_star=min(iota, h1)
  corr_factor_rcd_num=1-alpha-rho+rho*alpha/kappa

  I = ((1-alpha)/r )*(h1/rho - h2/corr_factor_rcd_num)

  tau = h2*(1-alpha)/(iota_star*nu*eta*I*corr_factor_rcd_num)

  return(tau)
}


############################################

#' @title Transmission rate calculation
#'
#' @description Solve the equation for lambda, vivax with CM, vector control, RCD and importation
#'
#' @param h daily incidence
#' @param h1 daily incidence, directly detected cases
#' @param p proportion of imported cases
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param rho observation rate
#' @param omega vector control intensity. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#' @param iota RCD parameter: maximum number of index cases investigated per population per unit of time
#' @param nu RCD parameter: number of secondary cases investigated per index case
#' @param tau RCD parameter: targeting ratio of the RCD
#' @param eta RCD parameter: probability that an investigated cases is detected (symptoms, test sensitivity, etc.)
#' @param return.all if TRUE, returns lambda, I and delta. If FALSE, returns only lambda
#'
#' @return If return.all=TRUE, returns a list with lambda, I and delta. If return.all=FALSE, returns only a lambda value
#'
#' @export
#'
solve_lambda_vivax_rcd <- function(h, h1, r,  f, gamma, alpha, beta, rho, p, omega, iota, nu, tau, eta, return.all=FALSE){


  a_1 <- function(I, r,  f, gamma, alpha, beta, rho, delta){
    return(I^3*(1-I))
  }

  a_2 <- function(I, r,  f, gamma, alpha, beta, rho, delta){
    return((I^2)*(1-I)*(2*gamma+ 3*delta + r + f) - r*(I^3)/(1-alpha))
  }

  a_3 <- function(I, r,  f, gamma, alpha, beta, rho, delta, rcd_term){
    return(I*(1-I)*((gamma +delta + r)*(gamma +delta + f)+delta*(2*gamma+r+f+2*delta) )
           - r*(I^2)*(2*gamma +r+f*alpha*beta +2*delta)/(1-alpha)
           - beta*f*(I^2)*rcd_term)
  }

  a_4 <- function(I, r,  f, gamma, alpha, beta, rho, delta, rcd_term){
    return((gamma+delta+r)*(delta*(1-I)*(delta+f+gamma)-r*I*(gamma+delta+(1-alpha*(1-beta))*f)/(1-alpha))+f*r*(delta+r)*I
           - beta*f*I*rcd_term*(delta+r))
  }

  iota_star=min(iota, h1)
  rcd_term=iota_star*nu*tau*eta
  I=h1*(1-alpha)/(r+rcd_term)/rho
  delta=p * h/(1-I)/rho

  # solving the equation for lambda
  lambda_complex =polyroot(c(a_4(I, r= r+rcd_term,  f, gamma, alpha, beta, rho, delta, rcd_term),
                             a_3(I, r= r+rcd_term,  f, gamma, alpha, beta, rho, delta, rcd_term),
                             a_2(I, r= r+rcd_term,  f, gamma, alpha, beta, rho, delta),
                             a_1(I, r= r+rcd_term,  f, gamma, alpha, beta, rho, delta)))



  lambda_real=as.double(lambda_complex[abs(Im(lambda_complex)) < 1e-10])
  #print(lambda_complex)
  lambda_real=lambda_real[lambda_real>=0]

  if(return.all){
    return(list("lambda_real"=lambda_real/omega, "I"=I, "delta"=delta, "iota_star"=iota_star))
  } else{
    return(lambda_real/omega)
  }

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
#' @param kappa common part in rho and alpha
#' @param omega vector control intensity. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#' @param delta importation rate
#' @param iota_star RCD parameter: maximum number of index cases investigated per population per unit of time at equilibrium
#' @param nu RCD parameter: number of secondary cases investigated per index case
#' @param tau RCD parameter: targeting ratio of the RCD
#' @param eta RCD parameter: probability that an investigated cases is detected (symptoms, test sensitivity, etc.)
#'
#' @return A list with the equilibrium states (I0, Il, S0 and Sl)
#' @export
#'
get_equilibrium_states_vivax_rcd <- function(I, lambda, r, gamma, f, alpha, beta, rho, delta, omega, iota_star, nu, tau, eta,kappa){
  lambda=lambda*omega
  rcd_term=iota_star*nu*tau*eta
  corr_factor_rcd=(1-alpha-rho+rho*alpha/kappa)/(1-alpha)

  Il=I*(lambda*I+delta+r+rcd_term)/(lambda*I+delta+gamma+r+rcd_term)
  I0=I-Il
  Sl=((r+(1-beta)*rcd_term)*Il+alpha*(1-beta)*(lambda*I+delta)*(1-I))/(lambda*I+delta+(1-alpha*(1-beta))*f+gamma)
  S0=1-I-Sl
  h=rho*((lambda*(I0+Il)+delta)*(S0+Sl)+f*Sl)+rcd_term*I*corr_factor_rcd
  hr=rho*f*Sl
  hl=rho*((lambda*(I0+Il))*(S0+Sl)+f*Sl)+rcd_term*I*corr_factor_rcd
  hh=(lambda*(I0+Il)+delta)*(S0+Sl)+f*Sl
  hhl=(lambda*(I0+Il))*(S0+Sl)+f*Sl
  return(list("Il"=Il, "I0"=I0, "Sl"=Sl, "S0"=S0,
              "h"=h, "hr"=hr, "hl"=hl,"hh"=hh,"hhl"=hhl ))
}



#' @title Transmission rate calculation
#'
#' @description Solve the equation for lambda, vivax with delay in treatment, vector control, RCD  (no referral) and importation
#'
#' @param h daily incidence
#' @param h1 daily incidence, directly detected cases
#' @param p proportion of imported cases
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param sigma delay to treatment (1/duration of the delay)
#' @param rho observation rate
#' @param omega vector control intensity. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#' @param iota RCD parameter: maximum number of index cases investigated per population per unit of time
#' @param nu RCD parameter: number of secondary cases investigated per index case
#' @param tau RCD parameter: targeting ratio of the RCD
#' @param eta RCD parameter: probability that an investigated cases is detected (symptoms, test sensitivity, etc.)
#' @param return.all if TRUE, returns lambda, I and delta. If FALSE, returns only lambda
#'
#' @return If return.all=TRUE, returns a list with lambda, I, delta and iota_star. If return.all=FALSE, returns only a lambda value
#'
#' @export
#'
solve_lambda_vivax_rcd_no_referral <- function(h, h1, r,  f, gamma, alpha, beta, sigma, rho, p, omega, iota, nu, tau, eta, return.all=FALSE){


  a_1 <- function(I, r,  f, gamma, alpha, beta, rho,sigma,  delta){
    return(I^4*(1-I))
  }

  a_2 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta){
    return((I^3)*(1-I)*(3*gamma+ 4*delta + 2*r +sigma+ f) - r*(I^4)*(r+sigma)/(r+(1-alpha)*sigma))
  }

  a_3 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term){
    return((I^2)*(1-I)*((gamma +delta + r)*(3*gamma +3*delta +r+ 2*f)+delta*(3*gamma+2*r+f+3*delta)+sigma*(3*delta+2*gamma+r+f) )
           - r*(I^3)*((3*gamma +2*r+sigma +3*delta)*(r+sigma)+alpha*beta*f*sigma)/(r+(1-alpha)*sigma)
           -(I^3)*beta*rcd_term*f)
  }

  a_4 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term){
    return(I*(1-I)*((delta+gamma+f)*(delta+gamma+r)*(delta+gamma+r+sigma)+
                      delta*((delta+gamma+f)*(delta+gamma+r)+(delta+gamma+r)*(delta+gamma+r+sigma)+(delta+gamma+f)*(delta+gamma+r+sigma)))+
             I*I*r*(f*(r+(1-alpha*beta)*sigma)*(2*delta+2*r+sigma+gamma)-(r+sigma)*((delta+gamma+f)*(delta+gamma+r)+(delta+gamma+r)*(delta+gamma+r+sigma)+(delta+gamma+f)*(delta+gamma+r+sigma)))/(r+(1-alpha)*sigma)
           - (I^2)*beta*rcd_term*f*(2*delta+gamma+2*r+sigma))

  }

  a_5 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term){
    return((1-I)*delta*(delta+gamma+f)*(delta+gamma+r)*(delta+gamma+r+sigma)-
             I*r*(delta+gamma+r)*((delta+gamma)*(delta+gamma+r)*(r+sigma)+f*gamma*(r+sigma)+alpha*beta*f*sigma*(delta+r))/(r+(1-alpha)*sigma)-
             I*r*sigma*((delta+gamma)*(delta+gamma+r)*(r+sigma)+(1-alpha)*f*gamma*(r+sigma)+alpha*beta*f*sigma*(delta+gamma+r))/(r+(1-alpha)*sigma)
           - I*beta*rcd_term*f*((delta+r)*(delta+gamma+r+sigma)+alpha*r*sigma*gamma/(r+(1-alpha)*sigma)))
  }

  iota_star=min(iota, h1)
  rcd_term=iota_star*nu*tau*eta
  UU=h1*(1-alpha)/(r+rcd_term)/rho
  TT=UU*alpha*(r+rcd_term)/((1-alpha)*(r+sigma))
  I=UU+TT
  delta=p * h/(1-I)/rho
  # solving the equation for lambda
  lambda_complex =polyroot(c(a_5(I, r= r+rcd_term,  f, gamma, alpha, beta, rho, sigma-rcd_term, delta, rcd_term),
                             a_4(I, r= r+rcd_term,  f, gamma, alpha, beta, rho, sigma-rcd_term, delta, rcd_term),
                             a_3(I, r= r+rcd_term,  f, gamma, alpha, beta, rho, sigma-rcd_term, delta, rcd_term),
                             a_2(I, r= r+rcd_term,  f, gamma, alpha, beta, rho, sigma-rcd_term, delta),
                             a_1(I, r= r+rcd_term,  f, gamma, alpha, beta, rho, sigma-rcd_term, delta)))



  lambda_real=as.double(lambda_complex[abs(Im(lambda_complex)) < 1e-10])
  #print(lambda_complex)
  lambda_real=lambda_real[lambda_real>=0]

  if(return.all){
    return(list("lambda_real"=lambda_real/omega, "I"=I, "delta"=delta, "iota_star"=iota_star))
  } else{
    return(lambda_real/omega)
  }

}



#' @title Calculate equilibrium state variables, RCD model (no referral)
#'
#' @description Calculates equilibrium state variables U0, Ul, S0, Sl, T0 and Tl based on I and model parameters
#'
#' @param I proportion of infectious individuals at equilibrium (U0+Ul+T0+Tl)
#' @param lambda transmission rate
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param sigma delay to treatment (1/duration of the delay)
#' @param rho observation rate
#' @param kappa common part in rho and alpha
#' @param omega vector control intensity. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#' @param delta importation rate
#' @param iota_star RCD parameter: maximum number of index cases investigated per population per unit of time at equilibrium
#' @param nu RCD parameter: number of secondary cases investigated per index case
#' @param tau RCD parameter: targeting ratio of the RCD
#' @param eta RCD parameter: probability that an investigated cases is detected (symptoms, test sensitivity, etc.)
#'
#' @return A list with the equilibrium states (U0, Ul, T0, Tl, S0 and Sl)
#' @export
#'
get_equilibrium_states_vivax_rcd_no_referral <- function(I, lambda, r, gamma, f, alpha, beta, rho, sigma, delta, omega, iota_star, nu, tau, eta, kappa){
  lambda=lambda*omega
  rcd_term=iota_star*nu*tau*eta
  corr_factor_rcd=(1-alpha-rho+rho*alpha/kappa)/(1-alpha)

  TT=I*alpha*(r+rcd_term)/(r+(1-alpha)*sigma+alpha*rcd_term)  # T0+Tl
  UU=I*(1-alpha)*(r+sigma)/(r+(1-alpha)*sigma+alpha*rcd_term) #U0 + Ul
  T0=TT*gamma/(lambda*I+delta+gamma+r+sigma)
  Tl=TT-T0
  U0=UU*gamma/(lambda*I+delta+gamma+r+rcd_term)
  Ul=UU-U0
  Sl=(r*Ul+(1-beta)*rcd_term*Ul+(r+(1-beta)*sigma)*Tl)/(lambda*I+delta+gamma+f)
  S0=1-I-Sl

  h=rho*((lambda*(I)+delta)*(1-I)+f*Sl)+rcd_term*UU*corr_factor_rcd
  hl=rho*((lambda*(I))*(1-I)+f*Sl)+rcd_term*UU*corr_factor_rcd
  hh= ((lambda*I+delta)*(1-I) +f*Sl)
  hhl= (lambda*I*(1-I) +f*Sl)
  return(list("Ul"=Ul, "U0"=U0, "Tl"=Tl, "T0"=T0, "Sl"=Sl, "S0"=S0,
              "h"=h, "hl"=hl, "hh"=hh, "hhl"=hhl))
}


#' @title Transmission rate calculation
#'
#' @description Solve the equation for lambda, vivax with delay in treatment, vector control, RCD  (referral) and importation
#'
#' @param h daily incidence
#' @param h1 daily incidence, directly detected cases
#' @param p proportion of imported cases
#' @param f relapse frequency
#' @param r blood clearance rate
#' @param gamma liver clearance rate
#' @param alpha effective treatment probability
#' @param beta probability of liver clearance
#' @param sigma delay to treatment (1/duration of the delay)
#' @param rho observation rate
#' @param omega vector control intensity. If vector control is not the purpose of the analysis, omega can be fixed to 1.
#' @param iota RCD parameter: maximum number of index cases investigated per population per unit of time
#' @param nu RCD parameter: number of secondary cases investigated per index case
#' @param tau RCD parameter: targeting ratio of the RCD
#' @param eta RCD parameter: probability that an investigated cases is detected (symptoms, test sensitivity, etc.)
#' @param return.all if TRUE, returns lambda, I and delta. If FALSE, returns only lambda
#'
#' @return If return.all=TRUE, returns a list with lambda, I, delta and iota_star. If return.all=FALSE, returns only a lambda value
#'
#' @export
#'
solve_lambda_vivax_rcd_referral <- function(h, h1, r,  f, gamma, alpha, beta, sigma, rho, p, omega, iota, nu, tau, eta, return.all=FALSE){


  a_1 <- function(I, r,  f, gamma, alpha, beta, rho,sigma,  delta){
    return(I^4*(1-I))
  }

  a_2 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term){
    return((I^3)*(1-I)*(3*gamma+ 4*delta + 2*r +sigma+ f + rcd_term) - (I^4)*(r+rcd_term)*(r+sigma)/(r+(1-alpha)*sigma+rcd_term))
  }

  a_3 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term){
    return((I^2)*(1-I)*((delta+gamma+f)*(delta+gamma+r+rcd_term)+(delta+gamma+f)*(delta+gamma+r+sigma)+
                          (delta+gamma+r+rcd_term)*(delta+gamma+r+sigma)+delta*(3*delta+3*gamma+2*r+rcd_term+sigma+f))
           - (I^3)*((r+rcd_term)*(r+sigma)*(3*delta+3*gamma+2*r+rcd_term+sigma+f)-
                      f*(r*r+(1-alpha*beta)*r*sigma+rcd_term*r+(1-beta)*sigma*rcd_term))/(r+(1-alpha)*sigma+rcd_term))
  }

  a_4 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term){
    return(I*(1-I)*((delta+gamma+f)*(delta+gamma+r+rcd_term)*(delta+gamma+r+sigma)+
                      delta*((delta+gamma+f)*(2*delta+2*gamma+2*r+rcd_term+sigma)+(delta+gamma+r+rcd_term)*(delta+gamma+r+sigma)))
           - I*I*(r+rcd_term)*(r+sigma)*((delta+gamma+f)*(2*delta+2*gamma+2*r+rcd_term+sigma)+(delta+gamma+r+rcd_term)*(delta+gamma+r+sigma))/(r+(1-alpha)*sigma+rcd_term)
           + I*I*f*(2*delta+gamma+2*r+rcd_term+sigma)*(r*r+(1-alpha*beta)*r*sigma+r*rcd_term+(1-beta)*sigma*rcd_term)/(r+(1-alpha)*sigma+rcd_term)
    )

  }

  a_5 <- function(I, r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term){
    return((1-I)*delta*(delta+gamma+f)*(delta+gamma+r+rcd_term)*(delta+gamma+r+sigma)
           -  I*(r+rcd_term)*(r+sigma)*(delta+gamma+f)*(delta+gamma+r+rcd_term)*(delta+gamma+r+sigma)/(r+(1-alpha)*sigma+rcd_term)
           +  I*f*(delta+gamma+r+rcd_term)*(delta+r+sigma)*(r*r+(1-alpha*beta)*r*sigma+r*rcd_term+(1-beta)*rcd_term*sigma)/(r+(1-alpha)*sigma+rcd_term)
           -  I*f*gamma*sigma*(1-alpha)*(r+(1-beta)*rcd_term)*(r+sigma)/(r+(1-alpha)*sigma+rcd_term)
    )
  }

  iota_star=min(iota, h1)
  rcd_term=iota_star*nu*tau*eta
  UU=h1*(1-alpha)/(r+rcd_term)/rho
  TT=UU*(alpha*r+rcd_term)/((1-alpha)*(r+sigma))
  I=UU+TT
  delta=p * h/(1-I)/rho
  # solving the equation for lambda
  lambda_complex =polyroot(c(a_5(I, r= r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term),
                             a_4(I, r= r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term),
                             a_3(I, r= r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term),
                             a_2(I, r= r,  f, gamma, alpha, beta, rho, sigma, delta, rcd_term),
                             a_1(I, r= r,  f, gamma, alpha, beta, rho, sigma, delta)))



  lambda_real=as.double(lambda_complex[abs(Im(lambda_complex)) < 1e-10])
  #print(lambda_complex)
  lambda_real=lambda_real[lambda_real>=0]

  if(return.all){
    return(list("lambda_real"=lambda_real/omega, "I"=I, "delta"=delta, "iota_star"=iota_star))
  } else{
    return(lambda_real/omega)
  }

}



#' @title Calculate equilibrium state variables, RCD model (referral)
#'
#' @description Calculates equilibrium state variables I0, Il, S0, Sl, T0 and Tl based on I and model parameters
#'
#' @param I proportion of infectious individuals at equilibrium (I0+Il)
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
#' @param iota_star RCD parameter: maximum number of index cases investigated per population per unit of time at equilibrium
#' @param nu RCD parameter: number of secondary cases investigated per index case
#' @param kappa common part in rho and alpha
#' @param tau RCD parameter: targeting ratio of the RCD
#' @param eta RCD parameter: probability that an investigated cases is detected (symptoms, test sensitivity, etc.)
#'
#' @return A list with the equilibrium states (I0, Il, T0, Tl, S0 and Sl)
#' @export
#'
get_equilibrium_states_vivax_rcd_referral <- function(I, lambda, r, gamma, f, alpha, beta, rho, sigma, delta, omega, iota_star, nu, tau, eta, kappa){
  lambda=lambda*omega
  rcd_term=iota_star*nu*tau*eta
  corr_factor_rcd=(1-alpha-rho+rho*alpha/kappa)/(1-alpha)

  UU=I*(1-alpha)*(r+sigma)/(r+(1-alpha)*sigma+rcd_term) #I0 + Il
  TT=I-UU  # T0+Tl
  U0=UU*gamma/(lambda*I+delta+gamma+r+rcd_term)
  Ul=UU-U0
  T0=(TT*gamma+rcd_term*U0)/(lambda*I+delta+gamma+r+sigma)
  Tl=TT-T0
  Sl=(r*Ul+(r+(1-beta)*sigma)*Tl)/(lambda*I+delta+gamma+f)
  S0=1-I-Sl

  h=rho*((lambda*(I)+delta)*(1-I)+f*Sl)+rcd_term*UU*corr_factor_rcd
  hl=rho*((lambda*(I))*(1-I)+f*Sl)+rcd_term*UU*corr_factor_rcd
  hh= ((lambda*I+delta)*(1-I) +f*Sl)
  hhl= (lambda*I*(1-I) +f*Sl)
  return(list("Ul"=Ul, "U0"=U0, "Tl"=Tl, "T0"=T0, "Sl"=Sl, "S0"=S0,
              "h"=h, "hl"=hl, "hh"=hh, "hhl"=hhl))
}



