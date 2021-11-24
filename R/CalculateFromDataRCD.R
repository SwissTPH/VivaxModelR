#' @title Calculate R0 and Rc on an incidence dataset with RCD
#'
#' @description Uses the compartmental model to calculate R0 and Rc using data on
#' incidence and the proportion of imported cases. RCD is accounted for in the calculation of the transmission parameter lambda
#' (but R0 and Rc are the same as in the model without RCD)

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
#' @export
#'
calculate_r0_rc_fromdata_rcd=function(df, f=1/72, gamma=1/223, r=1/60,
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
  if((!"iota" %in% names(df) | !"tau" %in% names(df) | !"nu" %in% names(df) | !"eta" %in% names(df)| !"kappa" %in% names(df)) ){
    stop("RCD parameters are missing, please add them or use model without RCD")
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
    my.iota=dataTransform$iota[i]
    my.nu=dataTransform$nu[i]
    my.eta=dataTransform$eta[i]
    my.tau=dataTransform$tau[i]
    my.kappa=dataTransform$kappa[i]
    if(is.na(my.h)==FALSE & is.infinite(my.h)==FALSE  &  (my.h)>=h.cutoff){
      my.h1=calculate_h1_rcd(h=my.h, alpha=my.alpha,rho=my.rho,
                          iota=my.iota, nu=my.nu,eta=my.eta, tau=my.tau,kappa=my.kappa,r=r )


      my.lambda.solve=solve_lambda_vivax_rcd(h=my.h, h1=my.h1, r=r,  gamma=gamma, f=f,
                                         alpha=my.alpha, beta=my.beta,
                                         rho=my.rho,
                                         p=my.p, omega=my.omega,
                                         iota=my.iota, nu=my.nu,eta=my.eta, tau=my.tau,
                                         return.all = T)
      my.lambda=my.lambda.solve$lambda_real
      my.lambda=ifelse( (length(my.lambda)==0), -2, my.lambda)

    } else{my.lambda =-3}
    dataTransform$delta[i]=ifelse( my.lambda>0, my.lambda.solve$delta , NA)
    dataTransform$I[i]=ifelse( my.lambda>0, my.lambda.solve$I , NA)
    dataTransform$iota_star[i]=ifelse( my.lambda>0, my.lambda.solve$iota_star , NA)

    dataTransform$lambda[i]=my.lambda
    dataTransform$R0[i]=ifelse(my.lambda>0, get_r0_vivax(my.lambda,  f=f,r=r, gamma=gamma), NA)
    dataTransform$Rc[i]=ifelse(my.lambda>0, get_rc_vivax(my.lambda,  f=f,r=r, gamma=gamma, alpha=my.alpha, beta=my.beta, omega=my.omega), NA)
    utils::setTxtProgressBar(pb, i)
  }

  output_small=dataTransform[,! names(dataTransform) %in% c("delta", "I", "iota_star")]

  if(return.all){
    return(dataTransform)
  } else{
    return(output_small)
  }

}



#' @title Calculate R0 and Rc on an incidence dataset with delays and RCD
#'
#' @description Uses the compartmental model to calculate R0 and Rc using data on
#' incidence and the proportion of imported cases. RCD is accounted for in the calculation of the transmission parameter lambda
#' (but R0 and Rc are the same as in the model without RCD)
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
#' @param referral a boolean indicating if the rcd model includes referral. Default (FALSE) is the model without referral for RCD. This parameter is used only if rcd==TRUE.
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
#' @export
#'
calculate_r0_rc_fromdata_delay_rcd=function(df, f=1/72, gamma=1/223, r=1/60,
                                        return.all=F, h.cutoff=5e-08, referral=FALSE){

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
  if((!"iota" %in% names(df) | !"tau" %in% names(df) | !"nu" %in% names(df) | !"eta" %in% names(df)| !"kappa" %in% names(df)) ){
    stop("RCD parameters are missing, please add them or use model without RCD")
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
    my.iota=dataTransform$iota[i]
    my.nu=dataTransform$nu[i]
    my.eta=dataTransform$eta[i]
    my.tau=dataTransform$tau[i]
    my.kappa=dataTransform$kappa[i]
    if(is.na(my.h)==FALSE & is.infinite(my.h)==FALSE  &  (my.h)>=h.cutoff){
      my.h1=calculate_h1_rcd(h=my.h, alpha=my.alpha,rho=my.rho,
                             iota=my.iota, nu=my.nu,eta=my.eta, tau=my.tau,kappa=my.kappa,r=r )

      if(referral){
        my.lambda.solve=solve_lambda_vivax_rcd_referral(h=my.h, h1=my.h1, r=r,  gamma=gamma, f=f,
                                                        alpha=my.alpha, beta=my.beta,sigma=my.sigma,
                                                        rho=my.rho,
                                                        p=my.p, omega=my.omega,
                                                        iota=my.iota, nu=my.nu,eta=my.eta, tau=my.tau,
                                                        return.all = T)
      } else{
        my.lambda.solve=solve_lambda_vivax_rcd_no_referral(h=my.h, h1=my.h1, r=r,  gamma=gamma, f=f,
                                                        alpha=my.alpha, beta=my.beta,sigma=my.sigma,
                                                        rho=my.rho,
                                                        p=my.p, omega=my.omega,
                                                        iota=my.iota, nu=my.nu,eta=my.eta, tau=my.tau,
                                                        return.all = T)
      }

      my.lambda=my.lambda.solve$lambda_real
      my.lambda=ifelse( (length(my.lambda)==0), -2, my.lambda)

    } else{my.lambda =-3}
    dataTransform$delta[i]=ifelse( my.lambda>0, my.lambda.solve$delta , NA)
    dataTransform$I[i]=ifelse( my.lambda>0, my.lambda.solve$I , NA)
    dataTransform$iota_star[i]=ifelse( my.lambda>0, my.lambda.solve$iota_star , NA)

    dataTransform$lambda[i]=my.lambda
    dataTransform$R0[i]=ifelse(my.lambda>0, get_r0_vivax(my.lambda,  f=f,r=r, gamma=gamma), NA)
    dataTransform$Rc[i]=ifelse(my.lambda>0, get_rc_vivax_delay(my.lambda,  f=f,r=r, gamma=gamma, alpha=my.alpha, beta=my.beta, sigma=my.sigma, omega=my.omega), NA)
    utils::setTxtProgressBar(pb, i)
  }

  output_small=dataTransform[,! names(dataTransform) %in% c("delta", "I", "iota_star")]

  if(return.all){
    return(dataTransform)
  } else{
    return(output_small)
  }

}
