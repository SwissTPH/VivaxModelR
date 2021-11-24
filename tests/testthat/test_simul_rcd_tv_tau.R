test_that("test compare RCD incl. time varying tau in non delay model with non RCD model", {


  my_tau=function(pr){
    return(varying_tau(nu=5, pr=pr, N=10000))
  }

  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.00155531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=1,"omega"=1,
                      "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0,
                      "hl"=0, "hh"=0, "hhl"=0,"kappa"=1,
                      "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1)

  parameters_rcd_tv=parameters_rcd
  parameters_rcd_tv$tau=my_tau


  parameters_rcd_0iota=parameters_rcd_tv; parameters_rcd_0iota$iota=0
  parameters_rcd_0eta=parameters_rcd_tv; parameters_rcd_0eta$eta=0

  vivax_0=simulate_vivax_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_cm_importation_vc, maxtime=1465, year=FALSE)
  rcd_0iota=simulate_vivax_ode(parameters=parameters_rcd_0iota , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  rcd_0eta=simulate_vivax_ode(parameters=parameters_rcd_0eta , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  vivax_rcd=simulate_vivax_ode(parameters=parameters_rcd_tv , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  vivax_rcd_f=simulate_vivax_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)

  expect_equal(vivax_0[1465,], rcd_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax_0[1465,], rcd_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_gt(vivax_0$I[1465], vivax_rcd$I[1465], label = "rcd better than no rcd")
  expect_gt(vivax_rcd_f$I[1465], vivax_rcd$I[1465], label = "rcd tv better than rcd fixed")



  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.0055531,"delta"=0.0001,
                       "alpha"=0.2, "beta"=0.7, "rho"=1,"omega"=0.9,
                       "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1,
                       "h"=0, "hr"=0,"hl"=0,  "hh"=0, "hhl"=0,"kappa"=1,
                       "tau"=1.5, "nu"=5, "iota"=5/7/10000, "eta"=1)

  parameters2_rcd_tv=parameters2_rcd
  parameters2_rcd_tv$tau=my_tau

  parameters2_rcd_0iota=parameters2_rcd_tv; parameters2_rcd_0iota$iota=0
  parameters2_rcd_0eta=parameters2_rcd_tv; parameters2_rcd_0eta$eta=0

  vivax2_0=simulate_vivax_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_cm_importation_vc, maxtime=1465, year=FALSE)
  rcd2_0iota=simulate_vivax_ode(parameters=parameters2_rcd_0iota , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  rcd2_0eta=simulate_vivax_ode(parameters=parameters2_rcd_0eta , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  vivax2_rcd=simulate_vivax_ode(parameters=parameters2_rcd_tv , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  vivax2_rcd_f=simulate_vivax_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)

  expect_equal(vivax2_0[1465,], rcd2_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax2_0[1465,], rcd2_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_gt(vivax2_0$I[1465], vivax2_rcd$I[1465], label = "rcd")
  expect_gt(vivax2_rcd_f$I[1465], vivax2_rcd$I[1465], label = "rcd")


})

test_that("test compare RCD incl. time varying tau  in delay model (non referral) with non RCD model", {

  my_tau=function(pr){
    return(varying_tau(nu=5, pr=pr, N=10000))
  }

  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.00155531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=1,"omega"=1,
                      "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0,
                      "sigma"=1/15, "T0"=0, "Tl"=0,"kappa"=1,
                      "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1)


  parameters_rcd_tv=parameters_rcd
  parameters_rcd_tv$tau=my_tau

  parameters_rcd_0iota=parameters_rcd_tv; parameters_rcd_0iota$iota=0
  parameters_rcd_0eta=parameters_rcd_tv; parameters_rcd_0eta$eta=0
  parameters_rcd_0sigma=parameters_rcd_tv; parameters_rcd_0sigma$sigma=0

  vivax_0=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay, maxtime=1465, year=FALSE)
  rcd_0iota=simulate_vivax_delay_ode(parameters=parameters_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  rcd_0eta=simulate_vivax_delay_ode(parameters=parameters_rcd_0eta , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  rcd_0sigma=simulate_vivax_delay_ode(parameters=parameters_rcd_0sigma , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  vivax_rcd=simulate_vivax_delay_ode(parameters=parameters_rcd_tv , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  vivax_rcd_f=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)

  expect_equal(vivax_0[1465,], rcd_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax_0[1465,], rcd_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_equal(vivax_rcd, rcd_0sigma, tolerance = 1e-09, label = "sigma=0 and normal rcd")
  expect_gt(vivax_0$I[1465], vivax_rcd$I[1465], label = "rcd")
  expect_gt(vivax_rcd_f$I[1465], vivax_rcd$I[1465], label = "rcd")


  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.00155531,"delta"=0.0001,
                       "alpha"=0.2, "beta"=0.7, "rho"=1,"omega"=0.9,
                       "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0,"hl"=0,"hh"=0, "hhl"=0,
                       "sigma"=1/15, "T0"=0, "Tl"=0,"kappa"=1,
                       "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1)

  parameters2_rcd_tv=parameters2_rcd
  parameters2_rcd_tv$tau=my_tau

  parameters2_rcd_0iota=parameters2_rcd_tv; parameters2_rcd_0iota$iota=0
  parameters2_rcd_0eta=parameters2_rcd_tv; parameters2_rcd_0eta$eta=0
  parameters2_rcd_0sigma=parameters2_rcd_tv; parameters2_rcd_0sigma$sigma=0
  parameters2_rcd_infsigma=parameters2_rcd_tv; parameters2_rcd_infsigma$sigma=1/0.4

  vivax2_0=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay, maxtime=1465, year=FALSE)
  rcd2_0iota=simulate_vivax_delay_ode(parameters=parameters2_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  rcd2_0eta=simulate_vivax_delay_ode(parameters=parameters2_rcd_0eta , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  vivax2_rcd=simulate_vivax_delay_ode(parameters=parameters2_rcd_tv , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  vivax2_rcd_f=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)

  expect_equal(vivax2_0[1465,], rcd2_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax2_0[1465,], rcd2_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_gt(vivax2_0$I[1465], vivax2_rcd$I[1465], label = "rcd")
  expect_gt(vivax2_rcd_f$I[1465], vivax2_rcd$I[1465], label = "rcd")

})


test_that("test compare RCD in delay model (referral) with non RCD model", {

  my_tau=function(pr){
    return(varying_tau(nu=5, pr=pr, N=10000))
  }

  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.00455531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=1,"omega"=1,
                      "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0,"hl"=0,"hh"=0, "hhl"=0,
                      "sigma"=1/15, "T0"=0, "Tl"=0,"kappa"=1,
                      "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1)

  parameters_rcd_tv=parameters_rcd
  parameters_rcd_tv$tau=my_tau

  parameters_rcd_0iota=parameters_rcd_tv; parameters_rcd_0iota$iota=0
  parameters_rcd_0eta=parameters_rcd_tv; parameters_rcd_0eta$eta=0
  parameters_rcd_0sigma=parameters_rcd_tv; parameters_rcd_0sigma$sigma=0

  vivax_0=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay, maxtime=1465, year=FALSE)
  rcd_0iota=simulate_vivax_delay_ode(parameters=parameters_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  rcd_0eta=simulate_vivax_delay_ode(parameters=parameters_rcd_0eta , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  rcd_0sigma=simulate_vivax_delay_ode(parameters=parameters_rcd_0sigma , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  vivax_rcd=simulate_vivax_delay_ode(parameters=parameters_rcd_tv , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  vivax_rcd_f=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)

  expect_equal(vivax_0[1465,], rcd_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax_0[1465,], rcd_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_equal(vivax_0[1465,c("Sl", "S0", "I")], rcd_0sigma[1465,c("Sl", "S0", "I")], tolerance = 1e-09, label = "sigma=0")
  expect_gt(vivax_0$I[1465], vivax_rcd$I[1465], label = "rcd")
  expect_lt(vivax_rcd$I[1465], rcd_0sigma$I[1465],  label = "delay increases eq. prevalence")
  expect_lt(vivax_rcd$I[1465], vivax_rcd_f$I[1465],  label = "fixed increases eq. prevalence")


  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.00155531,"delta"=0.001,
                       "alpha"=0.2, "beta"=0.7, "rho"=1,"omega"=0.9,
                       "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0,"hl"=0,"hh"=0, "hhl"=0,
                       "sigma"=1/15, "T0"=0, "Tl"=0,"kappa"=1,
                       "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1)

  parameters2_rcd_tv=parameters2_rcd
  parameters2_rcd_tv$tau=my_tau

  parameters2_rcd_0iota=parameters2_rcd_tv; parameters2_rcd_0iota$iota=0
  parameters2_rcd_0eta=parameters2_rcd_tv; parameters2_rcd_0eta$eta=0
  parameters2_rcd_0sigma=parameters2_rcd_tv; parameters2_rcd_0sigma$sigma=0
  parameters2_rcd_infsigma=parameters2_rcd_tv; parameters2_rcd_infsigma$sigma=1/0.4

  vivax2_0=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay, maxtime=1465, year=FALSE)
  rcd2_0iota=simulate_vivax_delay_ode(parameters=parameters2_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  rcd2_0eta=simulate_vivax_delay_ode(parameters=parameters2_rcd_0eta , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  rcd2_0sigma=simulate_vivax_delay_ode(parameters=parameters2_rcd_0sigma , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  vivax2_rcd=simulate_vivax_delay_ode(parameters=parameters2_rcd_tv , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  vivax2_rcd_f=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)

  expect_equal(vivax2_0[1465,], rcd2_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax2_0[1465,], rcd2_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_gt(vivax2_0$I[1465], vivax2_rcd$I[1465], label = "rcd")
  expect_lt(vivax2_rcd$I[1465], rcd2_0sigma$I[1465],  label = "delay increases eq. prevalence")
  expect_lt(vivax2_rcd$I[1465], vivax2_rcd_f$I[1465],  label = "fixed tau increases eq. prevalence")

})

test_that("test compare the 3 RCD models", {

  my_tau=function(pr){
    return(varying_tau(nu=5, pr=pr, N=10000))
  }

  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.0155531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=1,"omega"=1,
                      "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0,
                      "sigma"=1/15, "T0"=0, "Tl"=0,"kappa"=1,
                      "tau"=my_tau, "nu"=5, "iota"=5/7/10000, "eta"=1)


  parameters_rcd_infsigma=parameters_rcd; parameters_rcd_infsigma$sigma=1/0.4

  vivax_rcd=simulate_vivax_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  vivax_delay_non_ref=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  vivax_delay_ref=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  vivax_delay_ref_siginf=simulate_vivax_delay_ode(parameters=parameters_rcd_infsigma , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)

  expect_equal(vivax_delay_non_ref$I[1465], vivax_rcd$I[1465], tolerance = 1e-09, label = "rcd non ref = rcd no delay")
  expect_lt(vivax_delay_non_ref$I[1465], vivax_delay_ref$I[1465], label = "rcd non ref is more effective")
  expect_equal(vivax_delay_non_ref$I[1465], vivax_delay_ref_siginf$I[1465], tolerance = 1e-03, label = "rcd non ref approx rcd ref if sigma =inf")


  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.0155531,"delta"=0,
                       "alpha"=0.2, "beta"=0.7, "rho"=1,"omega"=0.9,
                       "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0,
                       "sigma"=1/15, "T0"=0, "Tl"=0,"kappa"=1,
                       "tau"=my_tau, "nu"=5, "iota"=5/7/10000, "eta"=1)
  parameters2_rcd_infsigma=parameters2_rcd; parameters2_rcd_infsigma$sigma=1/0.4
  parameters2_rcd_0sigma=parameters2_rcd; parameters2_rcd_0sigma$sigma=0

  vivax2_rcd=simulate_vivax_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_rcd, maxtime=5000, year=FALSE)
  vivax2_delay_non_ref=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  vivax2_delay_ref=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)
  vivax2_delay_non_ref_sig0=simulate_vivax_delay_ode(parameters=parameters2_rcd_0sigma , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  vivax2_delay_non_ref_siginf=simulate_vivax_delay_ode(parameters=parameters2_rcd_infsigma , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  vivax2_delay_ref_siginf=simulate_vivax_delay_ode(parameters=parameters2_rcd_infsigma , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)

  expect_gt(vivax2_delay_non_ref$I[5000], vivax2_rcd$I[5000], label = "no delay is more effective")
  expect_lt(vivax2_delay_non_ref$I[5000], vivax2_delay_ref$I[5000], label = "rcd non ref is more effective")

})



test_that("test simulation of future scenarios, with RCD and time varying tau", {

  mydata=data.frame(incidence=c(23,112),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$prop_import=c(0,0.01)
  mydata$omega=c(1,1)
  f=1/72
  gamma=1/223
  r=1/60

  mydata2=calculate_r0_rc_fromdata(mydata,f=f, gamma=gamma, r=r, return.all = T )

  my_tau=function(pr){
    return(varying_tau(nu=5, pr=pr, N=10000))
  }

  int_A=list(intervention_name="A", "alpha.new"=0.3, "beta.new"=0.6, "omega.new"=0.9, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  int_B=list(intervention_name="B", "alpha.new"=0.3, "beta.new"=0.6, "omega.new"=0.9, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=my_tau, "rho.new"=0.3, "kappa.new"=0.3)
  my_intervention_list=list(int_A,int_B)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*3, rcd=T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*3, rcd=T)

  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,f=f, gamma=gamma, r=r, year=F, maxtime = 365*6, rcd=T)
  row.names(simul2)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2 %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step) ,
               tolerance = 1e-07, label = "chaining from equilibrium is not the same as simulating from equilibrium")

  expect_gt(simul_1_2[simul_1_2$id==1 & simul_1_2$intervention=="A",]$h[1000], simul_1_2[simul_1_2$id==1 & simul_1_2$intervention=="B",]$h[1000], label = "with importation: h>hl")
  expect_gt(simul_1_2[simul_1_2$id==2 & simul_1_2$intervention=="A",]$h[1000], simul_1_2[simul_1_2$id==2 & simul_1_2$intervention=="B",]$h[1000], label = "with importation: h>hl")

})




test_that("test simulation of future scenarios, with MDA and RCD incl. time varying tau", {

  mydata=data.frame(incidence=c(273,312),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$prop_import=c(0.1,0.001)
  mydata$omega=c(1,1)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata(mydata,f=f, gamma=gamma, r=r, return.all = T )

  my_tau=function(pr){
    return(varying_tau(nu=5, pr=pr, N=10000))
  }

  int_A=list(intervention_name="A","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "kappa.new"=0.18)
  int_B=list(intervention_name="B","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=my_tau, "kappa.new"=0.18)
  int_AM=list(intervention_name="A+MDA","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "kappa.new"=0.18)
  int_BM=list(intervention_name="B+MDA","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=my_tau, "kappa.new"=0.18)
  my_intervention_list=list(int_A, int_B,int_AM,int_BM)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*1,mda = F, rcd=T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*10, mda = T, rcd=T)

  row.names(simul2)=NULL
  expect_equal(simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A", "B")) %>% dplyr::select(-intervention),
               simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A+MDA", "B+MDA")) %>% dplyr::select(-intervention),
               tolerance = 5e-05, label = "MDA effect disappears in the long term")

  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="B+MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="B" & simul2$id==1],
            label =  "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="B+MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="B" & simul2$id==2],
            label =  "MDA better than no MDA")


  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="B+MDA" & simul2$id==1],
            label =  "A better than B for high PR")
  expect_lt(varying_tau(int_B$nu.new, simul2$I[simul2$time==4000 & simul2$intervention =="B+MDA" & simul2$id==1], 10000),
            int_A$tau.new,
            label =  "check the tau values")

  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="B+MDA" & simul2$id==2],
            label =  "A better than B for high PR")
  expect_lt(varying_tau(int_B$nu.new, simul2$I[simul2$time==4000 & simul2$intervention =="B+MDA" & simul2$id==2], 10000),
            int_A$tau.new,
            label =  "check the tau values")



})




test_that("test simulation of future scenarios, with RCD incl. time varying tau and delay", {

  mydata=data.frame(incidence=c(323,312),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0.01,0.01)
  mydata$omega=c(1,1)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata_delay(mydata,f=f, gamma=gamma, r=r, return.all = T )

  my_tau=function(pr){
    return(varying_tau(nu=5, pr=pr, N=10000))
  }


  int_A=list(intervention_name="A","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA,"MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "kappa.new"=0.18)
  int_B=list(intervention_name="B","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA,"MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=my_tau, "kappa.new"=0.18)
  int_AM=list(intervention_name="A+MDA","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA,"sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "kappa.new"=0.18)
  int_BM=list(intervention_name="B+MDA","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA,"sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=my_tau, "kappa.new"=0.18)
  my_intervention_list=list(int_A, int_B,int_AM,int_BM)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*1,mda = F, rcd=T, delay = T, referral = F)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*10, mda = T, rcd=T, delay = T, referral = F)

  simul3=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*1,mda = F, rcd=T, delay = T, referral =T)

  simul4=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*10, mda = T, rcd=T, delay = T, referral = T)


  row.names(simul2)=NULL
  expect_equal(simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A", "B")) %>% dplyr::select(-intervention),
               simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A+MDA", "B+MDA")) %>% dplyr::select(-intervention),
               tolerance = 5e-05, label = "MDA effect disappears in the long term")

  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="B+MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="B" & simul2$id==1],
            label =  "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="B+MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="B" & simul2$id==2],
            label =  "MDA better than no MDA")


  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="B+MDA" & simul2$id==1],
            label =  "A better than B for high PR")
  expect_lt(varying_tau(int_B$nu.new, simul2$I[simul2$time==4000 & simul2$intervention =="B+MDA" & simul2$id==1], 10000),
            int_A$tau.new,
            label =  "check the tau values")

  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="B+MDA" & simul2$id==2],
            label =  "A better than B for high PR")
  expect_lt(varying_tau(int_B$nu.new, simul2$I[simul2$time==4000 & simul2$intervention =="B+MDA" & simul2$id==2], 10000),
            int_A$tau.new,
            label =  "check the tau values")


  row.names(simul4)=NULL
  expect_equal(simul4 %>% dplyr::filter(time==11*365, intervention %in% c("A", "B")) %>% dplyr::select(-intervention),
               simul4 %>% dplyr::filter(time==11*365, intervention %in% c("A+MDA", "B+MDA")) %>% dplyr::select(-intervention),
               tolerance = 5e-05, label = "MDA effect disappears in the long term")

  expect_lt(simul4$I[simul4$time==366 & simul4$intervention =="B+MDA" & simul4$id==1],
            simul4$I[simul4$time==366 & simul4$intervention =="B" & simul4$id==1],
            label =  "MDA better than no MDA")
  expect_lt(simul4$I[simul4$time==366 & simul4$intervention =="B+MDA" & simul4$id==2],
            simul4$I[simul4$time==366 & simul4$intervention =="B" & simul4$id==2],
            label =  "MDA better than no MDA")


  expect_lt(simul4$I[simul4$time==366 & simul4$intervention =="A+MDA" & simul4$id==1],
            simul4$I[simul4$time==366 & simul4$intervention =="B+MDA" & simul4$id==1],
            label =  "A better than B for high PR")
  expect_lt(varying_tau(int_B$nu.new, simul4$I[simul4$time==4000 & simul4$intervention =="B+MDA" & simul4$id==1], 10000),
            int_A$tau.new,
            label =  "check the tau values")

  expect_lt(simul4$I[simul4$time==366 & simul4$intervention =="A+MDA" & simul4$id==2],
            simul4$I[simul4$time==366 & simul4$intervention =="B+MDA" & simul4$id==2],
            label =  "A better than B for high PR")
  expect_lt(varying_tau(int_B$nu.new, simul4$I[simul4$time==4000 & simul4$intervention =="B+MDA" & simul4$id==2], 10000),
            int_A$tau.new,
            label =  "check the tau values")
})
