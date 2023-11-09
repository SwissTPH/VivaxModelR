test_that("test compare RCD in non delay model with non RCD model", {

  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.0155531,"delta"=0,
                  "alpha"=0, "beta"=1, "rho"=1,"omega"=1,
                  "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0,
                  "hl"=0, "hh"=0, "hhl"=0,
                  "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters_rcd_0iota=parameters_rcd; parameters_rcd_0iota$iota=0
  parameters_rcd_0nu=parameters_rcd; parameters_rcd_0nu$nu=0
  parameters_rcd_0eta=parameters_rcd; parameters_rcd_0eta$eta=0

  vivax_0=simulate_vivax_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_cm_importation_vc, maxtime=1465, year=FALSE)
  rcd_0iota=simulate_vivax_ode(parameters=parameters_rcd_0iota , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  rcd_0nu=simulate_vivax_ode(parameters=parameters_rcd_0nu , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  rcd_0eta=simulate_vivax_ode(parameters=parameters_rcd_0eta , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  vivax_rcd=simulate_vivax_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)

  expect_equal(vivax_0[1465,], rcd_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax_0[1465,], rcd_0nu[1465,], tolerance = 1e-09, label = "nu=0")
  expect_equal(vivax_0[1465,], rcd_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_gt(vivax_0$I[1465], vivax_rcd$I[1465], label = "rcd")



  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.0155531,"delta"=0.01,
                      "alpha"=0.2, "beta"=0.7, "rho"=1,"omega"=0.9,
                      "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1,
                      "h"=0, "hr"=0,"hl"=0,  "hh"=0, "hhl"=0,
                      "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters2_rcd_0iota=parameters2_rcd; parameters2_rcd_0iota$iota=0
  parameters2_rcd_0nu=parameters2_rcd; parameters2_rcd_0nu$nu=0
  parameters2_rcd_0eta=parameters2_rcd; parameters2_rcd_0eta$eta=0

  vivax2_0=simulate_vivax_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_cm_importation_vc, maxtime=1465, year=FALSE)
  rcd2_0iota=simulate_vivax_ode(parameters=parameters2_rcd_0iota , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  rcd2_0nu=simulate_vivax_ode(parameters=parameters2_rcd_0nu , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  rcd2_0eta=simulate_vivax_ode(parameters=parameters2_rcd_0eta , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  vivax2_rcd=simulate_vivax_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)

  expect_equal(vivax2_0[1465,], rcd2_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax2_0[1465,], rcd2_0nu[1465,], tolerance = 1e-09, label = "nu=0")
  expect_equal(vivax2_0[1465,], rcd2_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_gt(vivax2_0$I[1465], vivax2_rcd$I[1465], label = "rcd")

})

test_that("test compare RCD in delay model (non referral) with non RCD model", {

  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.0155531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=1,"omega"=1,
                      "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0,
                      "sigma"=1/15, "T0"=0, "Tl"=0,
                      "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters_rcd_0iota=parameters_rcd; parameters_rcd_0iota$iota=0
  parameters_rcd_0nu=parameters_rcd; parameters_rcd_0nu$nu=0
  parameters_rcd_0eta=parameters_rcd; parameters_rcd_0eta$eta=0
  parameters_rcd_0sigma=parameters_rcd; parameters_rcd_0sigma$sigma=0

  vivax_0=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay, maxtime=1465, year=FALSE)
  rcd_0iota=simulate_vivax_delay_ode(parameters=parameters_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  rcd_0nu=simulate_vivax_delay_ode(parameters=parameters_rcd_0nu , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  rcd_0eta=simulate_vivax_delay_ode(parameters=parameters_rcd_0eta , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  rcd_0sigma=simulate_vivax_delay_ode(parameters=parameters_rcd_0sigma , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  vivax_rcd=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)

  expect_equal(vivax_0[1465,], rcd_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax_0[1465,], rcd_0nu[1465,], tolerance = 1e-09, label = "nu=0")
  expect_equal(vivax_0[1465,], rcd_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_equal(vivax_rcd, rcd_0sigma, tolerance = 1e-09, label = "sigma=0 and normal rcd")
  expect_gt(vivax_0$I[1465], vivax_rcd$I[1465], label = "rcd")

  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.0155531,"delta"=0.01,
                       "alpha"=0.2, "beta"=0.7, "rho"=1,"omega"=0.9,
                       "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hr"=0,"hl"=0,"hh"=0, "hhl"=0,
                       "sigma"=1/15, "T0"=0, "Tl"=0,
                       "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters2_rcd_nodelay=parameters2_rcd
  parameters2_rcd_nodelay$Il=parameters2_rcd$Ul
  parameters2_rcd_nodelay$I0=parameters2_rcd$U0
  parameters2_rcd_nodelay$Ul=NULL
  parameters2_rcd_nodelay$U0=NULL


  parameters2_rcd_0iota=parameters2_rcd; parameters2_rcd_0iota$iota=0
  parameters2_rcd_0nu=parameters2_rcd; parameters2_rcd_0nu$nu=0
  parameters2_rcd_0eta=parameters2_rcd; parameters2_rcd_0eta$eta=0
  parameters2_rcd_infsigma=parameters2_rcd; parameters2_rcd_infsigma$sigma=1/0.00001

  vivax2_0=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay, maxtime=1465, year=FALSE)
  rcd2_0iota=simulate_vivax_delay_ode(parameters=parameters2_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  rcd2_0nu=simulate_vivax_delay_ode(parameters=parameters2_rcd_0nu , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  rcd2_0eta=simulate_vivax_delay_ode(parameters=parameters2_rcd_0eta , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  vivax2_rcd=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  vivax2_rcd_infsigma=simulate_vivax_delay_ode(parameters=parameters2_rcd_infsigma , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  vivax2_rcd_nodel=simulate_vivax_ode(parameters=parameters2_rcd_nodelay , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)

  expect_equal(vivax2_0[1465,], rcd2_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax2_0[1465,], rcd2_0nu[1465,], tolerance = 1e-09, label = "nu=0")
  expect_equal(vivax2_0[1465,], rcd2_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_gt(vivax2_0$I[1465], vivax2_rcd$I[1465], label = "rcd")
  expect_equal(vivax2_rcd_infsigma[1465,]$I, vivax2_rcd_nodel[1465,]$I, tolerance = 1e-07, label = "sigma=inf")

})


test_that("test compare RCD in delay model (referral) with non RCD model", {

  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.0155531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=1,"omega"=1,
                      "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hr"=0,"hl"=0,"hh"=0, "hhl"=0,
                      "sigma"=1/15, "T0"=0, "Tl"=0,
                      "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters_rcd_0iota=parameters_rcd; parameters_rcd_0iota$iota=0
  parameters_rcd_0nu=parameters_rcd; parameters_rcd_0nu$nu=0
  parameters_rcd_0eta=parameters_rcd; parameters_rcd_0eta$eta=0
  parameters_rcd_0sigma=parameters_rcd; parameters_rcd_0sigma$sigma=0

  vivax_0=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay, maxtime=1465, year=FALSE)
  rcd_0iota=simulate_vivax_delay_ode(parameters=parameters_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  rcd_0nu=simulate_vivax_delay_ode(parameters=parameters_rcd_0nu , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  rcd_0eta=simulate_vivax_delay_ode(parameters=parameters_rcd_0eta , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  rcd_0sigma=simulate_vivax_delay_ode(parameters=parameters_rcd_0sigma , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  vivax_rcd=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)

  expect_equal(vivax_0[1465,], rcd_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax_0[1465,], rcd_0nu[1465,], tolerance = 1e-09, label = "nu=0")
  expect_equal(vivax_0[1465,], rcd_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_equal(vivax_0[1465,c("Sl", "S0", "I")], rcd_0sigma[1465,c("Sl", "S0", "I")], tolerance = 1e-09, label = "sigma=0")
  expect_gt(vivax_0$I[1465], vivax_rcd$I[1465], label = "rcd")
  expect_lt(vivax_rcd$I[1465], rcd_0sigma$I[1465],  label = "delay increases eq. prevalence")

  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.0155531,"delta"=0.01,
                       "alpha"=0.2, "beta"=0.7, "rho"=1,"omega"=0.9,
                       "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hr"=0,"hl"=0,"hh"=0, "hhl"=0,
                       "sigma"=1/15, "T0"=0, "Tl"=0,
                       "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters2_rcd_nodelay=parameters2_rcd
  parameters2_rcd_nodelay$Il=parameters2_rcd$Ul
  parameters2_rcd_nodelay$I0=parameters2_rcd$U0
  parameters2_rcd_nodelay$Ul=NULL
  parameters2_rcd_nodelay$U0=NULL

  parameters2_rcd_0iota=parameters2_rcd; parameters2_rcd_0iota$iota=0
  parameters2_rcd_0nu=parameters2_rcd; parameters2_rcd_0nu$nu=0
  parameters2_rcd_0eta=parameters2_rcd; parameters2_rcd_0eta$eta=0
  parameters2_rcd_0sigma=parameters2_rcd; parameters2_rcd_0sigma$sigma=0
  parameters2_rcd_infsigma=parameters2_rcd; parameters2_rcd_infsigma$sigma=1/0.00001

  vivax2_0=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay, maxtime=1465, year=FALSE)
  rcd2_0iota=simulate_vivax_delay_ode(parameters=parameters2_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  rcd2_0nu=simulate_vivax_delay_ode(parameters=parameters2_rcd_0nu , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  rcd2_0eta=simulate_vivax_delay_ode(parameters=parameters2_rcd_0eta , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  rcd2_0sigma=simulate_vivax_delay_ode(parameters=parameters2_rcd_0sigma , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  vivax2_rcd=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  vivax2_rcd_infsigma=simulate_vivax_delay_ode(parameters=parameters2_rcd_infsigma , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  vivax2_rcd_nodel=simulate_vivax_ode(parameters=parameters2_rcd_nodelay , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)

  expect_equal(vivax2_0[1465,], rcd2_0iota[1465,], tolerance = 1e-09, label = "iota=0")
  expect_equal(vivax2_0[1465,], rcd2_0nu[1465,], tolerance = 1e-09, label = "nu=0")
  expect_equal(vivax2_0[1465,], rcd2_0eta[1465,], tolerance = 1e-09, label = "eta=0")
  expect_gt(vivax2_0$I[1465], vivax2_rcd$I[1465], label = "rcd")
  expect_lt(vivax2_rcd$I[1465], rcd2_0sigma$I[1465],  label = "delay increases eq. prevalence")
  expect_equal(vivax2_rcd_infsigma[1465,]$I, vivax2_rcd_nodel[1465,]$I, tolerance = 1e-07,  label = "sigma=inf")

})

test_that("test compare the 3 RCD models", {


  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.0155531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=1,"omega"=1,
                      "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0,
                      "sigma"=1/15, "T0"=0, "Tl"=0,
                      "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters_rcd_nodelay=parameters_rcd
  parameters_rcd_nodelay$Il=parameters_rcd$Ul
  parameters_rcd_nodelay$I0=parameters_rcd$U0
  parameters_rcd_nodelay$Ul=NULL
  parameters_rcd_nodelay$U0=NULL

  parameters_rcd_infsigma=parameters_rcd; parameters_rcd_infsigma$sigma=1/0.0001

  vivax_rcd=simulate_vivax_ode(parameters=parameters_rcd_nodelay , ODEmodel=ode_vivax_rcd, maxtime=1465, year=FALSE)
  vivax_delay_non_ref=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=1465, year=FALSE)
  vivax_delay_ref=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)
  vivax_delay_ref_siginf=simulate_vivax_delay_ode(parameters=parameters_rcd_infsigma , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=1465, year=FALSE)

  expect_equal(vivax_delay_non_ref$I[1465], vivax_rcd$I[1465], tolerance = 1e-09, label = "rcd non ref = rcd no delay")
  expect_lt(vivax_delay_non_ref$I[1465], vivax_delay_ref$I[1465], label = "rcd non ref is more effective")
  expect_equal(vivax_delay_non_ref$I[1465], vivax_delay_ref_siginf$I[1465], tolerance = 1e-06, label = "rcd non ref approx rcd ref if sigma =inf")
  expect_equal(vivax_rcd$I[1465], vivax_delay_ref_siginf$I[1465], tolerance = 1e-06, label = "rcd non ref approx no delay if sigma =inf")


  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.0155531,"delta"=0,
                       "alpha"=0.2, "beta"=0.7, "rho"=1,"omega"=0.9,
                       "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0,
                       "sigma"=1/15, "T0"=0, "Tl"=0,
                       "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters2_rcd_nodelay=parameters2_rcd
  parameters2_rcd_nodelay$Il=parameters2_rcd$Ul
  parameters2_rcd_nodelay$I0=parameters2_rcd$U0
  parameters2_rcd_nodelay$Ul=NULL
  parameters2_rcd_nodelay$U0=NULL

  parameters2_rcd_infsigma=parameters2_rcd; parameters2_rcd_infsigma$sigma=1/0.00001
  parameters2_rcd_0sigma=parameters2_rcd; parameters2_rcd_0sigma$sigma=0
  parameters2_rcd_0sigma_iota0=parameters2_rcd; parameters2_rcd_0sigma_iota0$sigma=parameters2_rcd$tau*parameters2_rcd$iota*parameters2_rcd$nu*parameters2_rcd$eta
  parameters2_rcd_0sigma_iota0$alpha=1-parameters2_rcd$alpha


  vivax2_rcd=simulate_vivax_ode(parameters=parameters2_rcd_nodelay , ODEmodel=ode_vivax_rcd, maxtime=5000, year=FALSE)
  vivax2_delay_non_ref=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  vivax2_delay_ref=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)
  vivax2_delay_non_ref_sig0=simulate_vivax_delay_ode(parameters=parameters2_rcd_0sigma , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  vivax2_delay_non_ref_siginf=simulate_vivax_delay_ode(parameters=parameters2_rcd_infsigma , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  vivax2_delay_ref_siginf=simulate_vivax_delay_ode(parameters=parameters2_rcd_infsigma , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)
  vivax2_delay_ref_iota0=simulate_vivax_delay_ode(parameters=parameters2_rcd_0sigma_iota0 , ODEmodel=ode_vivax_delay, maxtime=5000, year=FALSE)
  names(vivax2_delay_ref_iota0)=c("time", "Tl",   "T0" ,  "Sl","S0"  , "Ul",   "U0" ,  "h" ,   "hl"  , "hh",   "hhl",  "I"  ,  "p") # swap IL and I0 with Tl and T0, and rename to U

  expect_gt(vivax2_delay_non_ref$I[5000], vivax2_rcd$I[5000], label = "no delay is more effective")
  expect_lt(vivax2_delay_non_ref$I[5000], vivax2_delay_ref$I[5000], label = "rcd non ref is more effective")
  expect_equal(vivax2_delay_non_ref_siginf$I[5000], vivax2_delay_ref_siginf$I[5000], tolerance = 1e-06, label = "rcd non ref approx rcd ref if sigma =inf")
  expect_equal(vivax2_rcd$I[5000], vivax2_delay_ref_siginf$I[5000], tolerance = 1e-06, label = "rcd ref approx no delay if sigma =inf")
  expect_equal(vivax2_delay_ref_iota0$I[5000], vivax2_delay_non_ref_sig0$I[5000], tolerance = 1e-06, label = "if sigma=0, equiv to non rcd with sigma =iota*nu*tau*eta and alpha replaced by 1-alpha")
  expect_equal(vivax2_delay_ref_iota0[5000,names(vivax2_delay_ref_iota0)[-c(8,9)]], vivax2_delay_non_ref_sig0[5000,names(vivax2_delay_ref_iota0)[-c(8,9)]], tolerance = 1e-06,
               label = "if sigma=0, equiv to non rcd with sigma =iota*nu*tau*eta and alpha replaced by 1-alpha. h and hl are not tested because Ul/U0 and Tl/T0 have been swapped")



  parameters3_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.0155531,"delta"=0.01,
                       "alpha"=0.2, "beta"=0, "rho"=1,"omega"=0.9,
                       "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.05, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0,
                       "sigma"=1/15, "T0"=0, "Tl"=0.05,
                       "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters3_rcd_rep=parameters3_rcd
  parameters3_rcd_rep$r=parameters3_rcd$r+parameters3_rcd$iota*parameters3_rcd$nu*parameters3_rcd$tau*parameters3_rcd$eta
  parameters3_rcd_rep$sigma=parameters3_rcd$sigma-parameters3_rcd$iota*parameters3_rcd$nu*parameters3_rcd$tau*parameters3_rcd$eta


  parameters4_rcd=parameters3_rcd
  parameters4_rcd$lambda=0.0955531
  parameters4_rcd$iota=parameters4_rcd$iota*50
  parameters4_rcd_rep=parameters4_rcd
  parameters4_rcd_rep$r=parameters4_rcd$r+parameters4_rcd$sigma
  parameters4_rcd_rep$sigma=-parameters4_rcd$sigma+parameters4_rcd$iota*parameters4_rcd$nu*parameters4_rcd$tau*parameters4_rcd$eta
  parameters4_rcd_rep$alpha=1-parameters4_rcd$alpha

  parameters5_rcd=parameters3_rcd
  parameters5_rcd$sigma=parameters5_rcd$iota*parameters5_rcd$nu*parameters5_rcd$tau*parameters5_rcd$eta
  parameters5_rcd_rep=parameters5_rcd
  parameters5_rcd_rep$r=parameters5_rcd$r+parameters5_rcd$sigma
  parameters5_rcd_rep$alpha=0
  parameters5_rcd_rep$Il=0.1
  parameters5_rcd_rep$I0=parameters5_rcd$U0

  vivax3_rcd=simulate_vivax_delay_ode(parameters=parameters3_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  vivax3_rcd_rep=simulate_vivax_delay_ode(parameters=parameters3_rcd_rep , ODEmodel=ode_vivax_delay, maxtime=5000, year=FALSE)

  vivax4_rcd=simulate_vivax_delay_ode(parameters=parameters4_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  vivax4_rcd_rep=simulate_vivax_delay_ode(parameters=parameters4_rcd_rep , ODEmodel=ode_vivax_delay, maxtime=5000, year=FALSE)
  names(vivax4_rcd_rep)=c("time", "Tl",   "T0" ,  "Sl","S0"  , "Ul",   "U0" ,  "h" ,   "hl"  , "hh",   "hhl",  "I"  ,  "p") # swap IL and I0 with Tl and T0

  vivax5_rcd=simulate_vivax_delay_ode(parameters=parameters5_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE) %>%
    dplyr::mutate(Il=Ul+Tl, I0=U0+T0) %>% dplyr::select(-T0,-Tl, -U0, -Ul)
  vivax5_rcd_rep=simulate_vivax_ode(parameters=parameters5_rcd_rep , ODEmodel=ode_vivax_cm_importation_vc, maxtime=5000, year=FALSE)%>% dplyr::select(-hr)

  expect_equal(vivax3_rcd_rep%>% dplyr::select(-h,-hl, -p), vivax3_rcd%>% dplyr::select(-h,-hl, -p), tolerance = 1e-08, label = "rcd non ref with beta=0 and sigma>rcd_term, fixed iota")
  expect_equal(vivax4_rcd_rep[names(vivax4_rcd)]%>% dplyr::select(-h,-hl, -p), vivax4_rcd%>% dplyr::select( -h,-hl, -p), tolerance = 1e-08, label = "rcd non ref with beta=0 and sigma<rcd_term, fixed iota")
  expect_equal(vivax5_rcd_rep%>% dplyr::select(-h,-hl, -p), vivax5_rcd[names(vivax5_rcd_rep)]%>% dplyr::select(-h,-hl, -p), tolerance = 1e-09, label = "rcd non ref with beta=0 and sigma=rcd_term, fixed iota")

})

test_that("test solve with RCD model, without any delays", {

  # first parameter set (no importation, no CM)
  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.0155531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=0.5,"omega"=1,
                      "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1,
                      "h"=0, "hr"=0,"hl"=0, "hh"=0, "hhl"=0,
                      "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters_rcd_0iota=parameters_rcd; parameters_rcd_0iota$iota=0

  rcd_0iota=simulate_vivax_ode(parameters=parameters_rcd_0iota , ODEmodel=ode_vivax_rcd, maxtime=5000, year=FALSE)
  vivax_rcd=simulate_vivax_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_rcd, maxtime=5000, year=FALSE)

  h1=calculate_h1_rcd(h=vivax_rcd$h[5000], alpha=parameters_rcd$alpha,
                      rho=parameters_rcd$rho,
                      iota=parameters_rcd$iota, nu=parameters_rcd$nu,
                      eta=parameters_rcd$eta, tau=parameters_rcd$tau,
                       r=parameters_rcd$r, rho2=parameters_rcd$rho2)
  expect_equal(h1, vivax_rcd$hh[5000]*parameters_rcd$rho)

  mylambda_calc=solve_lambda_vivax_rcd(h=vivax_rcd$h[5000],h1=h1,
                                       r=parameters_rcd$r,f=parameters_rcd$f,
                                       gamma=parameters_rcd$gamma,
                                       alpha=parameters_rcd$alpha, beta=parameters_rcd$beta,
                                       rho=parameters_rcd$rho, p=vivax_rcd$p[5000],
                                       omega=parameters_rcd$omega,
                                       iota=parameters_rcd$iota, nu=parameters_rcd$nu,
                                       eta=parameters_rcd$eta, tau=parameters_rcd$tau,
                                       return.all = T)

  h1_rcd_0iota=calculate_h1_rcd(h=rcd_0iota$h[5000], alpha=parameters_rcd_0iota$alpha,
                      rho=parameters_rcd_0iota$rho,
                      iota=parameters_rcd_0iota$iota, nu=parameters_rcd_0iota$nu,
                      eta=parameters_rcd_0iota$eta, tau=parameters_rcd_0iota$tau,
                      r=parameters_rcd_0iota$r,rho2=parameters_rcd_0iota$rho2 )
  mylambda_calc0=solve_lambda_vivax_rcd(h=rcd_0iota$h[5000],h1=h1_rcd_0iota,
                                        r=parameters_rcd_0iota$r,f=parameters_rcd_0iota$f,
                                        gamma=parameters_rcd_0iota$gamma,
                                        alpha=parameters_rcd_0iota$alpha, beta=parameters_rcd_0iota$beta,
                                        rho=parameters_rcd_0iota$rho, p=rcd_0iota$p[5000],
                                        omega=parameters_rcd_0iota$omega,
                                        iota=parameters_rcd_0iota$iota, nu=parameters_rcd_0iota$nu,
                                        eta=parameters_rcd_0iota$eta, tau=parameters_rcd_0iota$tau,
                                        return.all = T)


  expect_equal(mylambda_calc$lambda_real, parameters_rcd$lambda, tolerance = 1e-09, label = "lambda, rcd no delay no CM")
  expect_equal(mylambda_calc0$lambda_real, parameters_rcd_0iota$lambda, tolerance = 1e-09, label = "lambda, rcd=0 no delay no CM")


  eq_vivax_rcd=get_equilibrium_states_vivax_rcd(I=vivax_rcd$I[5000],
                                            lambda=parameters_rcd$lambda,
                                            r=parameters_rcd$r, gamma=parameters_rcd$gamma, f=parameters_rcd$f,
                                            alpha=parameters_rcd$alpha,  beta=parameters_rcd$beta, rho=parameters_rcd$rho,
                                            delta=parameters_rcd$delta, omega=parameters_rcd$omega,
                                            iota_star=mylambda_calc$iota_star, nu=parameters_rcd$nu, tau=parameters_rcd$tau,
                                            eta=parameters_rcd$eta, rho2=parameters_rcd$rho2)
  eq_vivax_rcd_0iota=get_equilibrium_states_vivax_rcd(I=rcd_0iota$I[5000],
                                                      lambda=parameters_rcd_0iota$lambda,
                                                      r=parameters_rcd_0iota$r, gamma=parameters_rcd_0iota$gamma, f=parameters_rcd_0iota$f,
                                                      alpha=parameters_rcd_0iota$alpha,  beta=parameters_rcd_0iota$beta, rho=parameters_rcd_0iota$rho,
                                                      delta=parameters_rcd_0iota$delta, omega=parameters_rcd_0iota$omega,
                                                      iota_star=mylambda_calc0$iota_star, nu=parameters_rcd_0iota$nu, tau=parameters_rcd_0iota$tau, eta=parameters_rcd_0iota$eta, rho2=parameters_rcd_0iota$rho2)
  expect_equal(eq_vivax_rcd, as.list(vivax_rcd[5000,names(eq_vivax_rcd)]), tolerance = 1e-09, label = "equilibrium values, rcd no delay")
  expect_equal(eq_vivax_rcd_0iota, as.list(rcd_0iota[5000,names(eq_vivax_rcd_0iota)]), tolerance = 1e-09, label = "equilibrium values, rcd=0, no delay")


  my_tau=calculate_tau_rcd(h1=vivax_rcd$hh[5000]*parameters_rcd$rho, h2=vivax_rcd$h[5000]-vivax_rcd$hh[5000]*parameters_rcd$rho,
                    alpha=parameters_rcd$alpha,
                    rho=parameters_rcd$rho,
                    iota=parameters_rcd$iota, nu=parameters_rcd$nu,
                    eta=parameters_rcd$eta,
                    r=parameters_rcd$r, rho2 = parameters_rcd$rho2)
  expect_equal(my_tau,parameters_rcd$tau, tolerance = 1e-09, label = "calculating tau with h2 and h1")

  # second parameter set (with CM and importation)
  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.01631,"delta"=0.01,
                       "alpha"=0.2, "beta"=0.7, "rho"=0.6,"omega"=0.9,
                       "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1,
                       "h"=0, "hr"=0,"hl"=0, "hh"=0, "hhl"=0,
                       "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters2_rcd_0iota=parameters2_rcd; parameters2_rcd_0iota$iota=0

  vivax2_rcd=simulate_vivax_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_rcd, maxtime=5000, year=FALSE)
  rcd2_0iota=simulate_vivax_ode(parameters=parameters2_rcd_0iota , ODEmodel=ode_vivax_rcd, maxtime=5000, year=FALSE)

  h1_2=calculate_h1_rcd(h=vivax2_rcd$h[5000], alpha=parameters2_rcd$alpha,
                      rho=parameters2_rcd$rho,
                      iota=parameters2_rcd$iota, nu=parameters2_rcd$nu,
                      eta=parameters2_rcd$eta, tau=parameters2_rcd$tau,
                      r=parameters2_rcd$r,rho2=parameters2_rcd$rho2)

  expect_equal(h1_2, vivax2_rcd$hh[5000]*parameters2_rcd$rho)
  mylambda_calc2=solve_lambda_vivax_rcd(h=vivax2_rcd$h[5000],h1=h1_2,
                                        r=parameters2_rcd$r,f=parameters2_rcd$f,
                                        gamma=parameters2_rcd$gamma,
                                        alpha=parameters2_rcd$alpha, beta=parameters2_rcd$beta,
                                        rho=parameters2_rcd$rho, p=vivax2_rcd$p[5000],
                                        omega=parameters2_rcd$omega,
                                        iota=parameters2_rcd$iota, nu=parameters2_rcd$nu,
                                        eta=parameters2_rcd$eta, tau=parameters2_rcd$tau,
                                        return.all = T)

  h1_2=calculate_h1_rcd(h=rcd2_0iota$h[5000], alpha=parameters2_rcd_0iota$alpha,
                        rho=parameters2_rcd_0iota$rho,
                        iota=parameters2_rcd_0iota$iota, nu=parameters2_rcd_0iota$nu,
                        eta=parameters2_rcd_0iota$eta, tau=parameters2_rcd_0iota$tau,
                        r=parameters2_rcd_0iota$r,rho2=parameters2_rcd_0iota$rho2 )

  mylambda_calc0_2=solve_lambda_vivax_rcd(h=rcd2_0iota$h[5000],h1=h1_2,
                                          r=parameters2_rcd_0iota$r,f=parameters2_rcd_0iota$f,
                                          gamma=parameters2_rcd_0iota$gamma,
                                          alpha=parameters2_rcd_0iota$alpha, beta=parameters2_rcd_0iota$beta,
                                          rho=parameters2_rcd_0iota$rho, p=rcd2_0iota$p[5000],
                                          omega=parameters2_rcd_0iota$omega,
                                          iota=parameters2_rcd_0iota$iota, nu=parameters2_rcd_0iota$nu,
                                          eta=parameters2_rcd_0iota$eta, tau=parameters2_rcd_0iota$tau,
                                          return.all = T)


  expect_equal(mylambda_calc2$lambda_real, parameters2_rcd$lambda, tolerance = 1e-09, label = "lambda, rcd no delay with CM")
  expect_equal(mylambda_calc0_2$lambda_real, parameters2_rcd_0iota$lambda, tolerance = 1e-09, label = "lambda, rcd=0 no delay with CM")

  eq_vivax_rcd2=get_equilibrium_states_vivax_rcd(I=vivax2_rcd$I[5000],
                                                lambda=parameters2_rcd$lambda,
                                                r=parameters2_rcd$r, gamma=parameters2_rcd$gamma, f=parameters2_rcd$f,
                                                alpha=parameters2_rcd$alpha,  beta=parameters2_rcd$beta, rho=parameters2_rcd$rho,
                                                delta=parameters2_rcd$delta, omega=parameters2_rcd$omega,
                                                iota_star=mylambda_calc2$iota_star, nu=parameters2_rcd$nu,
                                                tau=parameters2_rcd$tau, eta=parameters2_rcd$eta,rho2=parameters2_rcd$rho2)
  eq_vivax_rcd2_0iota=get_equilibrium_states_vivax_rcd(I=rcd2_0iota$I[5000],
                                                      lambda=parameters2_rcd_0iota$lambda,
                                                      r=parameters2_rcd_0iota$r, gamma=parameters2_rcd_0iota$gamma, f=parameters2_rcd_0iota$f,
                                                      alpha=parameters2_rcd_0iota$alpha,  beta=parameters2_rcd_0iota$beta, rho=parameters2_rcd_0iota$rho,
                                                      delta=parameters2_rcd_0iota$delta, omega=parameters2_rcd_0iota$omega,
                                                      iota_star=mylambda_calc0_2$iota_star, nu=parameters2_rcd_0iota$nu,
                                                      tau=parameters2_rcd_0iota$tau, eta=parameters2_rcd_0iota$eta, rho2=parameters2_rcd_0iota$rho2)
  expect_equal(eq_vivax_rcd2, as.list(vivax2_rcd[5000,names(eq_vivax_rcd2)]), tolerance = 1e-09, label = "equilibrium values, rcd no delay")
  expect_equal(eq_vivax_rcd2_0iota, as.list(rcd2_0iota[5000,names(eq_vivax_rcd2_0iota)]), tolerance = 1e-09, label = "equilibrium values, rcd=0, no delay")

  my_tau2=calculate_tau_rcd(h1=vivax2_rcd$hh[5000]*parameters2_rcd$rho, h2=vivax2_rcd$h[5000]-vivax2_rcd$hh[5000]*parameters2_rcd$rho,
                           alpha=parameters2_rcd$alpha,
                           rho=parameters2_rcd$rho,
                           iota=parameters2_rcd$iota, nu=parameters2_rcd$nu,
                           eta=parameters2_rcd$eta,
                           r=parameters2_rcd$r, rho2=parameters2_rcd$rho2)
  expect_equal(my_tau2,parameters2_rcd$tau, tolerance = 1e-09, label = "calculating tau with h2 and h1")
})

test_that("test solve RCD in delay model (non referral)", {

  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.0155531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=0.7,"omega"=1,
                      "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0, "hh"=0, "hhl"=0,
                      "sigma"=1/15, "T0"=0, "Tl"=0,
                      "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters_rcd_0iota=parameters_rcd; parameters_rcd_0iota$iota=0
  parameters_rcd_0nu=parameters_rcd; parameters_rcd_0nu$nu=0
  parameters_rcd_0sigma=parameters_rcd; parameters_rcd_0sigma$sigma=0

  rcd_0iota=simulate_vivax_delay_ode(parameters=parameters_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  rcd_0nu=simulate_vivax_delay_ode(parameters=parameters_rcd_0nu , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  rcd_0sigma=simulate_vivax_delay_ode(parameters=parameters_rcd_0sigma , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  vivax_rcd=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)

  h1=calculate_h1_rcd(h=vivax_rcd$h[5000], alpha=parameters_rcd$alpha,
                      rho=parameters_rcd$rho,
                      iota=parameters_rcd$iota, nu=parameters_rcd$nu,
                      eta=parameters_rcd$eta, tau=parameters_rcd$tau,
                      r=parameters_rcd$r, rho2=parameters_rcd$rho2 )
  expect_equal(h1, vivax_rcd$hh[5000]*parameters_rcd$rho)


  mylambda_calc=solve_lambda_vivax_rcd_no_referral(h=vivax_rcd$h[5000],h1=h1,
                                                   r=parameters_rcd$r,f=parameters_rcd$f,
                                                   gamma=parameters_rcd$gamma,
                                                   alpha=parameters_rcd$alpha, beta=parameters_rcd$beta,
                                                   rho=parameters_rcd$rho, sigma=parameters_rcd$sigma,
                                                   p=vivax_rcd$p[5000],
                                                   omega=parameters_rcd$omega,
                                                   iota=parameters_rcd$iota, nu=parameters_rcd$nu,
                                                   eta=parameters_rcd$eta, tau=parameters_rcd$tau, return.all = T)

  h1_rcd_0iota=calculate_h1_rcd(h=rcd_0iota$h[5000], alpha=parameters_rcd_0iota$alpha,
                                rho=parameters_rcd_0iota$rho,
                                iota=parameters_rcd_0iota$iota, nu=parameters_rcd_0iota$nu,
                                eta=parameters_rcd_0iota$eta, tau=parameters_rcd_0iota$tau,
                                r=parameters_rcd_0iota$r, rho2=parameters_rcd_0iota$rho2 )
  mylambda_calc0=solve_lambda_vivax_rcd_no_referral(h=rcd_0iota$h[5000], h1=h1_rcd_0iota,
                                                    r=parameters_rcd_0iota$r,f=parameters_rcd_0iota$f,
                                                    gamma=parameters_rcd_0iota$gamma,
                                                    alpha=parameters_rcd_0iota$alpha, beta=parameters_rcd_0iota$beta,
                                                    rho=parameters_rcd_0iota$rho, sigma=parameters_rcd$sigma,
                                                    p=rcd_0iota$p[5000],
                                                    omega=parameters_rcd_0iota$omega,
                                                    iota=parameters_rcd_0iota$iota, nu=parameters_rcd_0iota$nu,
                                                    eta=parameters_rcd_0iota$eta, tau=parameters_rcd_0iota$tau, return.all = T)

  h1_rcd_0nu=calculate_h1_rcd(h=rcd_0nu$h[5000], alpha=parameters_rcd_0nu$alpha,
                                rho=parameters_rcd_0nu$rho,
                                iota=parameters_rcd_0nu$iota, nu=parameters_rcd_0nu$nu,
                                eta=parameters_rcd_0nu$eta, tau=parameters_rcd_0nu$tau,
                                r=parameters_rcd_0nu$r , rho2=parameters_rcd_0nu$rho2 )
  mylambda_calc0nu=solve_lambda_vivax_rcd_no_referral(h=rcd_0nu$h[5000],h1=h1_rcd_0nu,
                                                      r=parameters_rcd_0nu$r,f=parameters_rcd_0nu$f,
                                                      gamma=parameters_rcd_0nu$gamma,
                                                      alpha=parameters_rcd_0nu$alpha, beta=parameters_rcd_0nu$beta,
                                                      rho=parameters_rcd_0nu$rho, sigma=parameters_rcd_0nu$sigma,
                                                      p=rcd_0nu$p[5000],
                                                      omega=parameters_rcd_0nu$omega,
                                                      iota=parameters_rcd_0nu$iota, nu=parameters_rcd_0nu$nu,
                                                      eta=parameters_rcd_0nu$eta, tau=parameters_rcd_0nu$tau, return.all = T)


  expect_equal(mylambda_calc$lambda_real, parameters_rcd$lambda, tolerance = 1e-09, label = "lambda, rcd no delay no CM")
  expect_equal(mylambda_calc0$lambda_real, parameters_rcd_0iota$lambda, tolerance = 1e-09, label = "lambda, rcd=0 no delay no CM")
  expect_equal(mylambda_calc0nu$lambda_real, parameters_rcd_0nu$lambda, tolerance = 1e-09, label = "lambda, rcd=0 no delay no CM")

  eq_vivax_rcd=get_equilibrium_states_vivax_rcd_no_referral(I=vivax_rcd$I[5000],
                                                            lambda=parameters_rcd$lambda, sigma=parameters_rcd$sigma,
                                                            r=parameters_rcd$r, gamma=parameters_rcd$gamma, f=parameters_rcd$f,
                                                            alpha=parameters_rcd$alpha,  beta=parameters_rcd$beta, rho=parameters_rcd$rho,
                                                            delta=parameters_rcd$delta, omega=parameters_rcd$omega,
                                                            iota_star=mylambda_calc$iota_star, nu=parameters_rcd$nu,
                                                            tau=parameters_rcd$tau, eta=parameters_rcd$eta, rho2=parameters_rcd$rho2)
  eq_vivax_rcd_0iota=get_equilibrium_states_vivax_rcd_no_referral(I=rcd_0iota$I[5000],
                                                                  lambda=parameters_rcd_0iota$lambda, sigma=parameters_rcd_0iota$sigma,
                                                                  r=parameters_rcd_0iota$r, gamma=parameters_rcd_0iota$gamma, f=parameters_rcd_0iota$f,
                                                                  alpha=parameters_rcd_0iota$alpha,  beta=parameters_rcd_0iota$beta, rho=parameters_rcd_0iota$rho,
                                                                  delta=parameters_rcd_0iota$delta, omega=parameters_rcd_0iota$omega,
                                                                  iota_star=mylambda_calc0$iota_star, nu=parameters_rcd_0iota$nu,
                                                                  tau=parameters_rcd_0iota$tau, eta=parameters_rcd_0iota$eta, rho2=parameters_rcd_0iota$rho2)
  eq_vivax_rcd_0nu=get_equilibrium_states_vivax_rcd_no_referral(I=rcd_0nu$I[5000],
                                                                  lambda=parameters_rcd_0nu$lambda, sigma=parameters_rcd_0nu$sigma,
                                                                  r=parameters_rcd_0nu$r, gamma=parameters_rcd_0nu$gamma, f=parameters_rcd_0nu$f,
                                                                  alpha=parameters_rcd_0nu$alpha,  beta=parameters_rcd_0nu$beta, rho=parameters_rcd_0nu$rho,
                                                                  delta=parameters_rcd_0nu$delta, omega=parameters_rcd_0nu$omega,
                                                                  iota_star=mylambda_calc0nu$iota_star, nu=parameters_rcd_0nu$nu,
                                                                  tau=parameters_rcd_0nu$tau, eta=parameters_rcd_0nu$eta, rho2=parameters_rcd_0nu$rho2)
  expect_equal(eq_vivax_rcd, as.list(vivax_rcd[5000,names(eq_vivax_rcd)]), tolerance = 1e-09, label = "equilibrium values, rcd no delay")
  expect_equal(eq_vivax_rcd_0iota, as.list(rcd_0iota[5000,names(eq_vivax_rcd_0iota)]), tolerance = 1e-09, label = "equilibrium values, rcd=0, no delay")
  expect_equal(eq_vivax_rcd_0nu, as.list(rcd_0nu[5000,names(eq_vivax_rcd_0nu)]), tolerance = 1e-09, label = "equilibrium values, rcd=0, no delay")
  expect_equal(eq_vivax_rcd_0nu, eq_vivax_rcd_0iota, tolerance = 1e-09, label = "equilibrium values, rcd=0, no delay")

  my_tau=calculate_tau_rcd(h1=vivax_rcd$hh[5000]*parameters_rcd$rho, h2=vivax_rcd$h[5000]-vivax_rcd$hh[5000]*parameters_rcd$rho,
                           alpha=parameters_rcd$alpha,
                           rho=parameters_rcd$rho,
                           iota=parameters_rcd$iota, nu=parameters_rcd$nu,
                           eta=parameters_rcd$eta,
                           r=parameters_rcd$r, rho2=parameters_rcd$rho2)
  expect_equal(my_tau,parameters_rcd$tau, tolerance = 1e-09, label = "calculating tau with h2 and h1")

  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.0155531,"delta"=0.01,
                       "alpha"=0.2, "beta"=0.7, "rho"=0.5,"omega"=0.9,
                       "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0, "hh"=0, "hhl"=0,
                       "sigma"=1/15, "T0"=0, "Tl"=0,
                       "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters2_rcd_0iota=parameters2_rcd; parameters2_rcd_0iota$iota=0
  parameters2_rcd_0sigma=parameters2_rcd; parameters2_rcd_0sigma$sigma=0
  parameters2_rcd_infsigma=parameters2_rcd; parameters2_rcd_infsigma$sigma=1/0.4

  rcd2_0iota=simulate_vivax_delay_ode(parameters=parameters2_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  vivax2_rcd=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)

  h1_2=calculate_h1_rcd(h=vivax2_rcd$h[5000], alpha=parameters2_rcd$alpha,
                      rho=parameters2_rcd$rho,
                      iota=parameters2_rcd$iota, nu=parameters2_rcd$nu,
                      eta=parameters2_rcd$eta, tau=parameters2_rcd$tau,
                      r=parameters2_rcd$r , rho2=parameters2_rcd$rho2 )
  expect_equal(h1_2, vivax2_rcd$hh[5000]*parameters2_rcd$rho)

  mylambda_calc2=solve_lambda_vivax_rcd_no_referral(h=vivax2_rcd$h[5000],h1=h1_2,
                                                    r=parameters2_rcd$r,f=parameters2_rcd$f,
                                                    gamma=parameters2_rcd$gamma,
                                                    alpha=parameters2_rcd$alpha, beta=parameters2_rcd$beta,
                                                    rho=parameters2_rcd$rho, sigma=parameters2_rcd$sigma,
                                                    p=vivax2_rcd$p[5000],
                                                    omega=parameters2_rcd$omega,
                                                    iota=parameters2_rcd$iota, nu=parameters2_rcd$nu,
                                                    eta=parameters2_rcd$eta, tau=parameters2_rcd$tau,
                                                    return.all = T)
  h1_0iota2=calculate_h1_rcd(h=rcd2_0iota$h[5000], alpha=parameters2_rcd_0iota$alpha,
                        rho=parameters2_rcd_0iota$rho,
                        iota=parameters2_rcd_0iota$iota, nu=parameters2_rcd_0iota$nu,
                        eta=parameters2_rcd_0iota$eta, tau=parameters2_rcd_0iota$tau,
                        r=parameters2_rcd_0iota$r , rho2=parameters2_rcd_0iota$rho2 )

  mylambda_calc0_2=solve_lambda_vivax_rcd_no_referral(h=rcd2_0iota$h[5000],h1=h1_0iota2,
                                                      r=parameters2_rcd_0iota$r,f=parameters2_rcd_0iota$f,
                                                      gamma=parameters2_rcd_0iota$gamma,
                                                      alpha=parameters2_rcd_0iota$alpha, beta=parameters2_rcd_0iota$beta,
                                                      rho=parameters2_rcd_0iota$rho, sigma=parameters2_rcd_0iota$sigma,
                                                      p=rcd2_0iota$p[5000],
                                                      omega=parameters2_rcd_0iota$omega,
                                                      iota=parameters2_rcd_0iota$iota, nu=parameters2_rcd_0iota$nu,
                                                      eta=parameters2_rcd_0iota$eta, tau=parameters2_rcd_0iota$tau,
                                                      return.all = T)


  expect_equal(mylambda_calc2$lambda_real, parameters2_rcd$lambda, tolerance = 1e-09, label = "lambda, rcd no delay with CM")
  expect_equal(mylambda_calc0_2$lambda_real, parameters2_rcd_0iota$lambda, tolerance = 1e-09, label = "lambda, rcd=0 no delay with CM")

  eq_vivax_rcd2=get_equilibrium_states_vivax_rcd_no_referral(I=vivax2_rcd$I[5000],
                                                             lambda=parameters2_rcd$lambda, sigma=parameters2_rcd$sigma,
                                                             r=parameters2_rcd$r, gamma=parameters2_rcd$gamma, f=parameters2_rcd$f,
                                                             alpha=parameters2_rcd$alpha,  beta=parameters2_rcd$beta, rho=parameters2_rcd$rho,
                                                             delta=parameters2_rcd$delta, omega=parameters2_rcd$omega,
                                                             iota_star=mylambda_calc2$iota_star, nu=parameters2_rcd$nu,
                                                             tau=parameters2_rcd$tau, eta=parameters2_rcd$eta, rho2=parameters2_rcd$rho2)
  eq_vivax_rcd2_0iota=get_equilibrium_states_vivax_rcd_no_referral(I=rcd2_0iota$I[5000],
                                                                   lambda=parameters2_rcd_0iota$lambda, sigma=parameters2_rcd_0iota$sigma,
                                                                   r=parameters2_rcd_0iota$r, gamma=parameters2_rcd_0iota$gamma, f=parameters2_rcd_0iota$f,
                                                                   alpha=parameters2_rcd_0iota$alpha,  beta=parameters2_rcd_0iota$beta, rho=parameters2_rcd_0iota$rho,
                                                                   delta=parameters2_rcd_0iota$delta, omega=parameters2_rcd_0iota$omega,
                                                                   iota_star=mylambda_calc0_2$iota_star, nu=parameters2_rcd_0iota$nu,
                                                                   tau=parameters2_rcd_0iota$tau, eta=parameters2_rcd_0iota$eta, rho2=parameters2_rcd_0iota$rho2)
  expect_equal(eq_vivax_rcd2, as.list(vivax2_rcd[5000,names(eq_vivax_rcd2)]), tolerance = 1e-09, label = "equilibrium values, rcd no delay")
  expect_equal(eq_vivax_rcd2_0iota, as.list(rcd2_0iota[5000,names(eq_vivax_rcd2_0iota)]), tolerance = 1e-09, label = "equilibrium values, rcd=0, no delay")

  my_tau2=calculate_tau_rcd(h1=vivax2_rcd$hh[5000]*parameters2_rcd$rho, h2=vivax2_rcd$h[5000]-vivax2_rcd$hh[5000]*parameters2_rcd$rho,
                            alpha=parameters2_rcd$alpha,
                            rho=parameters2_rcd$rho,
                            iota=parameters2_rcd$iota, nu=parameters2_rcd$nu,
                            eta=parameters2_rcd$eta,
                            r=parameters2_rcd$r, rho2=parameters2_rcd$rho2)
  expect_equal(my_tau2,parameters2_rcd$tau, tolerance = 1e-09, label = "calculating tau with h2 and h1")

})



test_that("test solve RCD in delay model (referral)", {

  parameters_rcd=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.0155531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=0.6,"omega"=1,
                      "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0, "hh"=0, "hhl"=0,
                      "sigma"=1/15, "T0"=0, "Tl"=0,
                      "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters_rcd_0iota=parameters_rcd; parameters_rcd_0iota$iota=0
  parameters_rcd_0nu=parameters_rcd; parameters_rcd_0nu$nu=0
  parameters_rcd_0sigma=parameters_rcd; parameters_rcd_0sigma$sigma=0

  rcd_0iota=simulate_vivax_delay_ode(parameters=parameters_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)
  rcd_0nu=simulate_vivax_delay_ode(parameters=parameters_rcd_0nu , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)
  rcd_0sigma=simulate_vivax_delay_ode(parameters=parameters_rcd_0sigma , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)
  vivax_rcd=simulate_vivax_delay_ode(parameters=parameters_rcd , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)

  h1=calculate_h1_rcd(h=vivax_rcd$h[5000], alpha=parameters_rcd$alpha,
                      rho=parameters_rcd$rho,
                      iota=parameters_rcd$iota, nu=parameters_rcd$nu,
                      eta=parameters_rcd$eta, tau=parameters_rcd$tau,
                      r=parameters_rcd$r, rho2=parameters_rcd$rho2 )

  expect_equal(h1, vivax_rcd$hh[5000]*parameters_rcd$rho)

  mylambda_calc=solve_lambda_vivax_rcd_referral(h=vivax_rcd$h[5000],h1=h1,
                                                   r=parameters_rcd$r,f=parameters_rcd$f,
                                                   gamma=parameters_rcd$gamma,
                                                   alpha=parameters_rcd$alpha, beta=parameters_rcd$beta,
                                                   rho=parameters_rcd$rho, sigma=parameters_rcd$sigma,
                                                   p=vivax_rcd$p[5000],
                                                   omega=parameters_rcd$omega,
                                                   iota=parameters_rcd$iota, nu=parameters_rcd$nu,
                                                   eta=parameters_rcd$eta, tau=parameters_rcd$tau, return.all = T)

  h1_rcd_0iota=calculate_h1_rcd(h=rcd_0iota$h[5000], alpha=parameters_rcd_0iota$alpha,
                                rho=parameters_rcd_0iota$rho,
                                iota=parameters_rcd_0iota$iota, nu=parameters_rcd_0iota$nu,
                                eta=parameters_rcd_0iota$eta, tau=parameters_rcd_0iota$tau,
                                r=parameters_rcd_0iota$r, rho2=parameters_rcd_0iota$rho2 )
  mylambda_calc0=solve_lambda_vivax_rcd_referral(h=rcd_0iota$h[5000], h1=h1_rcd_0iota,
                                                    r=parameters_rcd_0iota$r,f=parameters_rcd_0iota$f,
                                                    gamma=parameters_rcd_0iota$gamma,
                                                    alpha=parameters_rcd_0iota$alpha, beta=parameters_rcd_0iota$beta,
                                                    rho=parameters_rcd_0iota$rho, sigma=parameters_rcd$sigma,
                                                    p=rcd_0iota$p[5000],
                                                    omega=parameters_rcd_0iota$omega,
                                                    iota=parameters_rcd_0iota$iota, nu=parameters_rcd_0iota$nu,
                                                    eta=parameters_rcd_0iota$eta, tau=parameters_rcd_0iota$tau, return.all = T)

  h1_rcd_0nu=calculate_h1_rcd(h=rcd_0nu$h[5000], alpha=parameters_rcd_0nu$alpha,
                              rho=parameters_rcd_0nu$rho,
                              iota=parameters_rcd_0nu$iota, nu=parameters_rcd_0nu$nu,
                              eta=parameters_rcd_0nu$eta, tau=parameters_rcd_0nu$tau,
                              r=parameters_rcd_0nu$r , rho2=parameters_rcd_0nu$rho2 )
  mylambda_calc0nu=solve_lambda_vivax_rcd_referral(h=rcd_0nu$h[5000],h1=h1_rcd_0nu,
                                                      r=parameters_rcd_0nu$r,f=parameters_rcd_0nu$f,
                                                      gamma=parameters_rcd_0nu$gamma,
                                                      alpha=parameters_rcd_0nu$alpha, beta=parameters_rcd_0nu$beta,
                                                      rho=parameters_rcd_0nu$rho, sigma=parameters_rcd_0nu$sigma,
                                                      p=rcd_0nu$p[5000],
                                                      omega=parameters_rcd_0nu$omega,
                                                      iota=parameters_rcd_0nu$iota, nu=parameters_rcd_0nu$nu,
                                                      eta=parameters_rcd_0nu$eta, tau=parameters_rcd_0nu$tau, return.all = T)


  expect_equal(mylambda_calc$lambda_real, parameters_rcd$lambda, tolerance = 1e-09, label = "lambda, rcd no delay no CM")
  expect_equal(mylambda_calc0$lambda_real, parameters_rcd_0iota$lambda, tolerance = 1e-09, label = "lambda, rcd=0 no delay no CM")
  expect_equal(mylambda_calc0nu$lambda_real, parameters_rcd_0nu$lambda, tolerance = 1e-09, label = "lambda, rcd=0 no delay no CM")

  eq_vivax_rcd=get_equilibrium_states_vivax_rcd_referral(I=vivax_rcd$I[5000],
                                                            lambda=parameters_rcd$lambda, sigma=parameters_rcd$sigma,
                                                            r=parameters_rcd$r, gamma=parameters_rcd$gamma, f=parameters_rcd$f,
                                                            alpha=parameters_rcd$alpha,  beta=parameters_rcd$beta, rho=parameters_rcd$rho,
                                                            delta=parameters_rcd$delta, omega=parameters_rcd$omega,
                                                            iota_star=mylambda_calc$iota_star, nu=parameters_rcd$nu,
                                                            tau=parameters_rcd$tau, eta=parameters_rcd$eta, rho2=parameters_rcd$rho2)
  eq_vivax_rcd_0iota=get_equilibrium_states_vivax_rcd_referral(I=rcd_0iota$I[5000],
                                                                  lambda=parameters_rcd_0iota$lambda, sigma=parameters_rcd_0iota$sigma,
                                                                  r=parameters_rcd_0iota$r, gamma=parameters_rcd_0iota$gamma, f=parameters_rcd_0iota$f,
                                                                  alpha=parameters_rcd_0iota$alpha,  beta=parameters_rcd_0iota$beta, rho=parameters_rcd_0iota$rho,
                                                                  delta=parameters_rcd_0iota$delta, omega=parameters_rcd_0iota$omega,
                                                                  iota_star=mylambda_calc0$iota_star, nu=parameters_rcd_0iota$nu,
                                                                  tau=parameters_rcd_0iota$tau, eta=parameters_rcd_0iota$eta, rho2=parameters_rcd_0iota$rho2)
  eq_vivax_rcd_0nu=get_equilibrium_states_vivax_rcd_referral(I=rcd_0nu$I[5000],
                                                                lambda=parameters_rcd_0nu$lambda, sigma=parameters_rcd_0nu$sigma,
                                                                r=parameters_rcd_0nu$r, gamma=parameters_rcd_0nu$gamma, f=parameters_rcd_0nu$f,
                                                                alpha=parameters_rcd_0nu$alpha,  beta=parameters_rcd_0nu$beta, rho=parameters_rcd_0nu$rho,
                                                                delta=parameters_rcd_0nu$delta, omega=parameters_rcd_0nu$omega,
                                                                iota_star=mylambda_calc0nu$iota_star, nu=parameters_rcd_0nu$nu,
                                                                tau=parameters_rcd_0nu$tau, eta=parameters_rcd_0nu$eta, rho2=parameters_rcd_0nu$rho2)
  expect_equal(eq_vivax_rcd, as.list(vivax_rcd[5000,names(eq_vivax_rcd)]), tolerance = 1e-09, label = "equilibrium values, rcd no delay")
  expect_equal(eq_vivax_rcd_0iota, as.list(rcd_0iota[5000,names(eq_vivax_rcd_0iota)]), tolerance = 1e-09, label = "equilibrium values, rcd=0, no delay")
  expect_equal(eq_vivax_rcd_0nu, as.list(rcd_0nu[5000,names(eq_vivax_rcd_0nu)]), tolerance = 1e-09, label = "equilibrium values, rcd=0, no delay")
  expect_equal(eq_vivax_rcd_0nu, eq_vivax_rcd_0iota, tolerance = 1e-09, label = "equilibrium values, rcd=0, no delay")

  my_tau=calculate_tau_rcd(h1=vivax_rcd$hh[5000]*parameters_rcd$rho, h2=vivax_rcd$h[5000]-vivax_rcd$hh[5000]*parameters_rcd$rho,
                           alpha=parameters_rcd$alpha,
                           rho=parameters_rcd$rho,
                           iota=parameters_rcd$iota, nu=parameters_rcd$nu,
                           eta=parameters_rcd$eta,
                           r=parameters_rcd$r, rho2=parameters_rcd$rho2)
  expect_equal(my_tau,parameters_rcd$tau, tolerance = 1e-09, label = "calculating tau with h2 and h1")

  parameters2_rcd=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.0155531,"delta"=0.01,
                       "alpha"=0.2, "beta"=0.7, "rho"=0.4,"omega"=0.9,
                       "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0, "hh"=0, "hhl"=0,
                       "sigma"=1/15, "T0"=0, "Tl"=0,
                       "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "rho2"=1)

  parameters2_rcd_0iota=parameters2_rcd; parameters2_rcd_0iota$iota=0
  parameters2_rcd_0sigma=parameters2_rcd; parameters2_rcd_0sigma$sigma=0
  parameters2_rcd_infsigma=parameters2_rcd; parameters2_rcd_infsigma$sigma=1/0.4

  rcd2_0iota=simulate_vivax_delay_ode(parameters=parameters2_rcd_0iota , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)
  vivax2_rcd=simulate_vivax_delay_ode(parameters=parameters2_rcd , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)

  h1_2=calculate_h1_rcd(h=vivax2_rcd$h[5000], alpha=parameters2_rcd$alpha,
                        rho=parameters2_rcd$rho,
                        iota=parameters2_rcd$iota, nu=parameters2_rcd$nu,
                        eta=parameters2_rcd$eta, tau=parameters2_rcd$tau,
                        r=parameters2_rcd$r, rho2=parameters2_rcd$rho2 )

  mylambda_calc2=solve_lambda_vivax_rcd_referral(h=vivax2_rcd$h[5000],h1=h1_2,
                                                    r=parameters2_rcd$r,f=parameters2_rcd$f,
                                                    gamma=parameters2_rcd$gamma,
                                                    alpha=parameters2_rcd$alpha, beta=parameters2_rcd$beta,
                                                    rho=parameters2_rcd$rho, sigma=parameters2_rcd$sigma,
                                                    p=vivax2_rcd$p[5000],
                                                    omega=parameters2_rcd$omega,
                                                    iota=parameters2_rcd$iota, nu=parameters2_rcd$nu,
                                                    eta=parameters2_rcd$eta, tau=parameters2_rcd$tau,
                                                    return.all = T)
  h1_0iota2=calculate_h1_rcd(h=rcd2_0iota$h[5000], alpha=parameters2_rcd_0iota$alpha,
                             rho=parameters2_rcd_0iota$rho,
                             iota=parameters2_rcd_0iota$iota, nu=parameters2_rcd_0iota$nu,
                             eta=parameters2_rcd_0iota$eta, tau=parameters2_rcd_0iota$tau,
                             r=parameters2_rcd_0iota$r , rho2=parameters2_rcd_0iota$rho2 )

  mylambda_calc0_2=solve_lambda_vivax_rcd_referral(h=rcd2_0iota$h[5000],h1=h1_0iota2,
                                                      r=parameters2_rcd_0iota$r,f=parameters2_rcd_0iota$f,
                                                      gamma=parameters2_rcd_0iota$gamma,
                                                      alpha=parameters2_rcd_0iota$alpha, beta=parameters2_rcd_0iota$beta,
                                                      rho=parameters2_rcd_0iota$rho, sigma=parameters2_rcd_0iota$sigma,
                                                      p=rcd2_0iota$p[5000],
                                                      omega=parameters2_rcd_0iota$omega,
                                                      iota=parameters2_rcd_0iota$iota, nu=parameters2_rcd_0iota$nu,
                                                      eta=parameters2_rcd_0iota$eta, tau=parameters2_rcd_0iota$tau,
                                                      return.all = T)


  expect_equal(mylambda_calc2$lambda_real, parameters2_rcd$lambda, tolerance = 1e-09, label = "lambda, rcd no delay with CM")
  expect_equal(mylambda_calc0_2$lambda_real, parameters2_rcd_0iota$lambda, tolerance = 1e-09, label = "lambda, rcd=0 no delay with CM")

  eq_vivax_rcd2=get_equilibrium_states_vivax_rcd_referral(I=vivax2_rcd$I[5000],
                                                             lambda=parameters2_rcd$lambda, sigma=parameters2_rcd$sigma,
                                                             r=parameters2_rcd$r, gamma=parameters2_rcd$gamma, f=parameters2_rcd$f,
                                                             alpha=parameters2_rcd$alpha,  beta=parameters2_rcd$beta, rho=parameters2_rcd$rho,
                                                             delta=parameters2_rcd$delta, omega=parameters2_rcd$omega,
                                                             iota_star=mylambda_calc2$iota_star, nu=parameters2_rcd$nu,
                                                             tau=parameters2_rcd$tau, eta=parameters2_rcd$eta, rho2=parameters2_rcd$rho2)
  eq_vivax_rcd2_0iota=get_equilibrium_states_vivax_rcd_referral(I=rcd2_0iota$I[5000],
                                                                   lambda=parameters2_rcd_0iota$lambda, sigma=parameters2_rcd_0iota$sigma,
                                                                   r=parameters2_rcd_0iota$r, gamma=parameters2_rcd_0iota$gamma, f=parameters2_rcd_0iota$f,
                                                                   alpha=parameters2_rcd_0iota$alpha,  beta=parameters2_rcd_0iota$beta, rho=parameters2_rcd_0iota$rho,
                                                                   delta=parameters2_rcd_0iota$delta, omega=parameters2_rcd_0iota$omega,
                                                                   iota_star=mylambda_calc0_2$iota_star, nu=parameters2_rcd_0iota$nu,
                                                                   tau=parameters2_rcd_0iota$tau, eta=parameters2_rcd_0iota$eta, rho2=parameters2_rcd_0iota$rho2)
  expect_equal(eq_vivax_rcd2, as.list(vivax2_rcd[5000,names(eq_vivax_rcd2)]), tolerance = 1e-09, label = "equilibrium values, rcd no delay")
  expect_equal(eq_vivax_rcd2_0iota, as.list(rcd2_0iota[5000,names(eq_vivax_rcd2_0iota)]), tolerance = 1e-09, label = "equilibrium values, rcd=0, no delay")

  my_tau2=calculate_tau_rcd(h1=vivax2_rcd$hh[5000]*parameters2_rcd$rho, h2=vivax2_rcd$h[5000]-vivax2_rcd$hh[5000]*parameters2_rcd$rho,
                            alpha=parameters2_rcd$alpha,
                            rho=parameters2_rcd$rho,
                            iota=parameters2_rcd$iota, nu=parameters2_rcd$nu,
                            eta=parameters2_rcd$eta,
                            r=parameters2_rcd$r, rho2=parameters2_rcd$rho2)
  expect_equal(my_tau2,parameters2_rcd$tau, tolerance = 1e-09, label = "calculating tau with h2 and h1")

})



