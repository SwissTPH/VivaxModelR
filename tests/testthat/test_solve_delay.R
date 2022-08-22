test_that("test solve model with CM", {

  parameters=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.03,"delta"=0,
                  "alpha"=0.4, "beta"=0.7, "sigma"=1/15, "rho"=0.4,"omega"=1,
                  "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0,"hh"=0, "hhl"=0, "T0"=0, "Tl"=0)

  simul=simulate_vivax_delay_ode(parameters, ODEmodel =ode_vivax_delay , maxtime = 500000)

  true_R0=get_r0_vivax(lambda=parameters$lambda, f=parameters$f,
                       r=parameters$r, gamma=parameters$gamma)
  true_Rc=get_rc_vivax_delay(lambda=parameters$lambda, f=parameters$f,
                       r=parameters$r, gamma=parameters$gamma,
                       alpha=parameters$alpha, beta=parameters$beta, sigma=parameters$sigma)

  simul=simul%>%
    dplyr::mutate(I_calc=h*(parameters$r+(1-parameters$alpha)*parameters$sigma)/parameters$r/parameters$rho/(parameters$r+parameters$sigma),
                  delta_calc=p * h/(1-I_calc)/parameters$rho)

  lambda_calc=solve_lambda_vivax_delay(h=simul$h[500000],
                                 r=parameters$r,f=parameters$f,
                                 gamma=parameters$gamma,
                                 alpha=parameters$alpha, beta=parameters$beta, sigma=parameters$sigma,
                                 rho=parameters$rho, p=simul$p[500000],
                                 omega=parameters$omega)

  eq_theo=get_equilibrium_states_vivax_delay(simul$I[500000], parameters$lambda, parameters$r, parameters$gamma, parameters$f, parameters$alpha, parameters$beta, parameters$sigma, parameters$rho, parameters$delta, parameters$omega)

  expect_equal(simul$I_calc[500000], simul$I[500000], tolerance = 1e-09, label = "I")
  expect_equal(simul$delta_calc[500000], parameters$delta, tolerance = 1e-09, label = "delta")
  expect_equal(eq_theo, as.list(simul[500000,names(eq_theo)]), tolerance = 1e-09, label = "equilibrium values")
  expect_equal(lambda_calc, parameters$lambda, tolerance = 1e-09, label = "lambda")
  expect_equal(true_R0, get_r0_vivax(lambda_calc, parameters$f, parameters$r, parameters$gamma), tolerance = 1e-09, label = "R0")
  expect_equal(true_Rc, get_rc_vivax_delay(lambda_calc, parameters$f, parameters$r, parameters$gamma, parameters$alpha, parameters$beta, parameters$sigma), tolerance = 1e-09, label = "Rc")

})


test_that("test solve model with CM and importation", {

  parameters=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.02,"delta"=0.01,
                  "alpha"=0.4, "beta"=0.7, "sigma"=1/15, "rho"=0.4,"omega"=1,
                  "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0,"hh"=0, "hhl"=0, "T0"=0, "Tl"=0)

  simul=simulate_vivax_delay_ode(parameters, ODEmodel =ode_vivax_delay , maxtime = 500000)

  true_R0=get_r0_vivax(lambda=parameters$lambda, f=parameters$f,
                       r=parameters$r, gamma=parameters$gamma)
  true_Rc=get_rc_vivax_delay(lambda=parameters$lambda, f=parameters$f,
                       r=parameters$r, gamma=parameters$gamma,
                       alpha=parameters$alpha, beta=parameters$beta, sigma=parameters$sigma)

  simul=simul%>%
    dplyr::mutate(I_calc=h*(parameters$r+(1-parameters$alpha)*parameters$sigma)/parameters$r/parameters$rho/(parameters$r+parameters$sigma),
                  delta_calc=p * h/(1-I_calc)/parameters$rho)

  lambda_calc=solve_lambda_vivax_delay(h=simul$h[500000],
                                       r=parameters$r,f=parameters$f,
                                       gamma=parameters$gamma,
                                       alpha=parameters$alpha, beta=parameters$beta, sigma=parameters$sigma,
                                       rho=parameters$rho, p=simul$p[500000],
                                       omega=parameters$omega)

  eq_theo=get_equilibrium_states_vivax_delay(simul$I[500000], parameters$lambda, parameters$r, parameters$gamma, parameters$f, parameters$alpha, parameters$beta, parameters$sigma, parameters$rho, parameters$delta, parameters$omega)

  expect_equal(simul$I_calc[500000], simul$I[500000], tolerance = 1e-09, label = "I")
  expect_equal(simul$delta_calc[500000], parameters$delta, tolerance = 1e-09, label = "delta")
  expect_equal(eq_theo, as.list(simul[500000,names(eq_theo)]), tolerance = 1e-09, label = "equilibrium values")
  expect_equal(lambda_calc, parameters$lambda, tolerance = 1e-09, label = "lambda")
  expect_equal(true_R0, get_r0_vivax(lambda_calc, parameters$f, parameters$r, parameters$gamma), tolerance = 1e-09, label = "R0")
  expect_equal(true_Rc, get_rc_vivax_delay(lambda_calc, parameters$f, parameters$r, parameters$gamma, parameters$alpha, parameters$beta, parameters$sigma), tolerance = 1e-09, label = "Rc")

})


test_that("test solve model with CM, importation and VC", {

  parameters=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.02,"delta"=0.01,
                  "alpha"=0.4, "beta"=0.7, "sigma"=1/15, "rho"=0.4,"omega"=0.7,
                  "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0,"hh"=0, "hhl"=0, "T0"=0, "Tl"=0)

  simul=simulate_vivax_delay_ode(parameters, ODEmodel =ode_vivax_delay , maxtime = 500000)

  true_R0=get_r0_vivax(lambda=parameters$lambda, f=parameters$f,
                       r=parameters$r, gamma=parameters$gamma)
  true_Rc=get_rc_vivax_delay(lambda=parameters$lambda, f=parameters$f,
                       r=parameters$r, gamma=parameters$gamma,
                       alpha=parameters$alpha, beta=parameters$beta, sigma=parameters$sigma)

  simul=simul%>%
    dplyr::mutate(I_calc=h*(parameters$r+(1-parameters$alpha)*parameters$sigma)/parameters$r/parameters$rho/(parameters$r+parameters$sigma),
                  delta_calc=p * h/(1-I_calc)/parameters$rho)

  eq_theo=get_equilibrium_states_vivax_delay(simul$I[500000], parameters$lambda, parameters$r, parameters$gamma, parameters$f, parameters$alpha, parameters$beta, parameters$sigma, parameters$rho, parameters$delta, parameters$omega)

  lambda_calc=solve_lambda_vivax_delay(h=simul$h[500000],
                                       r=parameters$r,f=parameters$f,
                                       gamma=parameters$gamma,
                                       alpha=parameters$alpha, beta=parameters$beta, sigma=parameters$sigma,
                                       rho=parameters$rho, p=simul$p[500000],
                                       omega=parameters$omega)

  eq_theo=get_equilibrium_states_vivax_delay(simul$I[500000], parameters$lambda, parameters$r, parameters$gamma, parameters$f, parameters$alpha, parameters$beta, parameters$sigma, parameters$rho, parameters$delta, parameters$omega)

  expect_equal(simul$I_calc[500000], simul$I[500000], tolerance = 1e-09, label = "I")
  expect_equal(simul$delta_calc[500000], parameters$delta, tolerance = 1e-09, label = "delta")
  expect_equal(eq_theo, as.list(simul[500000,names(eq_theo)]), tolerance = 1e-09, label = "equilibrium values")
  expect_equal(lambda_calc, parameters$lambda, tolerance = 1e-09, label = "lambda")
  expect_equal(true_R0, get_r0_vivax(lambda_calc, parameters$f, parameters$r, parameters$gamma), tolerance = 1e-09, label = "R0")
  expect_equal(true_Rc, get_rc_vivax_delay(lambda_calc, parameters$f, parameters$r, parameters$gamma, parameters$alpha, parameters$beta, parameters$sigma), tolerance = 1e-09, label = "Rc")

})

test_that("test compare delay model when sigma=0 with non delay model", {

  parameters_d1=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.03,"delta"=0,
                  "alpha"=0.4, "beta"=0.7, "sigma"=0, "rho"=0.4,"omega"=1,
                  "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0,"hh"=0, "hhl"=0, "T0"=0, "Tl"=0)

  simul_delay1=simulate_vivax_delay_ode(parameters_d1, ODEmodel =ode_vivax_delay , maxtime = 5000)

  parameters_d2=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.03,"delta"=0,
                  "alpha"=0, "beta"=0.7, "sigma"=1/15, "rho"=0.4,"omega"=1,
                  "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0,"hh"=0, "hhl"=0, "T0"=0, "Tl"=0)

  simul_delay2=simulate_vivax_delay_ode(parameters_d2, ODEmodel =ode_vivax_delay , maxtime = 5000)

  parameters_d3=list("r"=1/60, "gamma"=1/223,
                     "f"=1/72,"lambda"=0.4,"delta"=0,
                     "alpha"=1, "beta"=0, "sigma"=1/15, "rho"=0.4,"omega"=1,
                     "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0,"hh"=0, "hhl"=0, "T0"=0, "Tl"=0)

  simul_delay3=simulate_vivax_delay_ode(parameters_d3, ODEmodel =ode_vivax_delay , maxtime = 5000)


  parameters=list("r"=1/60, "gamma"=1/223,
                   "f"=1/72,"lambda"=0.03,"delta"=0,
                   "alpha"=0, "beta"=0.7,  "rho"=0.4,"omega"=1,
                   "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0)
  simul=simulate_vivax_ode(parameters, ODEmodel =ode_vivax_cm_importation_vc , maxtime = 5000)


  parameters_3=list("r"=1/60+1/15, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.4,"delta"=0,
                  "alpha"=0, "beta"=0,  "rho"=0.4,"omega"=1,
                  "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0)
  simul3=simulate_vivax_ode(parameters_3, ODEmodel =ode_vivax_cm_importation_vc , maxtime = 5000)

  expect_equal(simul_delay1[5000,c("S0", "Sl", "I" ,"h")], simul[5000,c("S0", "Sl", "I" ,"h")], tolerance = 1e-09, label = "delay sigma=0 and non delay alpha=0")
  expect_equal(simul_delay1[5000,c("S0", "Sl", "I" ,"h")], simul_delay2[5000,c("S0", "Sl", "I" ,"h")], tolerance = 1e-09, label = "delay alpha=0 and non delay alpha=0")
  expect_equal(simul_delay3[5000,c("S0", "Sl", "I" ,"h")], simul3[5000,c("S0", "Sl", "I" ,"h")], tolerance = 1e-09, label = "delay alpha=1 and non delay alpha=0 and r=r+sigma, beta=1")

})



test_that("test compare delay model when sigma=0 with non delay model, with importation and vector control", {

  parameters_d1=list("r"=1/60, "gamma"=1/223,
                     "f"=1/72,"lambda"=0.03,"delta"=0.1,
                     "alpha"=0.4, "beta"=0.7, "sigma"=0, "rho"=0.4,"omega"=0.7,
                     "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0,"hh"=0, "hhl"=0, "T0"=0, "Tl"=0)

  simul_delay1=simulate_vivax_delay_ode(parameters_d1, ODEmodel =ode_vivax_delay , maxtime = 5000)

  parameters_d2=list("r"=1/60, "gamma"=1/223,
                     "f"=1/72,"lambda"=0.03,"delta"=0.1,
                     "alpha"=0, "beta"=0.7, "sigma"=1/15, "rho"=0.4,"omega"=0.7,
                     "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0,"hh"=0, "hhl"=0, "T0"=0, "Tl"=0)

  simul_delay2=simulate_vivax_delay_ode(parameters_d2, ODEmodel =ode_vivax_delay , maxtime = 5000)

  parameters_d3=list("r"=1/60, "gamma"=1/223,
                     "f"=1/72,"lambda"=0.4,"delta"=0,
                     "alpha"=1, "beta"=0, "sigma"=1/15, "rho"=0.4,"omega"=1,
                     "U0"=0, "S0"=0.9, "Sl"=0, "Ul"=0.1, "h"=0, "hl"=0,"hh"=0, "hhl"=0, "T0"=0, "Tl"=0)

  simul_delay3=simulate_vivax_delay_ode(parameters_d3, ODEmodel =ode_vivax_delay , maxtime = 5000)


  parameters=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.03,"delta"=0.1,
                  "alpha"=0, "beta"=0.7,  "rho"=0.4,"omega"=0.7,
                  "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0)
  simul=simulate_vivax_ode(parameters, ODEmodel =ode_vivax_cm_importation_vc , maxtime = 5000)


  parameters_3=list("r"=1/60+1/15, "gamma"=1/223,
                    "f"=1/72,"lambda"=0.4,"delta"=0,
                    "alpha"=0, "beta"=0,  "rho"=0.4,"omega"=1,
                    "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0)
  simul3=simulate_vivax_ode(parameters_3, ODEmodel =ode_vivax_cm_importation_vc , maxtime = 5000)

  expect_equal(simul_delay1[5000,c("S0", "Sl", "I" ,"h")], simul[5000,c("S0", "Sl", "I" ,"h")], tolerance = 1e-09, label = "delay sigma=0 and non delay alpha=0")
  expect_equal(simul_delay1[5000,c("S0", "Sl", "I" ,"h")], simul_delay2[5000,c("S0", "Sl", "I" ,"h")], tolerance = 1e-09, label = "delay alpha=0 and non delay alpha=0")
  expect_equal(simul_delay3[5000,c("S0", "Sl", "I" ,"h")], simul3[5000,c("S0", "Sl", "I" ,"h")], tolerance = 1e-09, label = "delay alpha=1 and non delay alpha=0 and r=r+sigma, beta=1")

})
