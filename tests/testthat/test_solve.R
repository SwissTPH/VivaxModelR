test_that("test solve model with CM", {

  parameters=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.0155531,"delta"=0,
                  "alpha"=0.4, "beta"=0.7, "rho"=0.4,"omega"=1,
                  "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0, "hh"=0, "hhl"=0)

  simul=simulate_vivax_ode(parameters, ODEmodel =ode_vivax_cm_importation_vc , maxtime = 500000)

  true_R0=get_r0_vivax(lambda=parameters$lambda, f=parameters$f,
                       r=parameters$r, gamma=parameters$gamma)
  true_Rc=get_rc_vivax(lambda=parameters$lambda, f=parameters$f,
                       r=parameters$r, gamma=parameters$gamma,
                       alpha=parameters$alpha, beta=parameters$beta)

  simul=simul%>%
    dplyr::mutate(I_calc=h*(1-parameters$alpha)/parameters$r/parameters$rho,
           delta_calc=p * h/(1-I_calc)/parameters$rho)

  lambda_calc=solve_lambda_vivax(h=simul$h[500000],
                                                   r=parameters$r,f=parameters$f,
                                                   gamma=parameters$gamma,
                                                   alpha=parameters$alpha, beta=parameters$beta,
                                                   rho=parameters$rho, p=simul$p[500000],
                                                   omega=parameters$omega)
  hr_calc=get_prop_relapse(lambda=lambda_calc, h=simul$h[500000], r=parameters$r,  f=parameters$f,
                                 gamma=parameters$gamma, alpha=parameters$alpha, beta=parameters$beta,
                                 rho=parameters$rho, p=simul$p[500000], omega=parameters$omega)

  eq_theo=get_equilibrium_states_vivax(simul$I[500000], parameters$lambda, parameters$r, parameters$gamma, parameters$f, parameters$alpha, parameters$beta, parameters$rho, parameters$delta, parameters$omega)

  eq_calc=get_equilibrium_states_vivax(simul$I[500000], lambda_calc, parameters$r, parameters$gamma, parameters$f, parameters$alpha, parameters$beta, parameters$rho, parameters$delta, parameters$omega)


  expect_equal(simul$I_calc[500000], simul$I[500000], tolerance = 1e-09, label = "I")
  expect_equal(simul$delta_calc[500000], parameters$delta, tolerance = 1e-09, label = "delta")
  expect_equal(eq_theo, as.list(simul[500000,names(eq_theo)]), tolerance = 1e-09, label = "equilibrium values")
  expect_equal(lambda_calc, parameters$lambda, tolerance = 1e-09, label = "lambda")
  expect_equal(true_R0, get_r0_vivax(lambda_calc, parameters$f, parameters$r, parameters$gamma), tolerance = 1e-09, label = "R0")
  expect_equal(true_Rc, get_rc_vivax(lambda_calc, parameters$f, parameters$r, parameters$gamma, parameters$alpha, parameters$beta), tolerance = 1e-09, label = "Rc")
  expect_equal(hr_calc, simul$hr[500000]/simul$h[500000], tolerance = 1e-09, label = "relapse prop")

})

test_that("test solve model with CM and importation", {

  parameters=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.02,"delta"=0.01,
                  "alpha"=0.4, "beta"=0.7, "rho"=0.4,"omega"=1,
                  "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0, "hh"=0, "hhl"=0)

  simul=simulate_vivax_ode(parameters, ODEmodel =ode_vivax_cm_importation_vc , maxtime = 500000)

  true_R0=get_r0_vivax(lambda=parameters$lambda, f=parameters$f,
                       r=parameters$r, gamma=parameters$gamma)
  true_Rc=get_rc_vivax(lambda=parameters$lambda, f=parameters$f,
                          r=parameters$r, gamma=parameters$gamma,
                          alpha=parameters$alpha, beta=parameters$beta)

  simul=simul%>%
    dplyr::mutate(I_calc=h*(1-parameters$alpha)/parameters$r/parameters$rho,
                  delta_calc=p * h/(1-I_calc)/parameters$rho)

  lambda_calc=solve_lambda_vivax(h=simul$h[500000],
                                                   r=parameters$r,f=parameters$f,
                                                   gamma=parameters$gamma,
                                                   alpha=parameters$alpha, beta=parameters$beta,
                                                   rho=parameters$rho, p=simul$p[500000],
                                                   omega=parameters$omega)
  hr_calc=get_prop_relapse(lambda=lambda_calc,
                                 h=simul$h[500000], r=parameters$r,  f=parameters$f,
                                 gamma=parameters$gamma, alpha=parameters$alpha, beta=parameters$beta,
                                 rho=parameters$rho, p=simul$p[500000], omega=parameters$omega)

  eq_theo=get_equilibrium_states_vivax(simul$I[500000], parameters$lambda, parameters$r, parameters$gamma, parameters$f, parameters$alpha, parameters$beta, parameters$rho, parameters$delta, parameters$omega)
  eq_calc=get_equilibrium_states_vivax(simul$I[500000], lambda_calc, parameters$r, parameters$gamma, parameters$f, parameters$alpha, parameters$beta, parameters$rho, parameters$delta, parameters$omega)


  expect_equal(simul$I_calc[500000], simul$I[500000], tolerance = 1e-09, label = "I")
  expect_equal(simul$delta_calc[500000], parameters$delta, tolerance = 1e-09, label = "delta")
  expect_equal(eq_theo, as.list(simul[500000,names(eq_theo)]), tolerance = 1e-09, label = "equilibrium values")
  expect_equal(lambda_calc, parameters$lambda, tolerance = 1e-09, label = "lambda")
  expect_equal(true_R0, get_r0_vivax(lambda_calc, parameters$f, parameters$r, parameters$gamma), tolerance = 1e-09, label = "R0")
  expect_equal(true_Rc, get_rc_vivax(lambda_calc, parameters$f, parameters$r, parameters$gamma, parameters$alpha, parameters$beta), tolerance = 1e-09, label = "Rc")
  expect_equal(hr_calc, simul$hr[500000]/simul$h[500000], tolerance = 1e-09, label = "relapse prop")
  expect_equal(eq_calc, as.list(simul[500000,names(eq_theo)]), tolerance = 1e-09, label = "equilibrium values, with estimated lambda")

})


test_that("test solve model with CM, importation and vector control", {

  parameters=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.02,"delta"=0.01,
                  "alpha"=0.4, "beta"=0.7, "rho"=0.4,"omega"=0.7,
                  "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0, "hh"=0, "hhl"=0)

  simul=simulate_vivax_ode(parameters, ODEmodel =ode_vivax_cm_importation_vc , maxtime = 500000)

  true_R0=get_r0_vivax(lambda=parameters$lambda, f=parameters$f,
                       r=parameters$r, gamma=parameters$gamma)
  true_Rc=get_rc_vivax(lambda=parameters$lambda, f=parameters$f,
                          r=parameters$r, gamma=parameters$gamma,
                          alpha=parameters$alpha, beta=parameters$beta,
                          omega = parameters$omega)

  simul=simul%>%
    dplyr::mutate(I_calc=h*(1-parameters$alpha)/parameters$r/parameters$rho,
                  delta_calc=p * h/(1-I_calc)/parameters$rho)

  lambda_calc=solve_lambda_vivax(h=simul$h[500000],
                                                   r=parameters$r,f=parameters$f,
                                                   gamma=parameters$gamma,
                                                   alpha=parameters$alpha, beta=parameters$beta,
                                                   rho=parameters$rho, p=simul$p[500000],
                                                   omega = parameters$omega)
  hr_calc=get_prop_relapse(lambda=lambda_calc,
                                 h=simul$h[500000], r=parameters$r,  f=parameters$f,
                                 gamma=parameters$gamma, alpha=parameters$alpha, beta=parameters$beta,
                                 rho=parameters$rho, p=simul$p[500000], omega=parameters$omega)

  eq_theo=get_equilibrium_states_vivax(simul$I[500000], parameters$lambda, parameters$r, parameters$gamma, parameters$f, parameters$alpha, parameters$beta, parameters$rho, parameters$delta, parameters$omega)
  eq_calc=get_equilibrium_states_vivax(simul$I[500000], lambda_calc, parameters$r, parameters$gamma, parameters$f, parameters$alpha, parameters$beta, parameters$rho, parameters$delta, parameters$omega)


  expect_equal(simul$I_calc[500000], simul$I[500000], tolerance = 1e-09, label = "I")
  expect_equal(simul$delta_calc[500000], parameters$delta, tolerance = 1e-09, label = "delta")
  expect_equal(eq_theo, as.list(simul[500000,names(eq_theo)]), tolerance = 1e-09, label = "equilibrium values")
  expect_equal(lambda_calc, parameters$lambda, tolerance = 1e-09, label = "lambda")
  expect_equal(true_R0, get_r0_vivax(lambda_calc, parameters$f, parameters$r, parameters$gamma), tolerance = 1e-09, label = "R0")
  expect_equal(true_Rc, get_rc_vivax(lambda_calc, parameters$f, parameters$r, parameters$gamma, parameters$alpha, parameters$beta, omega = parameters$omega), tolerance = 1e-09, label = "Rc")
  expect_equal(hr_calc, simul$hr[500000]/simul$h[500000], tolerance = 1e-09, label = "relapse prop")
  expect_equal(eq_calc, as.list(simul[500000,names(eq_theo)]), tolerance = 1e-09, label = "equilibrium values, with estimated lambda")

})
