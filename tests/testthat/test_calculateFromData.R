test_that("test calculation of R0/Rc on data, with default values", {

  mydata=data.frame(incidence=c(12, 45,3,12, 45,3),
                    prop_import=c(0,0,0, 0.1,0.3,0.02)) %>%
    dplyr::mutate(h=incidence_year2day(incidence))

  new_data =calculate_r0_rc_fromdata(df=mydata, f=1/72, gamma=1/223, r=1/60)
  new_data_all =calculate_r0_rc_fromdata(df=mydata, f=1/72, gamma=1/223, r=1/60, return.all = T)

  expect_equal(dim(new_data), c(6,10) , label = "dimension, return.all=F ")
  expect_equal(dim(new_data_all), c(6,12) , label = "dimension, return.all=T ")

  lambda_true=c(0.006757782, 0.006810128, 0.006743638, 0.005087821, 0.001772866, 0.006410140)
  R0_true=c(1.0027967, 1.0105644, 1.0006978, 0.7549889, 0.2630780, 0.9512096)
  Rc_true=c(1.0027967, 1.0105644, 1.0006978, 0.7549889, 0.2630780, 0.9512096)
  expect_equal(new_data_all$lambda, lambda_true ,tolerance = 1e-06, label = "lambda")
  expect_equal(new_data_all$R0, R0_true ,tolerance = 1e-06, label = "R0")
  expect_equal(new_data_all$Rc, Rc_true , tolerance = 1e-06, label = "Rc")

})


test_that("test calculation of R0/Rc on data, with custom CM", {

  mydata=data.frame(incidence=c(12, 45,3,12, 45,3),
                    prop_import=c(0,0,0, 0.1,0.3,0.02),
                    alpha=c(0,0,0,0.4,0.4,0.4),
                    beta=c(1,1,1,0.7,0.7,0.7),
                    rho=c(1,1,1,0.4,0.4,0.4)) %>%
    dplyr::mutate(h=incidence_year2day(incidence))

  new_data_all =calculate_r0_rc_fromdata(df=mydata, f=1/72, gamma=1/223, r=1/60, return.all = T)

  lambda_true=c(0.006757782, 0.006810128, 0.006743638, 0.012614889, 0.007171192, 0.014791913)
  R0_true=c(1.0027967, 1.0105644, 1.0006978, 1.871941, 1.064143, 2.194993)
  Rc_true=c(1.0027967, 1.0105644, 1.0006978, 0.8228744, 0.4677798, 0.9648826)
  expect_equal(new_data_all$lambda, lambda_true ,tolerance = 1e-06, label = "lambda")
  expect_equal(new_data_all$R0, R0_true ,tolerance = 1e-06, label = "R0")
  expect_equal(new_data_all$Rc, Rc_true , tolerance = 1e-06, label = "Rc")

})


test_that("test calculation of R0/Rc on data, with VC", {

  mydata=data.frame(incidence=c(12, 45,3,12, 45,3),
                    prop_import=c(0,0,0, 0.1,0.3,0.02),
                    alpha=c(0,0,0,0.4,0.4,0.4),
                    beta=c(1,1,1,0.7,0.7,0.7),
                    rho=c(1,1,1,0.4,0.4,0.4),
                    omega=c(1,1,1,0.8,0.8,0.8)) %>%
    dplyr::mutate(h=incidence_year2day(incidence))

  new_data_all =calculate_r0_rc_fromdata(df=mydata, f=1/72, gamma=1/223, r=1/60, return.all = T)

  lambda_true=c(0.006757782, 0.006810128, 0.006743638, 0.015768611, 0.008963990, 0.018489892)
  R0_true=c(1.0027967, 1.0105644, 1.0006978, 2.339926 , 1.330179 , 2.743741)
  Rc_true=c(1.0027967, 1.0105644, 1.0006978, 0.8228744, 0.4677798, 0.9648826)
  expect_equal(new_data_all$lambda, lambda_true ,tolerance = 1e-06, label = "lambda")
  expect_equal(new_data_all$R0, R0_true ,tolerance = 1e-06, label = "R0")
  expect_equal(new_data_all$Rc, Rc_true , tolerance = 1e-06, label = "Rc")

})


#############################################################################
test_that("test correct data frame dimension when adding uncertainty", {

  mydata=data.frame(cases=c(12, 45,3,12, 45,3),
                    population=c(10000,10000,10000,1500,1500,1500),
                    cases_local=c(10,3,1,10,3,1))

  data_uncertainty=sample_uncertainty_incidence_import(mydata, ndraw = 100)

  expect_equal(dim(data_uncertainty), c(600,5) , label = "dimension, return.all=F ")

})



####################################################################
test_that("test simulation of future scenarios", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(rho=c(0.18,0.13,0.08) ,
                  rho.old=rho ,
                  beta.old=c(0.431,0.429,0.422),
                  TQ_effect=c(0.59,0.605,0.619),
                  alpha.old=0.95*rho,
                  lambda=c(0.006428464 ,0.007131245 ,0.019061302),
                  I=c(0.01741279 ,0.12413235 ,0.50693425 ),
                  id=c(1,2,3))

  simul.PQ=simulate_from_data(mydata,
                                          f=1/69, gamma=1/383, r=1/60,
                                          maxtime=2000,year=T)
  names(simul.PQ)=c("time",paste0(names(simul.PQ)[! names(simul.PQ) %in% c("time","id")], "_PQ"), "id")

  simul.TQ=simulate_from_data(mydata %>% dplyr::mutate(beta.new=TQ_effect),
                                          f=1/69, gamma=1/383, r=1/60,
                                          maxtime=2000,year=T)
  names(simul.TQ)=c("time",paste0(names(simul.TQ)[! names(simul.TQ) %in% c("time","id")], "_TQ"), "id")

  db_compare=simul.PQ %>% dplyr::left_join(simul.TQ) %>%
    dplyr::mutate(effect_size=-(h_TQ-h_PQ)/h_PQ,year=time/365+2020)

  true_TQ_I=c(0.01741279, 0.01606616, 0.01471418, 0.01351327, 0.01244107, 0.01147929,
              0.12413235, 0.11741081, 0.11194112, 0.10796205, 0.10500922, 0.10278525,
              0.50693425, 0.50127738, 0.50099896, 0.50098923, 0.50098889, 0.50098887)

  true_TQ_S0=c(0.9660211, 0.9689859, 0.9715949, 0.9739123, 0.9759816, 0.9778378,
               0.7664995, 0.7809187, 0.7911016, 0.7985103, 0.8040091, 0.8081510,
               0.1908853, 0.2038417, 0.2042764, 0.2042916, 0.2042922, 0.2042922)

  expect_equal(db_compare$I_PQ[db_compare$year==2025], mydata$I,tolerance = 1e-07, label = "I PQ")
  expect_equal(db_compare$I_TQ, true_TQ_I,tolerance = 1e-07, label = "I TQ")
  expect_equal(db_compare$S0_TQ, true_TQ_S0,tolerance = 1e-07, label = "S0 TQ")
  expect_equal(db_compare$effect_size[db_compare$year==2025], c(0.32270753, 0.16637104, 0.01172809),tolerance = 1e-07, label = "effect size")
})


test_that("test simulation of future scenarios, with importation", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      TQ_effect=c(0.59,0.605,0.619),
      alpha.old=0.95*rho,
      lambda=c(0.004382382 ,0.004960252 ,0.015403064),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3))

  simul.PQ=simulate_from_data(mydata, f=1/69, gamma=1/383, r=1/60,
                                          maxtime=2000,year=T)
  names(simul.PQ)=c("time",paste0(names(simul.PQ)[! names(simul.PQ) %in% c("time","id")], "_PQ"), "id")

  simul.TQ=simulate_from_data(mydata %>% dplyr::mutate(beta.new=TQ_effect),
                                          f=1/69, gamma=1/383, r=1/60,
                                          maxtime=2000,year=T)
  names(simul.TQ)=c("time",paste0(names(simul.TQ)[! names(simul.TQ) %in% c("time","id")], "_TQ"), "id")

  db_compare=simul.PQ %>% dplyr::left_join(simul.TQ) %>%
    dplyr::mutate(effect_size=-(h_TQ-h_PQ)/h_PQ,year=time/365+2020)

  true_TQ_I=c(0.01741279, 0.01629352, 0.01551324, 0.01504596, 0.01476577, 0.01459762,
              0.12413235, 0.11846027, 0.11519186, 0.11357376, 0.11276799, 0.11236553,
              0.50693425, 0.50175677, 0.50155668, 0.50155132, 0.50155117, 0.50155117)

  true_TQ_S0=c(0.9660211, 0.9685962, 0.9701411, 0.9710662, 0.9716210, 0.9719539,
               0.7664995, 0.7792322, 0.7854942, 0.7885941, 0.7901379, 0.7909090,
               0.1908853, 0.2033691, 0.2036964, 0.2037052, 0.2037055, 0.2037055)

  expect_equal(simul.PQ$incidence_PQ, rep(mydata$incidence, each=6),tolerance = 1e-07, label = "incidence PQ")
  expect_equal(db_compare$I_PQ[db_compare$year==2025], mydata$I,tolerance = 1e-07, label = "I PQ")
  expect_equal(db_compare$I_TQ, true_TQ_I,tolerance = 1e-07, label = "I TQ")
  expect_equal(db_compare$S0_TQ, true_TQ_S0,tolerance = 1e-07, label = "S0 TQ")
  expect_equal(db_compare$effect_size[db_compare$year==2025], c(0.15884026, 0.09389020, 0.01061887),tolerance = 1e-07, label = "effect size")
})


test_that("test simulation of future scenarios, with importation and vector control", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      TQ_effect=c(0.59,0.605,0.619),
      alpha.old=0.95*rho,
      lambda=c(0.006260546  ,0.007086074  ,0.022004378 ),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3),
      omega.old=c(0.7,0.7,0.7))

  simul.PQ=simulate_from_data(mydata,
                                          f=1/69, gamma=1/383, r=1/60,
                                          maxtime=2000,year=T)
  names(simul.PQ)=c("time",paste0(names(simul.PQ)[! names(simul.PQ) %in% c("time","id")], "_PQ"), "id")

  simul.TQ=simulate_from_data(mydata %>% dplyr::mutate(beta.new=TQ_effect, omega.new=omega.old*0.9),
                                          f=1/69, gamma=1/383, r=1/60,
                                          maxtime=2000,year=T)
  names(simul.TQ)=c("time",paste0(names(simul.TQ)[! names(simul.TQ) %in% c("time","id")], "_TQ"), "id")

  db_compare=simul.PQ %>% dplyr::left_join(simul.TQ) %>%
    dplyr::mutate(effect_size=-(h_TQ-h_PQ)/h_PQ,year=time/365+2020)

  true_TQ_I=c(0.01741279, 0.01519767, 0.01391008, 0.01319473, 0.01279654, 0.01257465,
              0.12413235, 0.11117220, 0.10490999, 0.10190372, 0.10044583, 0.09973533,
              0.50693425, 0.48160057, 0.48080663, 0.48077663, 0.48077549, 0.48077545)

  true_TQ_S0=c(0.9660211, 0.9706025, 0.9731664, 0.9745907, 0.9753836, 0.9758254,
               0.7664995, 0.7918453, 0.8039268, 0.8097270, 0.8125400, 0.8139110,
               0.1908853, 0.2247011, 0.2260237, 0.2260737, 0.2260756, 0.2260756)

  expect_equal(simul.PQ$incidence_PQ, rep(mydata$incidence, each=6),tolerance = 1e-07, label = "I PQ")
  expect_equal(db_compare$I_PQ[db_compare$year==2025], mydata$I,tolerance = 1e-07, label = "I PQ")
  expect_equal(db_compare$I_TQ, true_TQ_I,tolerance = 1e-07, label = "I TQ")
  expect_equal(db_compare$S0_TQ, true_TQ_S0,tolerance = 1e-07, label = "S0 TQ")
  expect_equal(db_compare$effect_size[db_compare$year==2025], c(0.27418992, 0.19495853, 0.05160195),tolerance = 1e-07, label = "effect size")
})



test_that("test simulation of future scenarios starting from initial condition, at equilibrium", {

  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      TQ_effect=c(0.59,0.605,0.619),
      alpha.old=0.95*rho,
      lambda=c(0.006260546  ,0.007086074  ,0.022004378 ),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3),
      omega.old=c(0.7,0.7,0.7))

  simul_no_rcd=simulate_from_data(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=F, rcd=F)

  my_initial_state=simul_no_rcd[simul_no_rcd$time==2000,c("Il", "I0","Sl", "S0","h", "hr", "hl", "hh", "hhl", "id")]
  names(my_initial_state)=c("Il_init", "I0_init","Sl_init", "S0_init","h_init", "hr_init", "hl_init", "hh_init", "hhl_init", "id")
  simul_no_rcd_chain=simulate_from_data(df=mydata, from_equilibrium = FALSE,
                                        initial_states = my_initial_state,
                                        f=1/69, gamma=1/383, r=1/60,  maxtime=2000,year=F, rcd=F)

  expect_equal(simul_no_rcd, simul_no_rcd_chain,tolerance = 1e-07, label = "chaining from equilibrium is the same as simulating from equilibrium")

})


test_that("test simulation of future scenarios starting from initial condition, not at equilibrium", {

  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      TQ_effect=c(0.59,0.605,0.619),
      alpha.old=0.95*rho,
      lambda=c(0.007  ,0.009  ,0.04 ),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3),
      omega.old=c(0.7,0.7,0.7))

  simul_no_rcd=simulate_from_data(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=1000,year=F, rcd=F)

  my_initial_state=simul_no_rcd[simul_no_rcd$time==100,c("Il", "I0","Sl", "S0","h", "hr","hl","hh", "hhl", "id")]
  names(my_initial_state)=c("Il_init", "I0_init","Sl_init", "S0_init","h_init", "hr_init", "hl_init", "hh_init", "hhl_init", "id")
  simul_no_rcd_chain=simulate_from_data(df=mydata, from_equilibrium = FALSE,
                                        initial_states = my_initial_state,
                                        f=1/69, gamma=1/383, r=1/60,  maxtime=900,year=F, rcd=F) %>%
    dplyr::mutate(time=time+100)


  expect_equal(simul_no_rcd %>% dplyr::filter(time>=100), simul_no_rcd_chain,tolerance = 1e-07,
               label = "chaining does not change the trajectory")


  simul_no_rcd_chain_newalpha=simulate_from_data(df=mydata %>% dplyr::mutate(alpha.new=0.5), from_equilibrium = FALSE,
                                        initial_states = my_initial_state,
                                        f=1/69, gamma=1/383, r=1/60,  maxtime=900,year=F, rcd=F) %>%
    dplyr::mutate(time=time+100)

  expect_lt(simul_no_rcd_chain_newalpha$I[900], simul_no_rcd_chain$I[900],
               label = "cadding an intervention changes the trajectory")


})

