test_that("test calculation of R0/Rc on data, with custom CM", {

  mydata=data.frame(incidence=c(12, 45,3,12, 45,3),
                    prop_import=c(0,0,0, 0.1,0.3,0.02),
                    alpha=c(0,0,0,0.4,0.4,0.4),
                    beta=c(1,1,1,0.7,0.7,0.7),
                    sigma=c(1/15,1/15,1/15,1/10,1/10,1/10),
                    rho=c(1,1,1,0.4,0.4,0.4)) %>%
    dplyr::mutate(h=incidence_year2day(incidence))

  new_data_all =calculate_r0_rc_fromdata_delay(df=mydata, f=1/72, gamma=1/223, r=1/60, return.all = T)

  lambda_true=c(0.006757782, 0.006810128, 0.006743638, 0.010866174, 0.005895041, 0.012852846)
  R0_true=c(1.0027967, 1.0105644, 1.0006978, 1.6124467, 0.8747734, 1.9072518)
  Rc_true=c(1.0027967, 1.0105644, 1.0006978, 0.8143201, 0.4417793, 0.9632030)
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
                    sigma=c(1/15,1/15,1/15,1/10,1/10,1/10),
                    omega=c(1,1,1,0.8,0.8,0.8)) %>%
    dplyr::mutate(h=incidence_year2day(incidence))

  new_data_all =calculate_r0_rc_fromdata_delay(df=mydata, f=1/72, gamma=1/223, r=1/60, return.all = T)

  lambda_true=c(0.006757782, 0.006810128, 0.006743638, 0.013582718, 0.007368802, 0.016066058)
  R0_true=c(1.0027967, 1.0105644, 1.0006978, 2.015558, 1.093467, 2.384065)
  Rc_true=c(1.0027967, 1.0105644, 1.0006978, 0.8143201, 0.4417793, 0.9632030)
  expect_equal(new_data_all$lambda, lambda_true ,tolerance = 1e-06, label = "lambda")
  expect_equal(new_data_all$R0, R0_true ,tolerance = 1e-06, label = "R0")
  expect_equal(new_data_all$Rc, Rc_true , tolerance = 1e-06, label = "Rc")

})





####################################################################
test_that("test simulation of future scenarios", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(rho=c(0.18,0.13,0.08) ,
                  rho.old=rho,
                  beta.old=c(0.431,0.429,0.422),
                  sigma.old=c(1/15,1/15,1/15),
                  TQ_effect=c(0.59,0.605,0.619),
                  alpha.old=0.95*rho,
                  id=c(1,2,3))%>%
    dplyr::mutate(h=incidence_year2day(incidence),prop_import=0,
                  alpha=alpha.old, beta=beta.old, sigma=sigma.old)


  mydata2=calculate_r0_rc_fromdata_delay(df=mydata,   f=1/69, gamma=1/383, r=1/60, return.all = T)
  simul.PQ=simulate_from_data_delay(mydata2,
                                              f=1/69, gamma=1/383, r=1/60,
                                              maxtime=2000,year=T)
  names(simul.PQ)=c("time",paste0(names(simul.PQ)[! names(simul.PQ) %in% c("time","id")], "_PQ"), "id")

  simul.TQ=simulate_from_data_delay(mydata2 %>% dplyr::mutate(beta.new=TQ_effect),
                                              f=1/69, gamma=1/383, r=1/60,
                                              maxtime=2000,year=T)
  names(simul.TQ)=c("time",paste0(names(simul.TQ)[! names(simul.TQ) %in% c("time","id")], "_TQ"), "id")

  db_compare=simul.PQ %>% dplyr::left_join(simul.TQ) %>%
    dplyr::mutate(effect_size=-(h_TQ-h_PQ)/h_PQ,year=time/365+2020)

  true_TQ_I=c(0.01813114, 0.01702976, 0.01591025, 0.01489733, 0.01397706, 0.01313775,
              0.12763043, 0.12222864, 0.11782000, 0.11460526, 0.11222417, 0.11043990,
              0.51527342, 0.51093339, 0.51074055, 0.51073450, 0.51073431, 0.51073430)

  true_TQ_S0=c(0.9650693, 0.9674487, 0.9695882, 0.9715239, 0.9732826, 0.9748867,
               0.7618467, 0.7732983, 0.7814635, 0.7874176, 0.7918280, 0.7951331,
               0.1825342, 0.1924759, 0.1927769, 0.1927863, 0.1927866, 0.1927866)

  expect_equal(db_compare$I_PQ[db_compare$year==2025], mydata2$I,tolerance = 1e-07, label = "I PQ")
  expect_equal(db_compare$I_TQ, true_TQ_I,tolerance = 1e-07, label = "I TQ")
  expect_equal(db_compare$S0_TQ, true_TQ_S0,tolerance = 1e-07, label = "S0 TQ")
  expect_equal(db_compare$effect_size[db_compare$year==2025], c(0.259977861, 0.130279425, 0.008809147),tolerance = 1e-07, label = "effect size")
})


test_that("test simulation of future scenarios, with importation", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      sigma.old=c(1/15,1/15,1/15),
      TQ_effect=c(0.59,0.605,0.619),
      alpha.old=0.95*rho,
      prop_import=c(0.1,0.02,0.2),
      id=c(1,2,3))%>%
    dplyr::mutate(h=incidence_year2day(incidence),
                  alpha=alpha.old, beta=beta.old, sigma=sigma.old)

  mydata2=calculate_r0_rc_fromdata_delay(df=mydata,   f=1/69, gamma=1/383, r=1/60, return.all = T)
  simul.PQ=simulate_from_data_delay(mydata2,
                                              f=1/69, gamma=1/383, r=1/60,
                                              maxtime=2000,year=T)
  names(simul.PQ)=c("time",paste0(names(simul.PQ)[! names(simul.PQ) %in% c("time","id")], "_PQ"), "id")

  simul.TQ=simulate_from_data_delay(mydata2 %>% dplyr::mutate(beta.new=TQ_effect),
                                              f=1/69, gamma=1/383, r=1/60,
                                              maxtime=2000,year=T)
  names(simul.TQ)=c("time",paste0(names(simul.TQ)[! names(simul.TQ) %in% c("time","id")], "_TQ"), "id")

  db_compare=simul.PQ %>% dplyr::left_join(simul.TQ) %>%
    dplyr::mutate(effect_size=-(h_TQ-h_PQ)/h_PQ,year=time/365+2020)

  true_TQ_I=c(0.01813114, 0.01721811, 0.01657331 ,0.01617896, 0.01593754, 0.01578964,
              0.12763043, 0.12241630, 0.11844896, 0.11577284, 0.11394541, 0.11268692,
              0.51527342, 0.51159722, 0.51149278, 0.51149077, 0.51149073, 0.51149073)

  expect_equal(simul.PQ$incidence_PQ, rep(mydata$incidence, each=6),tolerance = 1e-07, label = "incidence PQ")
  expect_equal(db_compare$I_PQ[db_compare$year==2025], mydata2$I,tolerance = 1e-07, label = "I PQ")
  expect_equal(db_compare$I_TQ, true_TQ_I,tolerance = 1e-07, label = "I TQ")
  expect_equal(db_compare$effect_size[db_compare$year==2025], c(0.126696521, 0.114044062, 0.007341144),tolerance = 1e-07, label = "effect size")
})


test_that("test simulation of future scenarios, with importation and vector control", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      TQ_effect=c(0.59,0.605,0.619),
      sigma.old=c(1/15,1/15,1/15),
      alpha.old=0.95*rho,
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      prop_import=c(0.1,0.02,0.2),
      id=c(1,2,3),
      omega.old=c(0.7,0.7,0.7))%>%
    dplyr::mutate(h=incidence_year2day(incidence),
                  alpha=alpha.old, beta=beta.old, sigma=sigma.old, omega=omega.old)

  mydata2=calculate_r0_rc_fromdata_delay(df=mydata,   f=1/69, gamma=1/383, r=1/60, return.all = T)
  simul.PQ=simulate_from_data_delay(mydata2,
                                              f=1/69, gamma=1/383, r=1/60,
                                              maxtime=2000,year=T)
  names(simul.PQ)=c("time",paste0(names(simul.PQ)[! names(simul.PQ) %in% c("time","id")], "_PQ"), "id")

  simul.TQ=simulate_from_data_delay(mydata2 %>% dplyr::mutate(beta.new=TQ_effect, omega.new=omega.old*0.9),
                                              f=1/69, gamma=1/383, r=1/60,
                                              maxtime=2000,year=T)
  names(simul.TQ)=c("time",paste0(names(simul.TQ)[! names(simul.TQ) %in% c("time","id")], "_TQ"), "id")

  db_compare=simul.PQ %>% dplyr::left_join(simul.TQ) %>%
    dplyr::mutate(effect_size=-(h_TQ-h_PQ)/h_PQ,year=time/365+2020)

  expect_equal(simul.PQ$incidence_PQ, rep(mydata$incidence, each=6),tolerance = 1e-07, label = "I PQ")
  expect_equal(db_compare$I_PQ[db_compare$year==2025], mydata2$I,tolerance = 1e-07, label = "I PQ")
  expect_equal(db_compare$effect_size[db_compare$year==2025], c(0.24388972, 0.32497697, 0.03446925),tolerance = 1e-07, label = "effect size")
})



test_that("test simulation of future scenarios starting from initial condition, at equilibrium", {

  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      TQ_effect=c(0.59,0.605,0.619),
      alpha.old=0.95*rho,
      sigma.old=c(1/15,1/15,1/15),
      lambda=c(0.0057453432   ,0.0091275728 ,0.0167897768  ),
      I=c(0.018131142 ,0.127630432 ,0.515273425 ),
      delta=c(0.000035654059 ,5.4114206e-05,3.772780814e-03),
      id=c(1,2,3),
      omega.old=c(0.7,0.7,0.7))

  simul_no_rcd=simulate_from_data_delay(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=F, rcd=F)

  my_initial_state=simul_no_rcd[simul_no_rcd$time==2000,c("Il", "I0","Sl", "S0", "Tl", "T0","h", "hl","hh", "hhl", "id")]
  names(my_initial_state)=c("Il_init", "I0_init","Sl_init", "S0_init","Tl_init", "T0_init","h_init","hl_init","hh_init","hhl_init",  "id")
  simul_no_rcd_chain=simulate_from_data_delay(df=mydata, from_equilibrium = FALSE,
                                        initial_states = my_initial_state,
                                        f=1/69, gamma=1/383, r=1/60,  maxtime=2000,year=F, rcd=F)

  expect_equal(simul_no_rcd, simul_no_rcd_chain,tolerance = 1e-07, label = "chaining from equilibrium is the same as simulating from equilibrium")

})


test_that("test simulation of future scenarios starting from initial condition, not at equilibrium, with delay", {

  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      TQ_effect=c(0.59,0.605,0.619),
      alpha.old=0.95*rho,
      sigma.old=c(1/15,1/15,1/15),
      lambda=c(0.007  ,0.009  ,0.04 ),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3),
      omega.old=c(0.7,0.7,0.7))

  simul_no_rcd=simulate_from_data_delay(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=1000,year=F, rcd=F)

  my_initial_state=simul_no_rcd[simul_no_rcd$time==100,c("Il", "I0","Sl", "S0", "Tl", "T0","h", "hl", "hh", "hhl", "id")]
  names(my_initial_state)=c("Il_init", "I0_init","Sl_init", "S0_init","Tl_init", "T0_init","h_init", "hl_init", "hh_init", "hhl_init",  "id")
  simul_no_rcd_chain=simulate_from_data_delay(df=mydata, from_equilibrium = FALSE,
                                        initial_states = my_initial_state,
                                        f=1/69, gamma=1/383, r=1/60,  maxtime=900,year=F, rcd=F) %>%
    dplyr::mutate(time=time+100)


  expect_equal(simul_no_rcd %>% dplyr::filter(time>=100), simul_no_rcd_chain,tolerance = 1e-07,
               label = "chaining does not change the trajectory")


  simul_no_rcd_chain_newalpha=simulate_from_data_delay(df=mydata %>% dplyr::mutate(alpha.new=0.5), from_equilibrium = FALSE,
                                                 initial_states = my_initial_state,
                                                 f=1/69, gamma=1/383, r=1/60,  maxtime=900,year=F, rcd=F) %>%
    dplyr::mutate(time=time+100)

  expect_lt(simul_no_rcd_chain_newalpha$I[900], simul_no_rcd_chain$I[900],
            label = "adding an intervention changes the trajectory")


})

