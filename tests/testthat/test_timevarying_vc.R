test_that("test time varying vector control working", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      alpha.old=0.95*rho,
      lambda=c(0.006260546  ,0.007086074  ,0.022004378 ),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3))

  my_omega=data.frame(t=seq(0, 70*365)) %>%
    dplyr::mutate(vc=ifelse(t<10*365, 0.7, 1-0.3*exp(-(t-10*365)/365/1.5)))
  omega_t <- approxfun(my_omega)

  mydata_tv=mydata
  #mydata_tv$omega.old=c(omega_t(0), omega_t(0), omega_t(0))
  mydata_tv$omega.old=c(omega_t, omega_t, omega_t)
  mydata_tv$omega.new=c(omega_t, omega_t, omega_t)
  simul.tv=simulate_from_data(mydata_tv,  f=1/69, gamma=1/383, r=1/60,  maxtime=60*365,year=F)

  mydata_1=mydata
  mydata_1$omega.old=c(1, 1, 1)
  simul.1=simulate_from_data(mydata_1,  f=1/69, gamma=1/383, r=1/60,  maxtime=60*365,year=F)

  mydata_07=mydata
  mydata_07$omega.old=c(0.7, 0.7, 0.7)
  simul.07=simulate_from_data(mydata_07,  f=1/69, gamma=1/383, r=1/60,  maxtime=60*365,year=F)

  expect_equal(simul.tv[simul.tv$time<10*365,], simul.07[simul.07$time<10*365,],tolerance = 1e-07, label = "before time varying")
  expect_equal(simul.tv[simul.tv$time>58*365,], simul.1[simul.1$time>58*365,],tolerance = 1e-05, label = "after time varying")
  expect_lt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==1],simul.1$I[simul.1$time==5000  & simul.tv$id==1], label = "between time varying, lower than 1")
  expect_lt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==2],simul.1$I[simul.1$time==5000 & simul.tv$id==2], label = "between time varying, lower than 1")
  expect_lt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==3],simul.1$I[simul.1$time==5000 & simul.tv$id==3], label = "between time varying, lower than 1")
  expect_gt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==1],simul.07$I[simul.07$time==5000  & simul.07$id==1], label = "between time varying, higher than 0.7")
  expect_gt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==2],simul.07$I[simul.07$time==5000 & simul.07$id==2], label = "between time varying, higher than 0.7")
  expect_gt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==3],simul.07$I[simul.07$time==5000 & simul.07$id==3], label = "between time varying, higher than 0.7")
})


test_that("test time varying vector control working, with delays", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      alpha.old=0.95*rho,
      sigma.old=c(1/15,1/15,1/15),
      lambda=c(0.006260546  ,0.007086074  ,0.022004378 ),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3))

  my_omega=data.frame(t=seq(0, 70*365)) %>%
    dplyr::mutate(vc=ifelse(t<10*365, 0.7, 1-0.3*exp(-(t-10*365)/365/1.5)))
  omega_t <- approxfun(my_omega)

  mydata_tv=mydata
  #mydata_tv$omega.old=c(omega_t(0), omega_t(0), omega_t(0))
  mydata_tv$omega.old=c(omega_t, omega_t, omega_t)
  mydata_tv$omega.new=c(omega_t, omega_t, omega_t)
  simul.tv=simulate_from_data_delay(mydata_tv,  f=1/69, gamma=1/383, r=1/60,  maxtime=60*365,year=F)

  mydata_1=mydata
  mydata_1$omega.old=c(1, 1, 1)
  simul.1=simulate_from_data_delay(mydata_1,  f=1/69, gamma=1/383, r=1/60,  maxtime=60*365,year=F)

  mydata_07=mydata
  mydata_07$omega.old=c(0.7, 0.7, 0.7)
  simul.07=simulate_from_data_delay(mydata_07,  f=1/69, gamma=1/383, r=1/60,  maxtime=60*365,year=F)

  expect_equal(simul.tv[simul.tv$time<10*365,], simul.07[simul.07$time<10*365,],tolerance = 1e-07, label = "before time varying")
  expect_equal(simul.tv[simul.tv$time>58*365,], simul.1[simul.1$time>58*365,],tolerance = 1e-05, label = "after time varying")
  expect_lt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==1],simul.1$I[simul.1$time==5000  & simul.tv$id==1], label = "between time varying, lower than 1")
  expect_lt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==2],simul.1$I[simul.1$time==5000 & simul.tv$id==2], label = "between time varying, lower than 1")
  expect_lt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==3],simul.1$I[simul.1$time==5000 & simul.tv$id==3], label = "between time varying, lower than 1")
  expect_gt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==1],simul.07$I[simul.1$time==5000  & simul.07$id==1], label = "between time varying, higher than 0.7")
  expect_gt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==2],simul.07$I[simul.1$time==5000 & simul.07$id==2], label = "between time varying, higher than 0.7")
  expect_gt(simul.tv$I[simul.tv$time==5000 & simul.tv$id==3],simul.07$I[simul.1$time==5000 & simul.07$id==3], label = "between time varying, higher than 0.7")
})


test_that("test simulation from equilibrium and chained, with time varying vector control", {

  # simulate a vector control intervention with 3 years half life, distributed every 4 years
  my_omega=vector_control_exponential_decay(initial_omega = 0.7, half_life = 3, every_x_years = 4, maxtime = 70*365 )

  mydata=data.frame(incidence=c(23,112),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$prop_import=c(0,0)
  mydata$omega=c(0.7,0.7)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "rho.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=NA, "beta.new"=NA, "omega.new"=my_omega, "rho.new"=NA)

  my_intervention_list=list(int_0,int_A)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*3)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*3)

  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,f=f, gamma=gamma, r=r, year=F, maxtime = 365*6)
  row.names(simul2)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2 %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step) ,
               tolerance = 1e-05, label = "chaining from equilibrium is not the same as simulating from equilibrium")
})


test_that("test simulation from equilibrium and chained, with time varying vector control, with delay", {

  # simulate a vector control intervention with 3 years half life, distributed every 4 years
  my_omega=vector_control_exponential_decay(initial_omega = 0.7, half_life = 3, every_x_years = 4, maxtime = 70*365 )

  mydata=data.frame(incidence=c(23,112),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0,0)
  mydata$omega=c(0.7,0.7)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata_delay(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "rho.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=NA, "beta.new"=NA, "omega.new"=my_omega, "sigma.new"=NA, "rho.new"=NA)

  my_intervention_list=list(int_0,int_A)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*3, delay=T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*3, delay=T)

  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,f=f, gamma=gamma, r=r, year=F, maxtime = 365*6, delay=T)
  row.names(simul2)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2 %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step) ,
               tolerance = 1e-05, label = "chaining from equilibrium is not the same as simulating from equilibrium")
})

