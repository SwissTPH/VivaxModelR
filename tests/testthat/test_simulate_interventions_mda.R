
test_that("test simulation of future scenarios, with MDA", {

  mydata=data.frame(incidence=c(273,212),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$prop_import=c(0.1,0)
  mydata$omega=c(1,1)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA)
  int_A=list(intervention_name="A","rho.new"=NA, "alpha.new"=0.2, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA)
  int_0M=list(intervention_name="MDA", "rho.new"=NA,"alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0)
  int_AM=list(intervention_name="A+MDA","rho.new"=NA, "alpha.new"=0.2, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0)
  my_intervention_list=list(int_0,int_A, int_0M,int_AM)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*1,mda = F)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*10, mda = T)

  row.names(simul2)=NULL
  expect_equal(simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A", "baseline")) %>% dplyr::select(-intervention),
               simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A+MDA", "MDA")) %>% dplyr::select(-intervention),
               tolerance = 5e-05, label = "MDA effect disappears in the long term")

  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="baseline" & simul2$id==1],
            label = "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="baseline" & simul2$id==2],
            label =  "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="A" & simul2$id==1],
            label =  "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="A" & simul2$id==2],
            label =  "MDA better than no MDA")
})

test_that("test simulation of future scenarios, with MDA and delay", {

  mydata=data.frame(incidence=c(273,212),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0.1,0)
  mydata$omega=c(0.7,1)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata_delay(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA)
  int_A=list(intervention_name="A","rho.new"=NA, "alpha.new"=0.22, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA)
  int_0M=list(intervention_name="MDA","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0)
  int_AM=list(intervention_name="A+MDA","rho.new"=NA, "alpha.new"=0.22, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0)
  my_intervention_list=list(int_0,int_A, int_0M,int_AM)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*1, delay=T,mda = F)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*10, delay=T, mda = T)

  row.names(simul2)=NULL
  expect_equal(simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A", "baseline")) %>% dplyr::select(-intervention),
               simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A+MDA", "MDA")) %>% dplyr::select(-intervention),
               tolerance = 5e-05, label = "MDA effect disappears in the long term")

  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="baseline" & simul2$id==1],
            label = "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="baseline" & simul2$id==2],
            label =  "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="A" & simul2$id==1],
            label =  "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="A" & simul2$id==2],
            label =  "MDA better than no MDA")
})



test_that("test simulation of future scenarios, with MDA and RCD", {

  mydata=data.frame(incidence=c(273,212),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$prop_import=c(0.1,0)
  mydata$omega=c(1,1)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA, "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA)
  int_A=list(intervention_name="A","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5)
  int_0M=list(intervention_name="MDA","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0, "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA)
  int_AM=list(intervention_name="A+MDA","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5)
  my_intervention_list=list(int_0,int_A, int_0M,int_AM)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*1,mda = F, rcd=T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*10, mda = T, rcd=T)

  row.names(simul2)=NULL
  expect_equal(simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A", "baseline")) %>% dplyr::select(-intervention),
               simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A+MDA", "MDA")) %>% dplyr::select(-intervention),
               tolerance = 5e-05, label = "MDA effect disappears in the long term")

  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="baseline" & simul2$id==1],
            label = "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="baseline" & simul2$id==2],
            label =  "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="A" & simul2$id==1],
            label =  "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="A" & simul2$id==2],
            label =  "MDA better than no MDA")
})



test_that("test simulation of future scenarios, with MDA and delay and RCD", {

  mydata=data.frame(incidence=c(273,212),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0.1,0)
  mydata$omega=c(0.7,1)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata_delay(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline", "rho.new"=NA,"alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA, "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA)
  int_A=list(intervention_name="A","rho.new"=NA, "alpha.new"=0.22, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA, "iota.new"=3/7/10000, "nu.new"=3, "eta.new"=1, "tau.new"=5)
  int_0M=list(intervention_name="MDA", "rho.new"=NA,"alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0, "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA)
  int_AM=list(intervention_name="A+MDA", "rho.new"=NA,"alpha.new"=0.22, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=30, "MDArad_cure.new"=0, "iota.new"=3/7/10000, "nu.new"=3, "eta.new"=1, "tau.new"=5)
  my_intervention_list=list(int_0,int_A, int_0M,int_AM)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*1, delay=T,mda = F, rcd=T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*10, delay=T, mda = T, rcd=T)

  simul2_ref=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*10, delay=T, mda = T, rcd=T, referral = T)

  row.names(simul2)=NULL
  row.names(simul2_ref)=NULL
  expect_equal(simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A", "baseline")) %>% dplyr::select(-intervention),
               simul2 %>% dplyr::filter(time==11*365, intervention %in% c("A+MDA", "MDA")) %>% dplyr::select(-intervention),
               tolerance = 1e-04, label = "MDA effect disappears in the long term")

  expect_equal(simul2_ref %>% dplyr::filter(time==11*365, intervention %in% c("A", "baseline")) %>% dplyr::select(-intervention),
               simul2_ref %>% dplyr::filter(time==11*365, intervention %in% c("A+MDA", "MDA")) %>% dplyr::select(-intervention),
               tolerance = 1e-04, label = "MDA effect disappears in the long term")

  expect_equal(simul2 %>% dplyr::filter(time==11*365, intervention %in% c("baseline", "MDA")),
               simul2_ref %>% dplyr::filter(time==11*365, intervention %in% c("baseline", "MDA")),
               tolerance = 1e-09, label = "MDA effect disappears in the long term")


  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="baseline" & simul2$id==1],
            label = "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="baseline" & simul2$id==2],
            label =  "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==1],
            simul2$I[simul2$time==366 & simul2$intervention =="A" & simul2$id==1],
            label =  "MDA better than no MDA")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==2],
            simul2$I[simul2$time==366 & simul2$intervention =="A" & simul2$id==2],
            label =  "MDA better than no MDA")


  expect_lt(simul2_ref$I[simul2_ref$time==366 & simul2_ref$intervention =="MDA" & simul2$id==1],
            simul2_ref$I[simul2_ref$time==366 & simul2_ref$intervention =="baseline" & simul2$id==1],
            label = "MDA better than no MDA")
  expect_lt(simul2_ref$I[simul2_ref$time==366 & simul2_ref$intervention =="MDA" & simul2_ref$id==2],
            simul2_ref$I[simul2_ref$time==366 & simul2_ref$intervention =="baseline" & simul2_ref$id==2],
            label =  "MDA better than no MDA")
  expect_lt(simul2_ref$I[simul2_ref$time==366 & simul2_ref$intervention =="A+MDA" & simul2_ref$id==1],
            simul2_ref$I[simul2_ref$time==366 & simul2_ref$intervention =="A" & simul2_ref$id==1],
            label =  "MDA better than no MDA")
  expect_lt(simul2_ref$I[simul2_ref$time==366 & simul2_ref$intervention =="A+MDA" & simul2_ref$id==2],
            simul2_ref$I[simul2_ref$time==366 & simul2_ref$intervention =="A" & simul2_ref$id==2],
            label =  "MDA better than no MDA")


  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==1],
            simul2_ref$I[simul2_ref$time==366 & simul2_ref$intervention =="A+MDA" & simul2$id==1],
            label = "no referral better than referral")
  expect_lt(simul2$I[simul2$time==366 & simul2$intervention =="A+MDA" & simul2$id==2],
            simul2_ref$I[simul2_ref$time==366 & simul2_ref$intervention =="A+MDA" & simul2$id==2],
            label = "no referral better than referral")
})

