test_that("test simulation of future scenarios", {

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
  mydata2=calibrate_vivax_equilibrium(mydata,f=f, gamma=gamma, r=r, return.all = T)

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "rho.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=0.3, "beta.new"=0.6, "omega.new"=0.9, "rho.new"=0.3)
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
               tolerance = 1e-07, label = "chaining from equilibrium is not the same as simulating from equilibrium")

  simul_1_2$rho=ifelse(simul_1_2$intervention=="A" & simul_1_2$time>0, 0.3,ifelse(simul_1_2$id==1, 0.18,0.13))
  expect_equal(simul_1_2[simul_1_2$id==1,c("h", "hh")], simul_1_2[simul_1_2$id==1,c("h", "hh")], label = "no importation: h=hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$h[1000], simul_1_2$hl[1000], label = "with importation: h>hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$hh[1000], simul_1_2[simul_1_2$id==2,]$hhl[1000], label = "with importation: hh>hhl")
  expect_equal(simul_1_2$h/simul_1_2$rho, simul_1_2$hh, label = "h=rho*hh")

  expect_equal(simul1 %>% dplyr::filter(intervention=="baseline", time==0) %>% dplyr::select(-time),
               simul1 %>% dplyr::filter(intervention=="baseline", time==1000) %>% dplyr::select(-time),
               tolerance = 1e-07, label = "baseline stays constant")
})


test_that("test simulation of future scenarios, with delays", {

  mydata=data.frame(incidence=c(23,112),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0,0.01)
  mydata$omega=c(1,1)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calibrate_vivax_equilibrium(mydata,f=f, gamma=gamma, r=r, return.all = T , delay = T)

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA ,"rho.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=0.3, "beta.new"=0.6, "omega.new"=0.9, "sigma.new"=1/5, "rho.new"=0.3)
  my_intervention_list=list(int_0,int_A)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*3, delay = T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*3, delay=T)

  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                         f=f, gamma=gamma, r=r, year=F, maxtime = 365*6, delay = T)
  row.names(simul2)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2 %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step) ,
               tolerance = 1e-07, label = "chaining from equilibrium is not the same as simulating from equilibrium")

  simul_1_2$rho=ifelse(simul_1_2$intervention=="A" & simul_1_2$time>0, 0.3,ifelse(simul_1_2$id==1, 0.18,0.13))
  expect_equal(simul_1_2[simul_1_2$id==1,c("h", "hh")], simul_1_2[simul_1_2$id==1,c("h", "hh")], label = "no importation: h=hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$h[1000], simul_1_2$hl[1000], label = "with importation: h>hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$hh[1000], simul_1_2[simul_1_2$id==2,]$hhl[1000], label = "with importation: hh>hhl")
  expect_equal(simul_1_2$h/simul_1_2$rho, simul_1_2$hh, label = "h=rho*hh")

  expect_equal(simul1 %>% dplyr::filter(intervention=="baseline", time==0) %>% dplyr::select(-time),
               simul1 %>% dplyr::filter(intervention=="baseline", time==1000) %>% dplyr::select(-time),
               tolerance = 1e-07, label = "baseline stays constant")


})


test_that("test simulation of future scenarios, with RCD", {

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
  mydata2=calibrate_vivax_equilibrium(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA, "rho.new"=NA, "kappa.new"=0.18)
  int_A=list(intervention_name="A", "alpha.new"=0.3, "beta.new"=0.6, "omega.new"=0.9, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  my_intervention_list=list(int_0,int_A)
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

  simul_1_2$rho=ifelse(simul_1_2$intervention=="A" & simul_1_2$time>0, 0.3,ifelse(simul_1_2$id==1, 0.18,0.13))
  expect_equal(simul_1_2[simul_1_2$id==1,c("h", "hh")], simul_1_2[simul_1_2$id==1,c("h", "hh")], label = "no importation: h=hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$h[1000], simul_1_2$hl[1000], label = "with importation: h>hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$hh[1000], simul_1_2[simul_1_2$id==2,]$hhl[1000], label = "with importation: hh>hhl")
  expect_gt(simul_1_2$h[1000]/simul_1_2$rho[1000], simul_1_2$hh[1000], label = "h/rho>hh")

  expect_equal(simul1 %>% dplyr::filter(intervention=="baseline", time==0) %>% dplyr::select(-time),
               simul1 %>% dplyr::filter(intervention=="baseline", time==1000) %>% dplyr::select(-time),
               tolerance = 1e-07, label = "baseline stays constant")

})



test_that("test simulation of future scenarios, with RCD and delay", {

  mydata=data.frame(incidence=c(23,112),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0,0.01)
  mydata$omega=c(1,1)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calibrate_vivax_equilibrium(mydata,f=f, gamma=gamma, r=r, return.all = T , delay = T)

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA, "rho.new"=NA, "kappa.new"=0.18)
  int_A=list(intervention_name="A", "alpha.new"=0.3, "beta.new"=0.6, "omega.new"=0.9, "sigma.new"=1/5, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  my_intervention_list=list(int_0,int_A)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*3, rcd=T, delay=T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*3, rcd=T, delay=T)

  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                         f=f, gamma=gamma, r=r, year=F, maxtime = 365*6, rcd=T, delay=T)
  row.names(simul2)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2 %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step) ,
               tolerance = 1e-07, label = "chaining from equilibrium is not the same as simulating from equilibrium")

  simul_1_2$rho=ifelse(simul_1_2$intervention=="A" & simul_1_2$time>0, 0.3,ifelse(simul_1_2$id==1, 0.18,0.13))
  expect_equal(simul_1_2[simul_1_2$id==1,c("h", "hh")], simul_1_2[simul_1_2$id==1,c("h", "hh")], label = "no importation: h=hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$h[1000], simul_1_2$hl[1000], label = "with importation: h>hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$hh[1000], simul_1_2[simul_1_2$id==2,]$hhl[1000], label = "with importation: hh>hhl")
  expect_gt(simul_1_2$h[1000]/simul_1_2$rho[1000], simul_1_2$hh[1000], label = "h/rho>hh")

  expect_equal(simul1 %>% dplyr::filter(intervention=="baseline", time==0) %>% dplyr::select(-time),
               simul1 %>% dplyr::filter(intervention=="baseline", time==1000) %>% dplyr::select(-time),
               tolerance = 1e-07, label = "baseline stays constant")
})


test_that("test simulation of future scenarios, with RCD and delay and referral", {

  mydata=data.frame(incidence=c(23,220),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0,0.01)
  mydata$omega=c(1,1)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calibrate_vivax_equilibrium(mydata,f=f, gamma=gamma, r=r, return.all = T, delay = T )

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA, "rho.new"=NA, "kappa.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=NA, "beta.new"=0.5, "omega.new"=0.9, "sigma.new"=1/5, "iota.new"=5/7/10000, "nu.new"=2, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  my_intervention_list=list(int_0,int_A)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*3, rcd=T, delay=T, referral = T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*3, rcd=T, delay=T, referral = T)

  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                         f=f, gamma=gamma, r=r, year=F, maxtime = 365*6, rcd=T, delay=T, referral = T)

  simul_1_2_nr=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                         f=f, gamma=gamma, r=r, year=F, maxtime = 365*6, rcd=T, delay=T, referral = F)

  row.names(simul2)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2 %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step) ,
               tolerance = 1e-07, label = "chaining from equilibrium is not the same as simulating from equilibrium")
  expect_lt(simul_1_2_nr$I[365*6] ,simul_1_2$I[365*6], label = "no referral better than referral")

  simul_1_2$rho=ifelse(simul_1_2$intervention=="A" & simul_1_2$time>0, 0.3,ifelse(simul_1_2$id==1, 0.18,0.13))
  expect_equal(simul_1_2[simul_1_2$id==1,c("h", "hh")], simul_1_2[simul_1_2$id==1,c("h", "hh")], label = "no importation: h=hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$h[1000], simul_1_2$hl[1000], label = "with importation: h>hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$hh[1000], simul_1_2[simul_1_2$id==2,]$hhl[1000], label = "with importation: hh>hhl")
  expect_gt(simul_1_2$h[1000]/simul_1_2$rho[1000], simul_1_2$hh[1000], label = "h/rho>hh")

  expect_equal(simul1 %>% dplyr::filter(intervention=="baseline", time==0) %>% dplyr::select(-time),
               simul1 %>% dplyr::filter(intervention=="baseline", time==1000) %>% dplyr::select(-time),
               tolerance = 1e-07, label = "baseline stays constant")
})


test_that("test simulation of future scenarios, with delays", {

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
  mydata1=calibrate_vivax_equilibrium(mydata,f=f, gamma=gamma, r=r, return.all = T)
  mydata2=calculate_r0_rc_fromdata(mydata,f=f, gamma=gamma, r=r, return.all = T )
  expect_equal(mydata1, mydata2, tolerance = 1e-09, label = "test high level function for calibration")
  expect_error(calibrate_vivax_equilibrium(mydata,f=f, gamma=gamma, r=r, return.all = T, rcd=T))

  mydata_r=mydata
  mydata_r$iota=c(5/7/365, Inf)
  mydata_r$nu=c(5,8)
  mydata_r$tau=c(5,5)
  mydata_r$eta=c(0.9,1)
  mydata_r$kappa=c(0.18,0.13)
  mydata1_r=calibrate_vivax_equilibrium(mydata_r,f=f, gamma=gamma, r=r, return.all = T, rcd=T)
  expect_lt(mydata1$R0[1],mydata1_r$R0[1], label = "R0 higher when rcd")
  expect_lt(mydata1$R0[2],mydata1_r$R0[2], label = "R0 higher when rcd")

  mydata_d=mydata
  mydata_d$sigma=c(1/15, 1/15)
  mydata1d=calibrate_vivax_equilibrium(mydata_d,f=f, gamma=gamma, r=r, return.all = T, delay = TRUE )
  mydata2d=calculate_r0_rc_fromdata_delay(mydata_d,f=f, gamma=gamma, r=r, return.all = T )
  expect_equal(mydata1d, mydata2d, tolerance = 1e-09, label = "test high level function for calibration")
  expect_lt(mydata1d$R0[1],mydata1$R0[1], label = "R0 higher when no delays")
  expect_lt(mydata1d$R0[2],mydata1$R0[2], label = "R0 higher when no delays")

  mydata_dr=mydata_r
  mydata_dr$sigma=c(1/15, 1/15)
  expect_error(calibrate_vivax_equilibrium(mydata_d,f=f, gamma=gamma, r=r, return.all = T, rcd=T, delay = TRUE))
  mydata1dr=calibrate_vivax_equilibrium(mydata_dr,f=f, gamma=gamma, r=r, return.all = T, delay = TRUE, rcd=T )
  expect_lt(mydata1d$R0[1],mydata1dr$R0[1], label = "R0 higher when rcd")
  expect_lt(mydata2d$R0[2],mydata1dr$R0[2], label = "R0 higher when rcd")
  expect_lt(mydata1d$R0[1],mydata1_r$R0[1], label = "R0 higher when no delays")
  expect_lt(mydata2d$R0[2],mydata1_r$R0[2], label = "R0 higher when no delays")

})




test_that("test simulation of future scenarios, with RCD and delay (RCD at baseline)", {

  mydata=data.frame(incidence=c(23,112),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0,0.01)
  mydata$iota=c(5/7/10000, 5/7/10000)
  mydata$nu=c(5, 5)
  mydata$tau=c(2, 2)
  mydata$eta=c(1, 1)
  mydata$kappa=c(0.3, 0.3)
  mydata$omega=c(1,1)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calibrate_vivax_equilibrium(mydata,f=f, gamma=gamma, r=r, return.all = T, rcd=T , delay=TRUE)

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA, "rho.new"=NA, "kappa.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=0.3, "beta.new"=0.6, "omega.new"=0.9, "sigma.new"=1/5, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  my_intervention_list=list(int_0,int_A)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*3, rcd=T, delay=T, rcd_at_baseline = T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*3, rcd=T, delay=T, rcd_at_baseline = T)

  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                         f=f, gamma=gamma, r=r, year=F, maxtime = 365*6, rcd=T, delay=T, rcd_at_baseline = T)
  row.names(simul2)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2 %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step) ,
               tolerance = 1e-07, label = "chaining from equilibrium is not the same as simulating from equilibrium")

  simul_1_2$rho=ifelse(simul_1_2$intervention=="A" & simul_1_2$time>0, 0.3,ifelse(simul_1_2$id==1, 0.18,0.13))
  expect_equal(simul_1_2[simul_1_2$id==1,c("h", "hh")], simul_1_2[simul_1_2$id==1,c("h", "hh")], label = "no importation: h=hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$h[1000], simul_1_2$hl[1000], label = "with importation: h>hl")
  expect_gt(simul_1_2[simul_1_2$id==2,]$hh[1000], simul_1_2[simul_1_2$id==2,]$hhl[1000], label = "with importation: hh>hhl")
  expect_gt(simul_1_2$h[1000]/simul_1_2$rho[1000], simul_1_2$hh[1000], label = "h/rho>hh")

  expect_equal(simul1 %>% dplyr::filter(intervention=="baseline", time==0) %>% dplyr::select(-time),
               simul1 %>% dplyr::filter(intervention=="baseline", time==1000) %>% dplyr::select(-time),
               tolerance = 1e-07, label = "baseline stays constant")

})
