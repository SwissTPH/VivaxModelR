test_that("test time varying delta working", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      alpha.old=0.95*rho,
      lambda=c(0.02560546  ,0.017086074  ,0.022004378 ),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3))

  my_delta1=data.frame(t=c(0,365,1000, 4001), value=c(1,1,0,0)*3.562799e-05)
  my_delta2=data.frame(t=c(0,365,1000, 4001), value=c(1,1,0,0)*2.694904e-04)
  my_delta3=data.frame(t=c(0,365,1000, 4001), value=c(1,1,0,0)*1.854486e-03)
  my_delta_t1 <- approxfun(my_delta1)
  my_delta_t2 <- approxfun(my_delta2)
  my_delta_t3 <- approxfun(my_delta3)

  mydata_tv=mydata
  mydata_tv$delta.new=c(my_delta_t1, my_delta_t2, my_delta_t3)
  simul.tv=simulate_from_data(mydata_tv,  f=1/69, gamma=1/383, r=1/60,  maxtime=4000,year=F)

  mydata_1=mydata
  simul.1=simulate_from_data(mydata_1,  f=1/69, gamma=1/383, r=1/60,  maxtime=4000,year=F)

  mydata_05=mydata
  mydata_05$delta=mydata$delta*0
  simul.05=simulate_from_data(mydata_05,  f=1/69, gamma=1/383, r=1/60,  maxtime=4000,year=F)

  # ggplot(simul.tv)+
  #   geom_line(aes(x=time, y=I, color=as.factor(id), linetype="decay"), lwd=1)+
  #   geom_line(data=simul.1,aes(x=time, y=I, color=as.factor(id), linetype="1"), lwd=1)+
  #   geom_line(data=simul.05,aes(x=time, y=I, color=as.factor(id), linetype="0"), lwd=1)
  #

  expect_equal(simul.tv[simul.tv$time<1*365,], simul.1[simul.1$time<1*365,],tolerance = 1e-07, label = "before time varying")
  expect_equal(simul.tv[simul.tv$time>4990,], simul.05[simul.05$time>4990,],tolerance = 1e-05, label = "after time varying")
  expect_lt(simul.tv$I[simul.tv$time==500 & simul.tv$id==1],simul.1$I[simul.1$time==500  & simul.tv$id==1], label = "between time varying, lower than 1")
  expect_lt(simul.tv$I[simul.tv$time==500 & simul.tv$id==2],simul.1$I[simul.1$time==500 & simul.tv$id==2], label = "between time varying, lower than 1")
  expect_lt(simul.tv$I[simul.tv$time==500 & simul.tv$id==3],simul.1$I[simul.1$time==500 & simul.tv$id==3], label = "between time varying, lower than 1")
  expect_gt(simul.tv$I[simul.tv$time==500 & simul.tv$id==1],simul.05$I[simul.05$time==500  & simul.05$id==1], label = "between time varying, higher than 0.7")
  expect_gt(simul.tv$I[simul.tv$time==500 & simul.tv$id==2],simul.05$I[simul.05$time==500 & simul.05$id==2], label = "between time varying, higher than 0.7")
  expect_gt(simul.tv$I[simul.tv$time==500 & simul.tv$id==3],simul.05$I[simul.05$time==500 & simul.05$id==3], label = "between time varying, higher than 0.7")
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

  my_delta1=data.frame(t=c(0,365,1000, 4001), value=c(1,1,0,0)*3.562799e-05)
  my_delta2=data.frame(t=c(0,365,1000, 4001), value=c(1,1,0,0)*2.694904e-04)
  my_delta3=data.frame(t=c(0,365,1000, 4001), value=c(1,1,0,0)*1.854486e-03)
  my_delta_t1 <- approxfun(my_delta1)
  my_delta_t2 <- approxfun(my_delta2)
  my_delta_t3 <- approxfun(my_delta3)

  mydata_tv=mydata
  mydata_tv$delta.new=c(my_delta_t1, my_delta_t2, my_delta_t3)
  simul.tv=simulate_from_data_delay(mydata_tv,  f=1/69, gamma=1/383, r=1/60,  maxtime=4000,year=F)

  mydata_1=mydata
  simul.1=simulate_from_data_delay(mydata_1,  f=1/69, gamma=1/383, r=1/60,  maxtime=4000,year=F)

  mydata_05=mydata
  mydata_05$delta=mydata$delta*0
  simul.05=simulate_from_data_delay(mydata_05,  f=1/69, gamma=1/383, r=1/60,  maxtime=4000,year=F)

  # ggplot(simul.tv)+
  #   geom_line(aes(x=time, y=I, color=as.factor(id), linetype="decay"), lwd=1)+
  #   geom_line(data=simul.1,aes(x=time, y=I, color=as.factor(id), linetype="1"), lwd=1)+
  #   geom_line(data=simul.05,aes(x=time, y=I, color=as.factor(id), linetype="0"), lwd=1)
  #

  expect_equal(simul.tv[simul.tv$time<1*365,], simul.1[simul.1$time<1*365,],tolerance = 1e-07, label = "before time varying")
  expect_equal(simul.tv[simul.tv$time>4990,], simul.05[simul.05$time>4990,],tolerance = 1e-05, label = "after time varying")
  expect_lt(simul.tv$I[simul.tv$time==500 & simul.tv$id==1],simul.1$I[simul.1$time==500  & simul.tv$id==1], label = "between time varying, lower than 1")
  expect_lt(simul.tv$I[simul.tv$time==500 & simul.tv$id==2],simul.1$I[simul.1$time==500 & simul.tv$id==2], label = "between time varying, lower than 1")
  expect_lt(simul.tv$I[simul.tv$time==500 & simul.tv$id==3],simul.1$I[simul.1$time==500 & simul.tv$id==3], label = "between time varying, lower than 1")
  expect_gt(simul.tv$I[simul.tv$time==500 & simul.tv$id==1],simul.05$I[simul.05$time==500  & simul.05$id==1], label = "between time varying, higher than 0.7")
  expect_gt(simul.tv$I[simul.tv$time==500 & simul.tv$id==2],simul.05$I[simul.05$time==500 & simul.05$id==2], label = "between time varying, higher than 0.7")
  expect_gt(simul.tv$I[simul.tv$time==500 & simul.tv$id==3],simul.05$I[simul.05$time==500 & simul.05$id==3], label = "between time varying, higher than 0.7")
})


test_that("test simulation from equilibrium and chained, with time varying vector control", {

  # simulate a vector control intervention with 3 years half life, distributed every 4 years

  mydata=data.frame(incidence=c(23,112),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$prop_import=c(0.01,0.01)
  mydata$omega=c(1,0.7)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata(mydata,f=f, gamma=gamma, r=r, return.all = T )

  full=list()
  full$id=mydata2$id
  full$time=c(0,365,2*365,6*365+1)
  my_delta0=expand.grid( full )
  my_decay=data.frame(time=c(0,365,2*365,6*365+1), value_fact=c(1,1,1.5,0))

  my_delta=dplyr::left_join(my_delta0, mydata2[,c("id","delta")]) %>% dplyr::left_join(my_decay) %>%
    dplyr::mutate(value=value_fact*delta) %>%
    dplyr::select(id,time,value)

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "rho.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "rho.new"=NA, "delta.new"=my_delta)

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

  my_delta=data.frame(time=c(0,365,6*365+1), value=c(3.562875e-06,3.562875e-05,0))


  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA,"rho.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "rho.new"=NA, "sigma.new"=NA, "delta.new"=my_delta)

  my_intervention_list=list(int_0,int_A)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*3, delay = T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*3, delay = T)

  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,f=f, gamma=gamma, r=r, year=F, maxtime = 365*6, delay = T)
  row.names(simul2)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2 %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step) ,
               tolerance = 1e-05, label = "chaining from equilibrium is not the same as simulating from equilibrium")
})



