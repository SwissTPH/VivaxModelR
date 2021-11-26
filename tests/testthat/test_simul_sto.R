test_that("test time varying vector control working", {

  mySTOmodel=model_sto_vivax_delay()

  parameters=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.03,"delta"=1.854486e-05,
                  "alpha"=0.4, "beta"=0.7, "sigma"=1/15, "rho"=0.4,"omega"=1,
                  "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hl"=0,"hh"=0, "hhl"=0, "T0"=0, "Tl"=0, "N"=10000)

  # fixed parameters
  simul_sto=simulate_vivax_delay_sto(parameters, STOmodel=mySTOmodel, runs = 10, maxtime = 1465, year=T)
  simul_sto_mean=simul_sto %>% dplyr::group_by(.data$time) %>% dplyr::summarise_all(.funs = "mean") %>% data.frame()
  simul_sto_mean[,c(-1,-2,-14)]=simul_sto_mean[,c(-1,-2,-14)]/parameters$N

  simul=simulate_vivax_delay_ode(parameters, ODEmodel =ode_vivax_delay , maxtime = 1465, year=T)
  simul$run=5.5
  row.names(simul)=1:5
  expect_equal(simul_sto_mean, simul[,names(simul_sto_mean)],tolerance = 3e-02, label = "sto close to det")

  simul_sto_minmax=simul_sto %>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -run)) %>%
    dplyr::group_by(.data$time, .data$name)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
      dplyr::left_join(simul%>% dplyr::select(-p,-run) %>% tidyr::pivot_longer(cols=c(-time))%>% data.frame())  %>%
      dplyr::mutate(is.max=(max/parameters$N>=value), is.min=(min/parameters$N<=value))

  expect_true(all(simul_sto_minmax$is.max ==TRUE), label = "max sto >= det")
  expect_true(all(simul_sto_minmax$is.min ==TRUE), label = "min sto <= det")

  # ggplot(simul_sto )+
  #   geom_line(aes(x=time, y=h/parameters$N, color="sto", group=as.factor(run)))+
  #   geom_line(data=simul,aes(x=time, y=h), color="black")
  #

  # time varying omega
  my_omega=vector_control_exponential_decay(initial_omega=0.5, half_life=2, every_x_years=2, maxtime=5*365)
  omega_t <- approxfun(my_omega)
  parameters_tv=parameters
  parameters_tv$omega=omega_t

  simultv_sto=simulate_vivax_delay_sto(parameters_tv, STOmodel=mySTOmodel, runs = 10, year=T)
  simultv_sto_mean=simultv_sto %>% dplyr::group_by(.data$time) %>% dplyr::summarise_all("mean") %>% data.frame()
  simultv_sto_mean[,c(-1,-2,-14)]=simultv_sto_mean[,c(-1,-2,-14)]/parameters$N

  simultv=simulate_vivax_delay_ode(parameters_tv, ODEmodel =ode_vivax_delay , maxtime = 1465, year=T)
  simultv$run=5.5
  row.names(simultv)=1:5
  expect_equal(simultv_sto_mean, simultv[,names(simultv_sto_mean)],tolerance = 3e-02, label = "sto close to det, omega tv")

  simultv_sto_minmax=simultv_sto %>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -run)) %>%
    dplyr::group_by(.data$time, .data$name)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(simultv%>% dplyr::select(-p,-run) %>% tidyr::pivot_longer(cols=c(-time))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max/parameters$N>=value), is.min=(min/parameters$N<=value))

  expect_true(all(simultv_sto_minmax$is.max ==TRUE), label = "max sto >= det")
  expect_true(all(simultv_sto_minmax$is.min ==TRUE), label = "min sto <= det")

   # ggplot(simultv_sto )+
   #   geom_line(aes(x=time, y=T0/parameters$N, color="sto", group=as.factor(run)))+
   #   geom_line(data=simultv,aes(x=time, y=T0), color="black")

  # time varying delta
  my_delta=data.frame(t=c(0,600,1000, 4001), value=c(1,0,1,0)*8e-03)
  my_delta_t <- approxfun(my_delta)
  parameters_delta=parameters
  parameters_delta$delta=my_delta_t

  simuldelta_sto=simulate_vivax_delay_sto(parameters_delta, STOmodel=mySTOmodel, runs = 10, year=T)
  simuldelta_sto_mean=simuldelta_sto %>% dplyr::group_by(.data$time) %>% dplyr::summarise_all("mean") %>% data.frame()
  simuldelta_sto_mean[,c(-1,-2,-14)]=simuldelta_sto_mean[,c(-1,-2,-14)]/parameters$N

  simuldelta=simulate_vivax_delay_ode(parameters_delta, ODEmodel =ode_vivax_delay , maxtime = 1465, year=T)
  simuldelta$run=5.5
  row.names(simuldelta)=1:5
  expect_equal(simuldelta_sto_mean[,c(-14)], simuldelta[,names(simuldelta_sto_mean)[1:13]],tolerance = 3e-02, label = "sto close to det, delta tv")

  simuldelta_sto_minmax=simuldelta_sto %>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -run)) %>%
    dplyr::group_by(.data$time, .data$name)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(simuldelta%>% dplyr::select(-p,-run) %>% tidyr::pivot_longer(cols=c(-time))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max/parameters$N>=value), is.min=(min/parameters$N<=value))

  expect_true(all(simuldelta_sto_minmax$is.max ==TRUE), label = "max sto >= det")
  expect_true(all(simuldelta_sto_minmax$is.min ==TRUE), label = "min sto <= det")

})


test_that("test simulation of future scenarios, with stochastic model", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(rho=c(0.18,0.13,0.08) ,
                  rho.old=rho,
                  beta.old=c(0.431,0.429,0.422),
                  sigma.old=c(1/15,1/15,1/15),
                  TQ_effect=c(0.59,0.605,0.619),
                  alpha.old=0.95*rho,
                  id=c(1,2,3),
                  N=c(10000, 2000, 3000),
                  prop_import=c(0, 0, 0.1))%>%
    dplyr::mutate(h=incidence_year2day(incidence),
                  alpha=alpha.old, beta=beta.old, sigma=sigma.old)


  mydata2=calculate_r0_rc_fromdata_delay(df=mydata,   f=1/69, gamma=1/383, r=1/60, return.all = T)
  simul.PQ=simulate_from_data_delay(mydata2,
                                    f=1/69, gamma=1/383, r=1/60,
                                    maxtime=2000,year=T, sto = T, runs=10)
  simul.TQ=simulate_from_data_delay(mydata2 %>% dplyr::mutate(beta.new=TQ_effect),
                                    f=1/69, gamma=1/383, r=1/60,
                                    maxtime=2000,year=T, sto = T, runs=10)
  simul.TQ.det=simulate_from_data_delay(mydata2 %>% dplyr::mutate(beta.new=TQ_effect),
                                    f=1/69, gamma=1/383, r=1/60,
                                    maxtime=2000,year=T)

  simul.PQ_mean=simul.PQ %>% dplyr::group_by(.data$time, .data$id) %>% dplyr::summarise_all("mean") %>% data.frame()

  average_later_years=simul.PQ_mean[simul.PQ_mean$time>0,] %>% dplyr::group_by(.data$id) %>%
    dplyr::summarise_all("mean") %>% data.frame()

  expect_equal(mydata[,c("id", "incidence")],
               average_later_years[,c("id", "incidence")],
               tolerance =2e-02, label = "constant simulation almost constant")

  # ggplot(simul.PQ)+
  #   geom_line(aes(x=time, y=Il, color="Il", group=run))+
  #   geom_line(aes(x=time, y=S0, color="S0", group=run))+
  #   facet_wrap(.~id)
  #
  # ggplot(simul.PQ)+
  #   geom_line(aes(x=time, y=incidence, color="incidence", group=run))+
  #   facet_wrap(.~id)

  simul.TQ_mean=simul.TQ %>% dplyr::group_by(.data$id, .data$time) %>% dplyr::summarise_all("mean") %>% data.frame() %>%
    dplyr::left_join(mydata[c("id", "N")]) %>%
    dplyr::mutate(Il=.data$Il/.data$N, I0=.data$I0/.data$N,  Sl=.data$Sl/.data$N, S0=.data$S0/.data$N, Tl=.data$Tl/.data$N, T0=.data$T0/.data$N,
                  I=.data$I/.data$N, hh=.data$hh/.data$N,hhl=.data$hhl/.data$N,h=.data$h/.data$N,hl=.data$hl/.data$N)

  row.names(simul.TQ.det)=NULL
  row.names(simul.TQ_mean)=NULL

  expect_equal(simul.TQ_mean[,names(simul.TQ.det)], simul.TQ.det,tolerance = 2.5e-02, label = "sto close to det")

  simul.TQ_minmax=simul.TQ %>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -run, -id, -incidence)) %>%
    dplyr::left_join(mydata[c("id", "N")])  %>% dplyr::mutate(value=value/N) %>% dplyr::select(-N) %>%
    dplyr::group_by(.data$time, .data$name, .data$id)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(simul.TQ.det%>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -id))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max>=value), is.min=(min<=value)) %>% dplyr::filter(time>0)

  expect_true(all(simul.TQ_minmax$is.max ==TRUE | simul.TQ_minmax$max==0), label = "max sto >= det")
  expect_true(all(simul.TQ_minmax$is.min ==TRUE), label = "min sto <= det")

  # ggplot(simul.TQ %>% dplyr::filter(id==1))+
  #   geom_line(aes(x=time, y=T0/10000, color="S0", group=run))+
  #   geom_line(data=simul.TQ.det%>% dplyr::filter(id==1),aes(x=time, y=S0), color="black")

})



test_that("test simulation of future scenarios", {

  mydata=data.frame(incidence=c(223,152),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0,0.01)
  mydata$omega=c(1,1)
  mydata$N=c(3000,9000)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2= calibrate_vivax_equilibrium(mydata,f=f, gamma=gamma, r=r, delay = T,return.all = T )

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA ,"rho.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=0.2, "beta.new"=0.5, "omega.new"=0.9, "sigma.new"=1/5, "rho.new"=0.2)
  my_intervention_list=list(int_0,int_A)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list, delay=T,
                                      f=f, gamma=gamma, r=r, year=T, maxtime = 365*5, sto=TRUE, runs = 15)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1, delay=T,
                                      f=f, gamma=gamma, r=r, year=T ,maxtime = 365*10, sto=TRUE, runs=15)

  simul2_mean=simul2 %>% dplyr::group_by(.data$time, .data$id, .data$intervention) %>% dplyr::summarise_all("mean") %>% data.frame()%>%
    dplyr::left_join(mydata[c("id", "N")]) %>%
    dplyr::mutate(Il=.data$Il/.data$N, I0=.data$I0/.data$N,  Sl=.data$Sl/.data$N, S0=.data$S0/.data$N, Tl=.data$Tl/.data$N, T0=.data$T0/.data$N,
                  I=.data$I/.data$N, hh=.data$hh/.data$N,hhl=.data$hhl/.data$N,h=.data$h/.data$N,hl=.data$hl/.data$N)


  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,f=f, gamma=gamma, r=r, year=T, maxtime = 365*15, delay=T)
  row.names(simul2)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2_mean %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step, -run, -N, -p, -incidence),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(time, id, intervention, Il, I0, Sl, S0, Tl, T0, hh, hhl, h, hl, I ) ,
               tolerance = 5e-02, label = "chaining from equilibrium is not the same as simulating from equilibrium")

  simul2_minmax=simul2 %>% dplyr::select(-p, -step) %>% tidyr::pivot_longer(cols=c(-time, -run, -id, -incidence, -intervention)) %>%
    dplyr::left_join(mydata[c("id", "N")])  %>% dplyr::mutate(value=value/N) %>% dplyr::select(-N) %>%
    dplyr::group_by(.data$time, .data$name, .data$id, .data$intervention)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(simul_1_2%>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -id, -intervention))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max>=value), is.min=(min<=value)) %>% dplyr::filter(time>0)

  expect_true(all(simul2_minmax$is.max ==TRUE | simul2_minmax$max==0), label = "max sto >= det")
  expect_true(all(simul2_minmax$is.min ==TRUE), label = "min sto <= det")

  # ggplot(simul2)+
  #   geom_line(aes(x=time, y=incidence, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul_1_2,aes(x=time, y=incidence, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id, scale="free")
  #


})



test_that("test compare MDA in delay model with non MDA model (no RCD)", {

  parameters_mda=list("r"=1/60, "gamma"=1/383,
                      "f"=1/69,"lambda"=0.013,"delta"=3.562799e-05,
                      "alpha"=0.18*0.95, "beta"=0.431, "rho"=0.18,"omega"=1, "sigma"=1/15,
                      "I0"=0.037, "S0"=0.38, "Sl"=0.2686, "Il"=0.3, "Tl"=0.014, "T0"=0.0004,
                      "h"=0.001, "hl"=0, "hh"=0, "hhl"=0,
                      "MDAcov"=0.5, "MDAp_length"=100, "MDArad_cure"=0.5, "N"=10000)

  parameters_mda_0c=parameters_mda; parameters_mda_0c$MDAcov=0
  parameters_mda_0p=parameters_mda; parameters_mda_0p$MDAp_length=1
  parameters_mda_0r=parameters_mda; parameters_mda_0r$MDArad_cure=0

  mySTOmodel=model_sto_vivax_delay()
  mySTOmodel_mda=model_sto_vivax_delay_mda()

  mda_cp=simulate_vivax_delay_mda_sto(parameters=parameters_mda , STOmodel =mySTOmodel, STOmodel_mda  =mySTOmodel_mda, maxtime=1000, year=T, runs =10)
  mda_cp_det=simulate_vivax_delay_mda_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay, ODEmodel_mda =ode_vivax_delay_mda, maxtime=1000, year=T)

  mda_0p=simulate_vivax_delay_mda_sto(parameters=parameters_mda_0p , STOmodel =mySTOmodel, STOmodel_mda  =mySTOmodel_mda, maxtime=1000, year=T, runs=15)
  mda_0p_det=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0p , ODEmodel=ode_vivax_delay, ODEmodel_mda =ode_vivax_delay_mda, maxtime=1000, year=T)

  mda_0r=simulate_vivax_delay_mda_sto(parameters=parameters_mda_0r , STOmodel =mySTOmodel, STOmodel_mda  =mySTOmodel_mda, maxtime=1000, year=T, runs=10)
  mda_0r_det=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0r , ODEmodel=ode_vivax_delay, ODEmodel_mda =ode_vivax_delay_mda, maxtime=1000, year=T)

  mda_cp_mean=mda_cp %>% dplyr::group_by(.data$time) %>% dplyr::summarise_all("mean") %>% data.frame() %>% dplyr::select(-run)
  mda_cp_mean[,c(-1,-13)]=mda_cp_mean[,c(-1,-13)]/parameters_mda$N

  mda_0p_mean=mda_0p %>% dplyr::group_by(.data$time) %>% dplyr::summarise_all("mean") %>% data.frame() %>% dplyr::select(-run)
  mda_0p_mean[,c(-1,-13)]=mda_0p_mean[,c(-1,-13)]/parameters_mda$N

  mda_0r_mean=mda_0r %>% dplyr::group_by(.data$time) %>% dplyr::summarise_all("mean") %>% data.frame() %>% dplyr::select(-run)
  mda_0r_mean[,c(-1,-13)]=mda_0r_mean[,c(-1,-13)]/parameters_mda$N

  row.names(mda_cp_mean)=NULL
  row.names(mda_0p_mean)=NULL
  row.names(mda_0r_mean)=NULL
  row.names(mda_cp_det)=NULL
  row.names(mda_0p_det)=NULL
  row.names(mda_0r_det)=NULL

  expect_equal(mda_cp_mean, mda_cp_det[names(mda_cp_mean)], tolerance = 5e-02, label = "sto and det are equal, with MDA")
  expect_equal(mda_0p_mean, mda_0p_det[names(mda_0p_mean)], tolerance = 5e-02, label = "sto and det are equal, with MDA, zero prophylaxis")
  expect_equal(mda_0r_mean, mda_0r_det[names(mda_0r_mean)], tolerance = 5e-02, label = "sto and det are equal, with MDA, zero radical cure")

 # test if min and max of stochastic simlations are lower and upper bound for deterministic run
  mda_cp_mean_minmax=mda_cp %>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -run)) %>%
    dplyr::group_by(.data$time, .data$name)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(mda_cp_det%>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max/parameters_mda$N>=value), is.min=(min/parameters_mda$N<=value)) %>% dplyr::filter(time>0)

  expect_true(all(mda_cp_mean_minmax$is.max ==TRUE), label = "max sto >= det")
  expect_true(all(mda_cp_mean_minmax$is.min ==TRUE), label = "min sto <= det")

  mda_0p_mean_minmax=mda_0p %>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -run)) %>%
    dplyr::group_by(.data$time, .data$name)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(mda_0p_det%>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max/parameters_mda$N>=value), is.min=(min/parameters_mda$N<=value)) %>% dplyr::filter(time>0)

  expect_true(all(mda_0p_mean_minmax$is.max ==TRUE), label = "max sto >= det")
  expect_true(all(mda_0p_mean_minmax$is.min ==TRUE), label = "min sto <= det")

  mda_0r_mean_minmax=mda_0r %>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -run)) %>%
    dplyr::group_by(.data$time, .data$name)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(mda_0r_det%>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max/parameters_mda$N>=value), is.min=(min/parameters_mda$N<=value)) %>% dplyr::filter(time>0)

  expect_true(all(mda_0r_mean_minmax$is.max ==TRUE), label = "max sto >= det")
  expect_true(all(mda_0r_mean_minmax$is.min ==TRUE), label = "min sto <= det")


    # ggplot(mda_cp )+
    #   geom_line(aes(x=time, y=I/parameters_mda$N, color="sto", group=as.factor(run)))+
    #   geom_line(data=mda_cp_det,aes(x=time, y=I), color="black")
    # ggplot(mda_0p)+
    #   geom_line(aes(x=time, y=I/parameters_mda$N, color="sto", group=as.factor(run)))+
    #   geom_line(data=mda_0p_det,aes(x=time, y=I), color="black")
    # ggplot(mda_0r)+
    #   geom_line(aes(x=time, y=I/parameters_mda$N, color="sto", group=as.factor(run)))+
    #   geom_line(data=mda_0r_det,aes(x=time, y=I), color="black")
    #
    # ggplot(mda_cp )+
    #   geom_line(aes(x=time, y=h/parameters_mda$N, color="sto", group=as.factor(run)))+
    #   geom_line(data=mda_cp_det,aes(x=time, y=h), color="black")
    # ggplot(mda_0p)+
    #   geom_line(aes(x=time, y=h/parameters_mda$N, color="sto", group=as.factor(run)))+
    #   geom_line(data=mda_0p_det,aes(x=time, y=h), color="black")
    # ggplot(mda_0r)+
    #   geom_line(aes(x=time, y=h/parameters_mda$N, color="sto", group=as.factor(run)))+
    #   geom_line(data=mda_0r_det,aes(x=time, y=h), color="black")

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
  mydata$N=c(5000,3000)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata_delay(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA)
  int_A=list(intervention_name="A","rho.new"=NA, "alpha.new"=0.22, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA)
  int_0M=list(intervention_name="MDA","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=100, "MDArad_cure.new"=0)
  int_AM=list(intervention_name="A+MDA","rho.new"=NA, "alpha.new"=0.22, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=100, "MDArad_cure.new"=0)
  my_intervention_list=list(int_0,int_A, int_0M,int_AM)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=T, maxtime = 365*1, delay=T,mda = F, sto=T, runs=10)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=T ,maxtime = 365*10, delay=T, mda = T, sto=T, runs=10)

  simul2_mean=simul2 %>% dplyr::group_by(.data$time, .data$id, .data$intervention) %>% dplyr::summarise_all("mean") %>%
    data.frame()%>% dplyr::select(-run)%>%
    dplyr::left_join(mydata[c("id", "N")]) %>%
    dplyr::mutate(Il=.data$Il/.data$N, I0=.data$I0/.data$N,  Sl=.data$Sl/.data$N, S0=.data$S0/.data$N, Tl=.data$Tl/.data$N, T0=.data$T0/.data$N,
                  I=.data$I/.data$N, hh=.data$hh/.data$N,hhl=.data$hhl/.data$N,h=.data$h/.data$N,hl=.data$hl/.data$N)

  row.names(simul2_mean)=NULL

  simul1_det=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=T, maxtime = 365*1, delay=T,mda = F)

  simul2_det=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1_det,
                                      f=f, gamma=gamma, r=r, year=T ,maxtime = 365*10, delay=T, mda = T)
  row.names(simul2_det)=NULL

  expect_equal(simul2_mean  %>% dplyr::arrange(id, intervention, time)%>% dplyr::select(-N, -p , - incidence),
               simul2_det  %>% dplyr::arrange(id, intervention, time)%>% dplyr::select(time, id, intervention, Il, I0, Sl, S0, Tl, T0, hh, hhl, h, hl, I, step ),
               tolerance = 5e-02, label = "MDA sto and non sto")

    # ggplot(simul2 )+
    #   geom_line(aes(x=time, y=I/2000, color=as.factor(intervention), group=as.factor(run)))+
    #   geom_line(data=simul2_det,aes(x=time, y=I, group=as.factor(intervention)), color="black")+
    #   facet_wrap(intervention~id)

    # ggplot(simul2 )+
    #   geom_line(aes(x=time, y=incidence, color=as.factor(intervention), group=as.factor(run)))+
    #   geom_line(data=simul2_det,aes(x=time, y=incidence, group=as.factor(intervention)), color="black")+
    #   facet_wrap(intervention~id)


  # simul2_minmax=simul2 %>% dplyr::select(-p, -step) %>% tidyr::pivot_longer(cols=c(-time, -run, -id, -incidence, -intervention)) %>%
  #   dplyr::left_join(mydata[c("id", "N")])  %>% dplyr::mutate(value=value/N) %>% dplyr::select(-N) %>%
  #   dplyr::group_by(.data$time, .data$name, .data$id, .data$intervention)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
  #   dplyr::left_join(simul2_det%>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -id, -intervention))%>% data.frame())  %>%
  #   dplyr::mutate(is.max=(max>=value), is.min=(min<=value)) %>% dplyr::filter(time>0)
  #
  # expect_true(all(simul2_minmax$is.max ==TRUE | simul2_minmax$max==0), label = "max sto >= det")
  # expect_true(all(simul2_minmax$is.min ==TRUE), label = "min sto <= det")

})


