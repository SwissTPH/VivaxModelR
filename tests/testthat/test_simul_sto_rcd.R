test_that("test compare RCD in delay model (non referral) with non RCD model", {

  parameters=list("r"=1/60, "gamma"=1/223,
                      "f"=1/72,"lambda"=0.0155531,"delta"=0,
                      "alpha"=0, "beta"=1, "rho"=0.5,"omega"=1,
                      "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0, "hl"=0,"hh"=0, "hhl"=0,
                      "sigma"=1/15, "T0"=0, "Tl"=0,"kappa"=1,
                      "tau"=5, "nu"=5, "iota"=50/7/10000, "eta"=1, "N"=10000)

  mySTOmodel=model_sto_vivax_delay_rcd_no_referral()

  simul_sto=simulate_vivax_delay_sto(parameters=parameters, STOmodel=mySTOmodel, runs = 10, maxtime = 1465, year=T)
  simul_sto_mean=simul_sto %>% dplyr::group_by(.data$time) %>% dplyr::summarise_all("mean") %>% data.frame()
  simul_sto_mean[,c(-1,-2,-14)]=simul_sto_mean[,c(-1,-2,-14)]/parameters$N

  simul=simulate_vivax_delay_ode(parameters, ODEmodel =ode_vivax_delay_rcd_no_referral , maxtime = 1465, year=T)
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
  #   geom_line(aes(x=time, y=I/parameters$N, color="sto", group=as.factor(run)))+
  #   geom_line(data=simul,aes(x=time, y=I), color="black")

  # ggplot(simul_sto )+
  #   geom_line(aes(x=time, y=h/parameters$N, color="sto", group=as.factor(run)))+
  #   geom_line(data=simul,aes(x=time, y=h), color="black")

  parameters=list("r"=1/60, "gamma"=1/223,
                       "f"=1/72,"lambda"=0.0155531,"delta"=0.01,
                       "alpha"=0.2, "beta"=0.7, "rho"=1,"omega"=0.9,
                       "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0,"hl"=0,"hh"=0, "hhl"=0,
                       "sigma"=1/15, "T0"=0, "Tl"=0,"kappa"=1,
                       "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "N"=10000)

  simul_sto=simulate_vivax_delay_sto(parameters=parameters, STOmodel=mySTOmodel, runs = 10, maxtime = 1465, year=T)
  simul_sto_mean=simul_sto %>% dplyr::group_by(.data$time) %>% dplyr::summarise_all("mean") %>% data.frame()
  simul_sto_mean[,c(-1,-2,-14)]=simul_sto_mean[,c(-1,-2,-14)]/parameters$N

  simul=simulate_vivax_delay_ode(parameters, ODEmodel =ode_vivax_delay_rcd_no_referral , maxtime = 1465, year=T)
  simul$run=5.5
  row.names(simul)=1:5
  expect_equal(simul_sto_mean, simul[,names(simul_sto_mean)],tolerance = 2e-02, label = "sto close to det")

  simul_sto_minmax=simul_sto %>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -run)) %>%
    dplyr::group_by(.data$time, .data$name)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(simul%>% dplyr::select(-p,-run) %>% tidyr::pivot_longer(cols=c(-time))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max/parameters$N>=value), is.min=(min/parameters$N<=value))

  expect_true(all(simul_sto_minmax$is.max ==TRUE), label = "max sto >= det")
  expect_true(all(simul_sto_minmax$is.min ==TRUE), label = "min sto <= det")

})

test_that("test compare RCD in delay model (referral) with non RCD model", {

  parameters=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.0155531,"delta"=0,
                  "alpha"=0, "beta"=1, "rho"=1,"omega"=1,
                  "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0,"hl"=0,"hh"=0, "hhl"=0,
                  "sigma"=1/15, "T0"=0, "Tl"=0,"kappa"=1,
                  "tau"=5, "nu"=5, "iota"=50/7/10000, "eta"=1, "N"=10000)

  mySTOmodel=model_sto_vivax_delay_rcd_referral()

  simul_sto=simulate_vivax_delay_sto(parameters=parameters, STOmodel=mySTOmodel, runs = 15, maxtime = 1465, year=T)
  simul_sto_mean=simul_sto %>% dplyr::group_by(.data$time) %>% dplyr::summarise_all("mean") %>% data.frame()
  simul_sto_mean[,c(-1,-2,-14)]=simul_sto_mean[,c(-1,-2,-14)]/parameters$N

  simul=simulate_vivax_delay_ode(parameters, ODEmodel =ode_vivax_delay_rcd_referral , maxtime = 1465, year=T)
  simul$run=8
  row.names(simul)=1:5
  expect_equal(simul_sto_mean, simul[,names(simul_sto_mean)],tolerance = 3e-02, label = "sto close to det")

  simul_sto_minmax=simul_sto %>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -run)) %>%
    dplyr::group_by(.data$time, .data$name)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(simul%>% dplyr::select(-p,-run) %>% tidyr::pivot_longer(cols=c(-time))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max/parameters$N>=value), is.min=(min/parameters$N<=value))

  expect_true(all(simul_sto_minmax$is.max ==TRUE), label = "max sto >= det")
  expect_true(all(simul_sto_minmax$is.min ==TRUE), label = "min sto <= det")

  # ggplot(simul_sto )+
  #   geom_line(aes(x=time, y=I/parameters$N, color="sto", group=as.factor(run)))+
  #   geom_line(data=simul,aes(x=time, y=I), color="black")

  # ggplot(simul_sto )+
  #   geom_line(aes(x=time, y=h/parameters$N, color="sto", group=as.factor(run)))+
  #   geom_line(data=simul,aes(x=time, y=h), color="black")

  parameters=list("r"=1/60, "gamma"=1/223,
                  "f"=1/72,"lambda"=0.0155531,"delta"=0.01,
                  "alpha"=0.2, "beta"=0.7, "rho"=1,"omega"=0.9,
                  "I0"=0, "S0"=0.9, "Sl"=0, "Il"=0.1, "h"=0, "hr"=0,"hl"=0,"hh"=0, "hhl"=0,
                  "sigma"=1/15, "T0"=0, "Tl"=0,"kappa"=1,
                  "tau"=5, "nu"=5, "iota"=5/7/10000, "eta"=1, "N"=10000)

  simul_sto=simulate_vivax_delay_sto(parameters=parameters, STOmodel=mySTOmodel, runs = 15, maxtime = 1465, year=T)
  simul_sto_mean=simul_sto %>% dplyr::group_by(.data$time) %>% dplyr::summarise_all("mean") %>% data.frame()
  simul_sto_mean[,c(-1,-2,-14)]=simul_sto_mean[,c(-1,-2,-14)]/parameters$N

  simul=simulate_vivax_delay_ode(parameters, ODEmodel =ode_vivax_delay_rcd_referral , maxtime = 1465, year=T)
  simul$run=8
  row.names(simul)=1:5
  expect_equal(simul_sto_mean, simul[,names(simul_sto_mean)],tolerance = 2e-02, label = "sto close to det")

  simul_sto_minmax=simul_sto %>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -run)) %>%
    dplyr::group_by(.data$time, .data$name)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(simul%>% dplyr::select(-p,-run) %>% tidyr::pivot_longer(cols=c(-time))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max/parameters$N>=value), is.min=(min/parameters$N<=value))

  expect_true(all(simul_sto_minmax$is.max ==TRUE), label = "max sto >= det")
  expect_true(all(simul_sto_minmax$is.min ==TRUE), label = "min sto <= det")

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
  mydata$N=c(10000,20000)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata_delay(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA, "rho.new"=NA, "kappa.new"=0.18)
  int_A=list(intervention_name="A","alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  my_intervention_list=list(int_0,int_A)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*3, rcd=T, delay=T, sto=T, runs=10)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*3, rcd=T, delay=T,  sto=T, runs=10)

  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                         f=f, gamma=gamma, r=r, year=F, maxtime = 365*6, rcd=T, delay=T)

  simul2_mean=simul2 %>% dplyr::group_by(.data$time, .data$id, .data$intervention) %>% dplyr::summarise_all("mean") %>% data.frame()%>%
    dplyr::left_join(mydata[c("id", "N")]) %>%
    dplyr::mutate(Il=.data$Il/.data$N, I0=.data$I0/.data$N,  Sl=.data$Sl/.data$N, S0=.data$S0/.data$N, Tl=.data$Tl/.data$N, T0=.data$T0/.data$N,
                  I=.data$I/.data$N, hh=.data$hh/.data$N,hhl=.data$hhl/.data$N,h=.data$h/.data$N,hl=.data$hl/.data$N)

  row.names(simul2_mean)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2_mean %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step, -run, -N, -p, -incidence),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(time, id, intervention, Il, I0, Sl, S0, Tl, T0, hh, hhl, h, hl, I ) ,
               tolerance = 4e-02, label = "sto and det are similar")

  expect_equal(simul2_mean %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(S0),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(S0) ,
               tolerance = 4e-03, label = "sto and det are similar")

    # ggplot(simul2 )+
    #   geom_line(aes(x=time, y=I/20000, color=as.factor(intervention), group=as.factor(run)))+  # one area has N=10000 and the other has 20000 that's why the plot looks odd for one facet
    #   geom_line(data=simul_1_2,aes(x=time, y=I, group=as.factor(intervention)), color="black")+
    #   facet_wrap(intervention~id)
    # ggplot(simul2 )+
    #   geom_line(aes(x=time, y=incidence, color=as.factor(intervention), group=as.factor(run)))+
    #   geom_line(data=simul_1_2,aes(x=time, y=incidence, group=as.factor(intervention)), color="black")+
    #   facet_wrap(intervention~id)

    # ggplot(simul2_mean )+
    #   geom_line(aes(x=time, y=I, color=as.factor(intervention)))+
    #   geom_line(data=simul_1_2,aes(x=time, y=I, group=as.factor(intervention)), color="black")+
    #   facet_wrap(intervention~id)

  # year=T
  simul1y=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                       f=f, gamma=gamma, r=r, year=T, maxtime = 365*3, rcd=T, delay=T, sto=T, runs=10)
  simul_1y_det=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                            f=f, gamma=gamma, r=r, year=T, maxtime = 365*3, rcd=T, delay=T)

  simul1y_mean=simul1y %>% dplyr::group_by(.data$time, .data$id, .data$intervention) %>% dplyr::summarise_all("mean") %>% data.frame()%>%
    dplyr::left_join(mydata[c("id", "N")]) %>%
    dplyr::mutate(Il=.data$Il/.data$N, I0=.data$I0/.data$N,  Sl=.data$Sl/.data$N, S0=.data$S0/.data$N, Tl=.data$Tl/.data$N, T0=.data$T0/.data$N,
                  I=.data$I/.data$N, hh=.data$hh/.data$N,hhl=.data$hhl/.data$N,h=.data$h/.data$N,hl=.data$hl/.data$N)

  row.names(simul1y_mean)=NULL
  row.names(simul_1y_det)=NULL

  expect_equal(simul1y_mean %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step, -run, -N, -p, -incidence),
               simul_1y_det%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(time, id, intervention, Il, I0, Sl, S0, Tl, T0, hh, hhl, h, hl, I ) ,
               tolerance = 3e-02, label = "sto and det are similar")

  expect_equal(simul1y_mean %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(S0),
               simul_1y_det%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(S0) ,
               tolerance = 3e-03, label = "sto and det are similar")


  simul1_minmax=simul1y %>% dplyr::select(-p, -step) %>% tidyr::pivot_longer(cols=c(-time, -run, -id, -incidence, -intervention)) %>%
    dplyr::left_join(mydata[c("id", "N")])  %>% dplyr::mutate(value=value/N) %>% dplyr::select(-N) %>%
    dplyr::group_by(.data$time, .data$name, .data$id, .data$intervention)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(simul_1y_det%>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -id, -intervention))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max>=value), is.min=(min<=value)) %>% dplyr::filter(time>0)

  expect_true(all(simul1_minmax$is.max ==TRUE | simul1_minmax$max==0), label = "max sto >= det")
  expect_true(all(simul1_minmax$is.min ==TRUE), label = "min sto <= det")
  # ggplot(simul1y )+
  #   geom_line(aes(x=time, y=I/20000, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul_1y_det,aes(x=time, y=I, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)
  # ggplot(simul1y )+
  #   geom_line(aes(x=time, y=incidence, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul_1y_det,aes(x=time, y=incidence, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)
  #
  # ggplot(simul1y_mean )+
  #   geom_line(aes(x=time, y=I, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul_1y_det,aes(x=time, y=I, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)
})


test_that("test simulation of future scenarios, with RCD and delay, referral", {

  mydata=data.frame(incidence=c(23,112),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0,0.01)
  mydata$omega=c(1,1)
  mydata$N=c(10000,20000)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata_delay(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA, "rho.new"=NA, "kappa.new"=0.18)
  int_A=list(intervention_name="A", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  my_intervention_list=list(int_0,int_A)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*3, rcd=T, delay=T, sto=T, runs=10, referral = T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*3, rcd=T, delay=T,  sto=T, runs=10, referral = T)

  simul_1_2=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                         f=f, gamma=gamma, r=r, year=F, maxtime = 365*6, rcd=T, delay=T, referral = T)

  simul2_mean=simul2 %>% dplyr::group_by(.data$time, .data$id, .data$intervention) %>% dplyr::summarise_all("mean") %>% data.frame()%>%
    dplyr::left_join(mydata[c("id", "N")]) %>%
    dplyr::mutate(Il=.data$Il/.data$N, I0=.data$I0/.data$N,  Sl=.data$Sl/.data$N, S0=.data$S0/.data$N, Tl=.data$Tl/.data$N, T0=.data$T0/.data$N,
                  I=.data$I/.data$N, hh=.data$hh/.data$N,hhl=.data$hhl/.data$N,h=.data$h/.data$N,hl=.data$hl/.data$N)


  row.names(simul2_mean)=NULL
  row.names(simul_1_2)=NULL

  expect_equal(simul2_mean %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step, -run, -N, -p, -incidence),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(time, id, intervention, Il, I0, Sl, S0, Tl, T0, hh, hhl, h, hl, I ) ,
               tolerance = 4e-02, label = "sto and det are similar")

  expect_equal(simul2_mean %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(S0),
               simul_1_2%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(S0) ,
               tolerance = 5e-03, label = "sto and det are similar")

  # ggplot(simul2_mean )+
  #   geom_line(aes(x=time, y=I, color=as.factor(intervention)))+
  #   geom_line(data=simul_1_2,aes(x=time, y=I, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)

  # year=T
  simul1y=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                       f=f, gamma=gamma, r=r, year=T, maxtime = 365*3, rcd=T, delay=T, sto=T, runs=15, referral = T)
  simul_1y_det=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                            f=f, gamma=gamma, r=r, year=T, maxtime = 365*3, rcd=T, delay=T, referral = T)

  simul1y_mean=simul1y %>% dplyr::group_by(.data$time, .data$id, .data$intervention) %>% dplyr::summarise_all("mean") %>% data.frame()%>%
    dplyr::left_join(mydata[c("id", "N")]) %>%
    dplyr::mutate(Il=.data$Il/.data$N, I0=.data$I0/.data$N,  Sl=.data$Sl/.data$N, S0=.data$S0/.data$N, Tl=.data$Tl/.data$N, T0=.data$T0/.data$N,
                  I=.data$I/.data$N, hh=.data$hh/.data$N,hhl=.data$hhl/.data$N,h=.data$h/.data$N,hl=.data$hl/.data$N)

  row.names(simul1y_mean)=NULL
  row.names(simul_1y_det)=NULL

  expect_equal(simul1y_mean %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(-step, -run, -N, -p, -incidence),
               simul_1y_det%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(time, id, intervention, Il, I0, Sl, S0, Tl, T0, hh, hhl, h, hl, I ) ,
               tolerance = 2e-02, label = "sto and det are similar")

  expect_equal(simul1y_mean %>% dplyr::arrange(id, intervention, time) %>% dplyr::select(S0),
               simul_1y_det%>% dplyr::arrange(id, intervention, time) %>% dplyr::select(S0) ,
               tolerance = 2e-03, label = "sto and det are similar")

  simul1_minmax=simul1y %>% dplyr::select(-p, -step) %>% tidyr::pivot_longer(cols=c(-time, -run, -id, -incidence, -intervention)) %>%
    dplyr::left_join(mydata[c("id", "N")])  %>% dplyr::mutate(value=value/N) %>% dplyr::select(-N) %>%
    dplyr::group_by(.data$time, .data$name, .data$id, .data$intervention)%>% dplyr::summarise(max=max(value, na.rm = T), min=min(value, na.rm = T) )%>% data.frame()%>%
    dplyr::left_join(simul_1y_det%>% dplyr::select(-p) %>% tidyr::pivot_longer(cols=c(-time, -id, -intervention))%>% data.frame())  %>%
    dplyr::mutate(is.max=(max>=value), is.min=(min<=value)) %>% dplyr::filter(time>0)

  expect_true(all(simul1_minmax$is.max ==TRUE | simul1_minmax$max==0), label = "max sto >= det")
  expect_true(all(simul1_minmax$is.min ==TRUE), label = "min sto <= det")

  # ggplot(simul1y )+
  #   geom_line(aes(x=time, y=I/20000, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul_1y_det,aes(x=time, y=I, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)
  # ggplot(simul1y )+
  #   geom_line(aes(x=time, y=incidence, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul_1y_det,aes(x=time, y=incidence, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)

  # ggplot(simul1y_mean )+
  #   geom_line(aes(x=time, y=I, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul_1y_det,aes(x=time, y=I, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)
})




test_that("test compare MDA in delay model with non MDA model (with RCD)", {

  parameters_mda=list("r"=1/60, "gamma"=1/383,
                      "f"=1/69,"lambda"=0.013,"delta"=3.562799e-05,
                      "alpha"=0.18*0.95, "beta"=0.431, "rho"=0.18,"omega"=1, "sigma"=1/15,
                      "I0"=0.037, "S0"=0.38, "Sl"=0.2686, "Il"=0.3, "Tl"=0.014, "T0"=0.0004,
                      "h"=0.001, "hl"=0, "hh"=0, "hhl"=0,
                      "MDAcov"=0.5, "MDAp_length"=100, "MDArad_cure"=0.5, "N"=10000,"tau"=5, "nu"=5, "iota"=50/7/10000, "eta"=1,"kappa"=0.18)

  parameters_mda_0c=parameters_mda; parameters_mda_0c$MDAcov=0
  parameters_mda_0p=parameters_mda; parameters_mda_0p$MDAp_length=1
  parameters_mda_0r=parameters_mda; parameters_mda_0r$MDArad_cure=0

  mySTOmodel=model_sto_vivax_delay_rcd_referral()
  mySTOmodel_mda=model_sto_vivax_delay_rcd_referral_mda()

  mda_cp=simulate_vivax_delay_mda_sto(parameters=parameters_mda , STOmodel =mySTOmodel, STOmodel_mda  =mySTOmodel_mda, maxtime=1000, year=T, runs =10)
  mda_cp_det=simulate_vivax_delay_mda_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay_rcd_referral, ODEmodel_mda =ode_vivax_delay_rcd_referral_mda, maxtime=1000, year=T)

  mda_0p=simulate_vivax_delay_mda_sto(parameters=parameters_mda_0p , STOmodel =mySTOmodel, STOmodel_mda  =mySTOmodel_mda, maxtime=1000, year=T, runs=10)
  mda_0p_det=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0p , ODEmodel=ode_vivax_delay_rcd_referral, ODEmodel_mda =ode_vivax_delay_rcd_referral_mda, maxtime=1000, year=T)

  mda_0r=simulate_vivax_delay_mda_sto(parameters=parameters_mda_0r , STOmodel =mySTOmodel, STOmodel_mda  =mySTOmodel_mda, maxtime=1000, year=T, runs=10)
  mda_0r_det=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0r , ODEmodel=ode_vivax_delay_rcd_referral, ODEmodel_mda =ode_vivax_delay_rcd_referral_mda, maxtime=1000, year=T)

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

  expect_equal(mda_cp_mean, mda_cp_det[names(mda_cp_mean)], tolerance = 2e-02, label = "sto and det are equal, with MDA")
  expect_equal(mda_0p_mean, mda_0p_det[names(mda_0p_mean)], tolerance = 3e-02, label = "sto and det are equal, with MDA, zero prophylaxis")
  expect_equal(mda_0r_mean, mda_0r_det[names(mda_0r_mean)], tolerance = 4e-02, label = "sto and det are equal, with MDA, zero radical cure")

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



test_that("test compare MDA in delay model with non MDA model (with RCD non referral)", {

  parameters_mda=list("r"=1/60, "gamma"=1/383,
                      "f"=1/69,"lambda"=0.013,"delta"=3.562799e-05,
                      "alpha"=0.18*0.95, "beta"=0.431, "rho"=0.18,"omega"=1, "sigma"=1/15,
                      "I0"=0.037, "S0"=0.38, "Sl"=0.2686, "Il"=0.3, "Tl"=0.014, "T0"=0.0004,
                      "h"=0.001, "hl"=0, "hh"=0, "hhl"=0,
                      "MDAcov"=0.5, "MDAp_length"=100, "MDArad_cure"=0.5, "N"=10000,"tau"=5, "nu"=5, "iota"=50/7/10000, "eta"=1,"kappa"=0.18)

  parameters_mda_0c=parameters_mda; parameters_mda_0c$MDAcov=0
  parameters_mda_0p=parameters_mda; parameters_mda_0p$MDAp_length=1
  parameters_mda_0r=parameters_mda; parameters_mda_0r$MDArad_cure=0

  mySTOmodel=model_sto_vivax_delay_rcd_no_referral()
  mySTOmodel_mda=model_sto_vivax_delay_rcd_no_referral_mda()

  mda_cp=simulate_vivax_delay_mda_sto(parameters=parameters_mda , STOmodel =mySTOmodel, STOmodel_mda  =mySTOmodel_mda, maxtime=1000, year=T, runs =10)
  mda_cp_det=simulate_vivax_delay_mda_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=1000, year=T)

  mda_0p=simulate_vivax_delay_mda_sto(parameters=parameters_mda_0p , STOmodel =mySTOmodel, STOmodel_mda  =mySTOmodel_mda, maxtime=1000, year=T, runs=10)
  mda_0p_det=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0p , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=1000, year=T)

  mda_0r=simulate_vivax_delay_mda_sto(parameters=parameters_mda_0r , STOmodel =mySTOmodel, STOmodel_mda  =mySTOmodel_mda, maxtime=1000, year=T, runs=10)
  mda_0r_det=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0r , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=1000, year=T)

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

  expect_equal(mda_cp_mean, mda_cp_det[names(mda_cp_mean)], tolerance = 3e-02, label = "sto and det are equal, with MDA")
  expect_equal(mda_0p_mean, mda_0p_det[names(mda_0p_mean)], tolerance = 2e-02, label = "sto and det are equal, with MDA, zero prophylaxis")
  expect_equal(mda_0r_mean, mda_0r_det[names(mda_0r_mean)], tolerance = 2e-02, label = "sto and det are equal, with MDA, zero radical cure")

  # test if min and max of stochastic simlations are lower and upper bound for deterministic run (sharper than previous test)
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
  #
})


test_that("test simulation of future scenarios, with MDA, RCD (non referral) and delay", {

  mydata=data.frame(incidence=c(273,212),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0.1,0)
  mydata$omega=c(0.7,1)
  mydata$N=c(10000,20000)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata_delay(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA,
             "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA, "rho.new"=NA, "kappa.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA,
             "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  int_0M=list(intervention_name="MDA","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=100, "MDArad_cure.new"=0,
              "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA, "rho.new"=NA, "kappa.new"=NA)
  int_AM=list(intervention_name="A+MDA","rho.new"=NA, "alpha.new"=0.22, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=100, "MDArad_cure.new"=0,
              "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  my_intervention_list=list(int_0,int_A, int_0M,int_AM)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*1, delay=T,mda = F, sto=T, runs=10, rcd=T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*5, delay=T, mda = T, sto=T, runs=10, rcd=T)

  simul2_mean=simul2 %>% dplyr::group_by(.data$time, .data$id, .data$intervention) %>% dplyr::summarise_all("mean") %>%
    data.frame()%>% dplyr::select(-run)%>%
    dplyr::left_join(mydata[c("id", "N")]) %>%
    dplyr::mutate(Il=.data$Il/.data$N, I0=.data$I0/.data$N,  Sl=.data$Sl/.data$N, S0=.data$S0/.data$N, Tl=.data$Tl/.data$N, T0=.data$T0/.data$N,
                  I=.data$I/.data$N, hh=.data$hh/.data$N,hhl=.data$hhl/.data$N,h=.data$h/.data$N,hl=.data$hl/.data$N)

  row.names(simul2_mean)=NULL

  simul1_det=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                          f=f, gamma=gamma, r=r, year=F, maxtime = 365*1, delay=T,mda = F, rcd=T)

  simul2_det=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1_det,
                                          f=f, gamma=gamma, r=r, year=F ,maxtime = 365*5, delay=T, mda = T, rcd=T)
  row.names(simul2_det)=NULL

  expect_equal(simul2_mean  %>% dplyr::arrange(id, intervention, time)%>% dplyr::select(-N, -p , - incidence),
               simul2_det  %>% dplyr::arrange(id, intervention, time)%>% dplyr::select(time, id, intervention, Il, I0, Sl, S0, Tl, T0, hh, hhl, h, hl, I, step ),
               tolerance = 4e-02, label = "MDA sto and non sto")

  # ggplot(simul2 )+
  #   geom_line(aes(x=time, y=I/20000, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul2_det,aes(x=time, y=I, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)
  #
  # ggplot(simul2 )+
  #   geom_line(aes(x=time, y=incidence, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul2_det,aes(x=time, y=incidence, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)

})



test_that("test simulation of future scenarios, with MDA, RCD (non referral) and delay", {

  mydata=data.frame(incidence=c(273,212),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0.1,0)
  mydata$omega=c(0.7,1)
  mydata$N=c(10000,20000)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata_delay(mydata,f=f, gamma=gamma, r=r, return.all = T )

  int_0=list(intervention_name="baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA,
             "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA, "rho.new"=NA, "kappa.new"=NA)
  int_A=list(intervention_name="A", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA,
             "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  int_0M=list(intervention_name="MDA","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=100, "MDArad_cure.new"=0,
              "iota.new"=NA, "nu.new"=NA, "eta.new"=NA, "tau.new"=NA, "rho.new"=NA, "kappa.new"=NA)
  int_AM=list(intervention_name="A+MDA","rho.new"=NA, "alpha.new"=0.22, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.3, "MDAp_length.new"=100, "MDArad_cure.new"=0,
              "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=5, "rho.new"=0.3, "kappa.new"=0.3)
  my_intervention_list=list(int_0,int_A, int_0M,int_AM)
  simul1=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*1, delay=T,mda = F, sto=T, runs=10, rcd=T, referral = T)

  simul2=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1,
                                      f=f, gamma=gamma, r=r, year=F ,maxtime = 365*10, delay=T, mda = T, sto=T, runs=10, rcd=T, referral = T)

  simul2_mean=simul2 %>% dplyr::group_by(.data$time, .data$id, .data$intervention) %>% dplyr::summarise_all("mean") %>%
    data.frame()%>% dplyr::select(-run)%>%
    dplyr::left_join(mydata[c("id", "N")]) %>%
    dplyr::mutate(Il=.data$Il/.data$N, I0=.data$I0/.data$N,  Sl=.data$Sl/.data$N, S0=.data$S0/.data$N, Tl=.data$Tl/.data$N, T0=.data$T0/.data$N,
                  I=.data$I/.data$N, hh=.data$hh/.data$N,hhl=.data$hhl/.data$N,h=.data$h/.data$N,hl=.data$hl/.data$N)

  row.names(simul2_mean)=NULL

  simul1_det=simulate_vivax_interventions(df=mydata2, intervention_list = my_intervention_list,
                                          f=f, gamma=gamma, r=r, year=F, maxtime = 365*1, delay=T,mda = F, rcd=T, referral = T)

  simul2_det=simulate_vivax_interventions(df=mydata2, intervention_list=my_intervention_list, previous_simulation=simul1_det,
                                          f=f, gamma=gamma, r=r, year=F ,maxtime = 365*10, delay=T, mda = T, rcd=T, referral = T)
  row.names(simul2_det)=NULL

  expect_equal(simul2_mean  %>% dplyr::arrange(id, intervention, time)%>% dplyr::select(-N, -p , - incidence),
               simul2_det  %>% dplyr::arrange(id, intervention, time)%>% dplyr::select(time, id, intervention, Il, I0, Sl, S0, Tl, T0, hh, hhl, h, hl, I, step ),
               tolerance = 4e-02, label = "MDA sto and non sto")

  # ggplot(simul2 )+
  #   geom_line(aes(x=time, y=I/2000, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul2_det,aes(x=time, y=I, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)
  #
  # ggplot(simul2 )+
  #   geom_line(aes(x=time, y=incidence, color=as.factor(intervention), group=as.factor(run)))+
  #   geom_line(data=simul2_det,aes(x=time, y=incidence, group=as.factor(intervention)), color="black")+
  #   facet_wrap(intervention~id)

})


test_that("test simulation of future scenarios, with RCD incl. time varying tau and delay", {

  mydata=data.frame(incidence=c(323,312),id=c(1,2))
  mydata$h=incidence_year2day(mydata$incidence)
  mydata$rho=c(0.18,0.13)
  mydata$beta=c(0.43,0.42)
  mydata$alpha=c(0.17, 0.12)
  mydata$sigma=c(1/15, 1/15)
  mydata$prop_import=c(0.01,0.01)
  mydata$omega=c(1,1)
  mydata$N=c(1000,1000)
  f=1/72
  gamma=1/223
  r=1/60
  mydata2=calculate_r0_rc_fromdata_delay(mydata,f=f, gamma=gamma, r=r, return.all = T )

  my_tau=function(pr){
    return(varying_tau(nu=5, pr=pr, N=10000))
  }

  int_B=list(intervention_name="B","rho.new"=NA, "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA,"MDAcov.new"=NA, "MDAp_length.new"=NA, "MDArad_cure.new"=NA, "iota.new"=5/7/10000, "nu.new"=5, "eta.new"=1, "tau.new"=my_tau, "kappa.new"=0.18)
  expect_error(simulate_vivax_interventions(df=mydata2, intervention_list = list(int_B),
                                      f=f, gamma=gamma, r=r, year=F, maxtime = 365*1,mda = F, rcd=T, delay = T, referral = F, sto=T))
  expect_error(simulate_vivax_interventions(df=mydata2, intervention_list = list(int_B),
                                            f=f, gamma=gamma, r=r, year=F, maxtime = 365*1,mda = T, rcd=T, delay = T, referral = F, sto=T))

  })
