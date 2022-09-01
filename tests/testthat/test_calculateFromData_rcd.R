test_that("test calculation of lambda on data, with RCD", {

  mydata=data.frame(incidence=c(12, 45,3,12, 45,3),
                    prop_import=c(0,0,0, 0.1,0.3,0.02)) %>%
    dplyr::mutate(h=incidence_year2day(incidence))

  expect_error(calculate_r0_rc_fromdata_rcd(df=mydata, f=1/72, gamma=1/223, r=1/60))

  mydata=data.frame(incidence=c(171.2512, 114.1977,25.54722,75.27592, 60.47808,46.25247),
                    prop_import=c(0,0,0, 0.1917024 ,0.2389345,0.3130773),
                    alpha=c(0,0,0,0.4,0.4,0.4),
                    beta=c(1,1,1,0.7,0.7,0.7),
                    rho=c(1,1,1,0.4,0.4,0.4),
                    omega=c(1,1,1,0.8,0.8,0.8),
                    iota=c(5/7/365,10000, Inf, 5/7/365,Inf, Inf),
                    nu=c(5,5,5,5,5,5), tau=c(5,5,5,5,5,5),eta=c(1,0.9,1,1,0.9,1))%>%
    dplyr::mutate(h=incidence_year2day(incidence))

  new_data_all =calculate_r0_rc_fromdata_rcd(df=mydata, f=1/72, gamma=1/223, r=1/60, return.all = T)

  lambda_true=c(0.015, 0.012, 0.0082, 0.015, 0.012, 0.0082)
  expect_equal(new_data_all$lambda, lambda_true ,tolerance = 1e-06, label = "lambda")

})

test_that("test calculation of lambda on data, with RCD and delays", {

  mydata=data.frame(incidence=c(12, 45,3,12, 45,3),
                    prop_import=c(0,0,0, 0.1,0.3,0.02)) %>%
    dplyr::mutate(h=incidence_year2day(incidence))

  expect_error(calculate_r0_rc_fromdata_delay_rcd(df=mydata, f=1/72, gamma=1/223, r=1/60))
  expect_error(calculate_r0_rc_fromdata_delay_rcd(df=mydata, f=1/72, gamma=1/223, r=1/60, referral = TRUE))

  mydata=data.frame(incidence=c(127.1534    , 74.20737,16.7957   ,94.70294  , 73.07549 ,52.75153 ),
                    prop_import=c(0,0,0, 0.1518257   ,0.1971761  ,0.2739934 ),
                    alpha=c(0.1,0.1,0.1,0.4,0.4,0.4),
                    beta=c(1,1,1,0.7,0.7,0.7),
                    rho=c(1,1,1,0.4,0.4,0.4),
                    omega=c(1,1,1,0.8,0.8,0.8),
                    sigma=c(1/15,1/15,1/15,1/10,1/10,1/10),
                    iota=c(5/7/365,10000, Inf, 5/7/365,Inf, Inf),
                    nu=c(5,5,5,5,5,5), tau=c(5,5,5,5,5,5),eta=c(1,0.9,1,1,0.9,1))%>%
    dplyr::mutate(h=incidence_year2day(incidence))

  mydata_ref=data.frame(incidence=c(187.1244   , 103.0747,21.98137   ,102.7095   , 76.89312 ,54.29549 ),
                        prop_import=c(0,0,0, 0.1397811     ,0.1872207  ,0.2660791 ),
                        alpha=c(0.1,0.1,0.1,0.4,0.4,0.4),
                        beta=c(1,1,1,0.7,0.7,0.7),
                        rho=c(1,1,1,0.4,0.4,0.4),
                        omega=c(1,1,1,0.8,0.8,0.8),
                        sigma=c(1/15,1/15,1/15,1/10,1/10,1/10),
                        iota=c(5/7/365,10000, Inf, 5/7/365,Inf, Inf),
                        nu=c(5,5,5,5,5,5), tau=c(5,5,5,5,5,5),eta=c(1,0.9,1,1,0.9,1))%>%
    dplyr::mutate(h=incidence_year2day(incidence))


  new_data_all =calculate_r0_rc_fromdata_delay_rcd(df=mydata, f=1/72, gamma=1/223, r=1/60, return.all = T)
  new_data_all_ref =calculate_r0_rc_fromdata_delay_rcd(df=mydata_ref, f=1/72, gamma=1/223, r=1/60, return.all = T, referral = T)

  lambda_true=c(0.015, 0.012, 0.0092, 0.015, 0.012, 0.0082)
  expect_equal(new_data_all$lambda, lambda_true ,tolerance = 1e-06, label = "lambda")
  expect_equal(new_data_all_ref$lambda, lambda_true ,tolerance = 1e-06, label = "lambda")

})

test_that("test simulation of future scenarios, with rcd", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      alpha.old=0.95*rho,
      lambda=c(0.008  ,0.008  ,0.03 ),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3),
      omega.old=c(0.7,0.7,0.7),
      nu.new=c(0,10,5), tau.new=c(5,5,5), eta.new=c(1,1,1), kappa.new=rho)

  expect_error(simulate_from_data(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=T, rcd=T),
               "RCD parameters are missing, please add them or use model without RCD", label = "error when missing RCD parameter")

  mydata=mydata %>% dplyr::mutate(iota.new=c(5/7/10000,5/7/10000,5/7/10000))

  simul_rcd=simulate_from_data(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=T, rcd=T)
  simul_no_rcd=simulate_from_data(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=T, rcd=F)

  expect_equal(simul_rcd[simul_rcd$id==1,], simul_no_rcd[simul_no_rcd$id==1,],tolerance = 1e-07, label = "no rcd is the same as rcd with nu=0")

  expect_lt(simul_rcd$I[simul_rcd$id==2 & simul_rcd$time==1825], simul_no_rcd$I[simul_no_rcd$id==2 & simul_no_rcd$time==1825], label = "rcd is better than no rcd")
  expect_lt(simul_rcd$I[simul_rcd$id==3 & simul_rcd$time==1825], simul_no_rcd$I[simul_no_rcd$id==3 & simul_no_rcd$time==1825], label = "rcd is better than no rcd")

  simul_rcd_true=data.frame(time=rep(1825, 3),
                            Il=c(0.02588585 ,0.07185492 ,0.48472036 ), I0=c(0.003978321 ,0.008957964 ,0.040390022 ),
                            Sl=c(0.02813485 ,0.08285706 ,0.30146630 ), S0=c(0.9420010 ,0.8363301 ,0.1734233 ),
                            h=c(0.03885228  ,0.19450692     ,0.64846559     ), hr=c(0.02615755 ,0.05721143,0.12757704  ),
                            hl=c(0.03657982   ,0.18275685    ,0.62274985    ),
                            hh=c(0.215846  ,0.682824 ,3.827577   ),hhl=c(0.2032212  ,0.5924389 ,3.5061300   ),
                            I=c(0.02986417,0.08081288 ,0.52511039  ), p=c(0.05848969      ,0.06040952    , 0.03965630     ),
                            incidence=c(38.85228  , 194.50692     ,648.46559      ), id=c(1,2,3))   #it's ok if we are not stable for id=1 because lambda, delta and omega do not correspond to an equilibrium

  row.names(simul_rcd_true)=c(1826, 18261, 18262)

  expect_equal(simul_rcd[simul_rcd$time==1825,], simul_rcd_true,tolerance = 1e-06, label = "check if the values are the same as pre-computed ones, with RCD")

})


test_that("test simulation of future scenarios starting from initial condition, at equilibrium (with 0 RCD)", {

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
      omega.old=c(0.7,0.7,0.7),
      nu.new=c(0,10,5), tau.new=c(5,5,5), eta.new=c(1,1,0),iota.new=c(5/7/10000,0,5/7/10000), kappa.new=rho)

  simul_rcd=simulate_from_data(mydata, f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=F, rcd=T)

  my_initial_state=simul_rcd[simul_rcd$time==2000,c("Il", "I0","Sl", "S0","h", "hr", "hl", "hh", "hhl", "id")]
  names(my_initial_state)=c("Il_init", "I0_init","Sl_init", "S0_init","h_init", "hr_init", "hl_init", "hh_init", "hhl_init", "id")
  simul_rcd_chain=simulate_from_data(df=mydata, from_equilibrium = FALSE,
                                        initial_states = my_initial_state,
                                        f=1/69, gamma=1/383, r=1/60,  maxtime=2000,year=F, rcd=T)

  expect_equal(simul_rcd, simul_rcd_chain,tolerance = 1e-07, label = "chaining from equilibrium is the same as simulating from equilibrium, with rcd")

})

test_that("test simulation of future scenarios starting from initial condition, not at equilibrium", {

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
      omega.old=c(0.7,0.7,0.7),
      nu.new=c(0,10,5), tau.new=c(5,5,5), eta.new=c(1,1,1),iota.new=c(5/7/10000,5/7/10000,5/7/10000), kappa.new=rho)

  simul_rcd=simulate_from_data(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=1000,year=F, rcd=T)

  my_initial_state=simul_rcd[simul_rcd$time==100,c("Il", "I0","Sl", "S0","h", "hr", "hl", "hh", "hhl", "id")]
  names(my_initial_state)=c("Il_init", "I0_init","Sl_init", "S0_init","h_init", "hr_init", "hl_init", "hh_init", "hhl_init", "id")
  simul_rcd_chain=simulate_from_data(df=mydata, from_equilibrium = FALSE,
                                        initial_states = my_initial_state,
                                        f=1/69, gamma=1/383, r=1/60,  maxtime=900,year=F, rcd=T) %>%
    dplyr::mutate(time=time+100)


  expect_equal(simul_rcd %>% dplyr::filter(time>=100), simul_rcd_chain,tolerance = 1e-07,
               label = "chaining does not change the trajectory, with rcd")


  simul_rcd_chain_newalpha=simulate_from_data(df=mydata %>% dplyr::mutate(alpha.new=0.5), from_equilibrium = FALSE,
                                                 initial_states = my_initial_state,
                                                 f=1/69, gamma=1/383, r=1/60,  maxtime=900,year=F, rcd=T) %>%
    dplyr::mutate(time=time+100)

  expect_lt(simul_rcd_chain_newalpha$I[900], simul_rcd_chain$I[900],
            label = "adding an intervention changes the trajectory, with rcd")


})

###########
# WITH DELAYS

test_that("test simulation of future scenarios, with rcd, with delays", {
  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      TQ_effect=c(0.59,0.605,0.619),
      alpha.old=0.95*rho,
      sigma.old=c(1/15  ,1/15   ,1/15  ),
      lambda=c(0.008  ,0.008  ,0.03 ),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3),
      omega.old=c(0.7,0.7,0.7),
      nu.new=c(0,10,5), tau.new=c(5,5,5), eta.new=c(1,1,1), kappa.new=c(0.18,0.18,0.18))

  expect_error(simulate_from_data_delay(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=T, rcd=T),
               "RCD parameters are missing, please add them or use model without RCD", label = "error when missing RCD parameter")

  mydata=mydata %>% dplyr::mutate(iota.new=c(5/7/10000,5/7/10000,5/7/10000))

  simul_rcd=simulate_from_data_delay(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=T, rcd=T)
  simul_rcd_ref=simulate_from_data_delay(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=T, rcd=T, referral = T)
  simul_no_rcd=simulate_from_data_delay(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=T, rcd=F)

  expect_equal(simul_rcd[simul_rcd$id==1,], simul_no_rcd[simul_no_rcd$id==1,],tolerance = 1e-07, label = "no rcd is the same as rcd with nu=0")
  expect_equal(simul_rcd[simul_rcd$id==1,], simul_no_rcd[simul_no_rcd$id==1,],tolerance = 1e-07, label = "no rcd is the same as rcd with nu=0")

  expect_lt(simul_rcd$I[simul_rcd$id==2 & simul_rcd$time==1825], simul_no_rcd$I[simul_no_rcd$id==2 & simul_no_rcd$time==1825], label = "rcd is better than no rcd")
  expect_lt(simul_rcd$I[simul_rcd$id==3 & simul_rcd$time==1825], simul_no_rcd$I[simul_no_rcd$id==3 & simul_no_rcd$time==1825], label = "rcd is better than no rcd")

  expect_lt(simul_rcd_ref$I[simul_rcd_ref$id==2 & simul_rcd_ref$time==1825], simul_no_rcd$I[simul_no_rcd$id==2 & simul_no_rcd$time==1825], label = "rcd ref is better than no rcd")
  expect_lt(simul_rcd_ref$I[simul_rcd_ref$id==3 & simul_rcd_ref$time==1825], simul_no_rcd$I[simul_no_rcd$id==3 & simul_no_rcd$time==1825], label = "rcd ref is better than no rcd")

  expect_lt(simul_rcd$I[simul_rcd$id==2 & simul_rcd$time==1825], simul_rcd_ref$I[simul_rcd_ref$id==2 & simul_no_rcd$time==1825], label = "no referral is better than referral")
  expect_lt(simul_rcd$I[simul_rcd$id==3 & simul_rcd$time==1825], simul_rcd_ref$I[simul_rcd_ref$id==3 & simul_no_rcd$time==1825], label = "no referral is better than referral")

  simul_rcd_true=data.frame(time=rep(1825, 3),
                            Ul=c(0.03164403  ,0.07705805  ,0.48330159  ), U0=c(0.004831643  ,0.009582254  ,0.040062752  ),
                            Sl=c(0.03465808  ,0.08950477  ,0.30071105  ), S0=c(0.9273484  ,0.8208915  ,0.1663927  ),
                            Tl=c(0.001471987   ,0.002874158   ,0.009280484   ), T0=c(4.588674e-05  ,8.923783e-05  ,2.514149e-04  ),
                            h=c(0.04706137   ,0.20836941    ,0.64630938     ),
                            hl=c(0.04480642    ,0.19673111      ,0.62101525  ),
                            hh=c(0.2614521   ,0.7316038   ,3.8148497   ),
                            hhl=c(0.2489246    ,0.6420784    ,3.4986731    ),
                            I=c(0.03799354 ,0.08960370  ,0.53289624   ), p=c(0.04791506      ,0.05585417      , 0.03913626      ),
                            incidence=c(47.06137    , 208.36941  ,646.30938   ), id=c(1,2,3)) #it's ok if we are not stable for id=1 because lambda, delta and omega do not correspond to an equilibrium

  row.names(simul_rcd_true)=c(1826, 18261, 18262)

  expect_equal(simul_rcd[simul_rcd$time==1825,], simul_rcd_true,tolerance = 1e-06, label = "check if the values are the same as pre-computed ones, with RCD")

  simul_rcd_ref_true=data.frame(time=rep(1825, 3),
                            Ul=c(0.03164403  ,0.08379253   ,0.48151826   ), U0=c(0.004831643  ,0.010384723   ,0.039677877   ),
                            Sl=c(0.03465808  ,0.09805748   ,0.29984063   ), S0=c(0.9273484  ,0.8005070   ,0.1583023   ),
                            Tl=c(0.001471987   ,0.006612432    ,0.019402570    ), T0=c(4.588674e-05  ,6.458097e-04  ,1.258329e-03  ),
                            h=c(0.04706137  ,0.22630157       ,0.64363184      ),
                            hl=c(0.04480642   ,0.21481370        ,0.61882294        ),
                            hh=c(0.2614521   ,0.7947005    ,3.7990455    ),
                            hhl=c(0.2489246    ,0.7063323     ,3.4889343  ),
                            I=c(0.03799354 ,0.10143549   ,0.54185704    ), p=c(0.04791506        ,0.05076355      , 0.03854517       ),
                            incidence=c(47.06137    , 226.30157    ,643.63184              ), id=c(1,2,3))

  row.names(simul_rcd_ref_true)=c(1826, 18261, 18262)
  expect_equal(simul_rcd_ref[simul_rcd_ref$time==1825,], simul_rcd_ref_true,tolerance = 1e-06, label = "check if the values are the same as pre-computed ones, with RCD referral")

})


test_that("test simulation of future scenarios starting from initial condition, with delays, at equilibrium (with 0 RCD)", {

  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      TQ_effect=c(0.59,0.605,0.619),
      alpha.old=0.95*rho,
      sigma.old=c(1/15  ,1/15   ,1/15  ),
      lambda=c(0.0057453432   ,0.0091275728 ,0.0167897768  ),
      I=c(0.018131142 ,0.127630432 ,0.515273425 ),
      delta=c(0.000035654059 ,5.4114206e-05,3.772780814e-03),
      id=c(1,2,3),
      omega.old=c(0.7,0.7,0.7),
      nu.new=c(0,10,5), tau.new=c(5,5,5), eta.new=c(1,1,0), kappa.new=rho,iota.new=c(5/7/10000,0,5/7/10000))

  simul_rcd=simulate_from_data_delay(mydata, f=1/69, gamma=1/383, r=1/60,      maxtime=2000,year=F, rcd=T)

  my_initial_state=simul_rcd[simul_rcd$time==2000,c("Ul", "U0","Sl", "S0","Tl", "T0","h", "hl","hh", "hhl", "id")]
  names(my_initial_state)=c("Ul_init", "U0_init","Sl_init", "S0_init","Tl_init", "T0_init","h_init","hl_init","hh_init","hhl_init", "id")
  simul_rcd_chain=simulate_from_data_delay(df=mydata, from_equilibrium = FALSE,
                                     initial_states = my_initial_state,
                                     f=1/69, gamma=1/383, r=1/60,  maxtime=2000,year=F, rcd=T)

  expect_equal(simul_rcd, simul_rcd_chain,tolerance = 1e-07, label = "chaining from equilibrium is the same as simulating from equilibrium, with delay, with rcd")

  simul_rcd_chain_ref=simulate_from_data_delay(df=mydata, from_equilibrium = FALSE,
                                           initial_states = my_initial_state,
                                           f=1/69, gamma=1/383, r=1/60,  maxtime=2000,year=F, rcd=T, referral = T)

  expect_equal(simul_rcd, simul_rcd_chain_ref,tolerance = 1e-07, label = "chaining from equilibrium is the same as simulating from equilibrium, with delay, with rcd, with referral")

})

test_that("test simulation of future scenarios starting from initial condition, not at equilibrium, with delays, with rcd", {

  mydata=data.frame(incidence=c(23,112,267)) %>%
    dplyr::mutate(
      rho=c(0.18,0.13,0.08) ,
      rho.old=rho,
      beta.old=c(0.431,0.429,0.422),
      TQ_effect=c(0.59,0.605,0.619),
      alpha.old=0.95*rho,
      sigma.old=c(1/15  ,1/15   ,1/15  ),
      lambda=c(0.006260546  ,0.007086074  ,0.022004378 ),
      I=c(0.01741279 ,0.12413235 ,0.50693425 ),
      delta=c(3.562799e-05,2.694904e-04,1.854486e-03),
      id=c(1,2,3),
      omega.old=c(0.7,0.7,0.7),
      nu.new=c(0,10,5), tau.new=c(5,5,5), eta.new=c(1,1,1), kappa.new=rho,iota.new=c(5/7/10000,5/7/10000,5/7/10000))

  simul_rcd=simulate_from_data_delay(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=1000,year=F, rcd=T)

  my_initial_state=simul_rcd[simul_rcd$time==100, c("Ul", "U0","Sl", "S0","Tl", "T0","h","hl","hh","hhl", "id")]
  names(my_initial_state)=c("Ul_init", "U0_init","Sl_init", "S0_init","Tl_init", "T0_init","h_init", "hl_init","hh_init", "hhl_init","id")
  simul_rcd_chain=simulate_from_data_delay(df=mydata, from_equilibrium = FALSE,
                                     initial_states = my_initial_state,
                                     f=1/69, gamma=1/383, r=1/60,  maxtime=900,year=F, rcd=T) %>%
    dplyr::mutate(time=time+100)


  expect_equal(simul_rcd %>% dplyr::filter(time>=100), simul_rcd_chain,tolerance = 1e-07,
               label = "chaining does not change the trajectory, with delay, with rcd")


  simul_rcd_chain_newalpha=simulate_from_data_delay(df=mydata %>% dplyr::mutate(alpha.new=0.5), from_equilibrium = FALSE,
                                              initial_states = my_initial_state,
                                              f=1/69, gamma=1/383, r=1/60,  maxtime=900,year=F, rcd=T) %>%
    dplyr::mutate(time=time+100)

  expect_lt(simul_rcd_chain_newalpha$I[900], simul_rcd_chain$I[900],
            label = "adding an intervention changes the trajectory, with delay, with rcd")


  ### with referral ###

  simul_rcd_ref=simulate_from_data_delay(mydata,    f=1/69, gamma=1/383, r=1/60,      maxtime=1000,year=F, rcd=T, referral = T)

  my_initial_state_ref=simul_rcd_ref[simul_rcd_ref$time==100, c("Ul", "U0","Sl", "S0","Tl", "T0","h", "hl","hh", "hhl", "id")]
  names(my_initial_state_ref)=c("Ul_init", "U0_init","Sl_init", "S0_init","Tl_init", "T0_init","h_init","hl_init","hh_init", "hhl_init", "id")
  simul_rcd_chain_ref=simulate_from_data_delay(df=mydata, from_equilibrium = FALSE,
                                           initial_states = my_initial_state_ref,
                                           f=1/69, gamma=1/383, r=1/60,  maxtime=900,year=F, rcd=T, referral = T) %>%
    dplyr::mutate(time=time+100)


  expect_equal(simul_rcd_ref %>% dplyr::filter(time>=100), simul_rcd_chain_ref,tolerance = 1e-07,
               label = "chaining does not change the trajectory, with delay, with rcd referral")


  simul_rcd_chain_newalpha_ref=simulate_from_data_delay(df=mydata %>% dplyr::mutate(alpha.new=0.5), from_equilibrium = FALSE,
                                                    initial_states = my_initial_state_ref,
                                                    f=1/69, gamma=1/383, r=1/60,  maxtime=900,year=F, rcd=T, referral = T) %>%
    dplyr::mutate(time=time+100)

  expect_lt(simul_rcd_chain_newalpha_ref$I[900], simul_rcd_chain_ref$I[900],
            label = "adding an intervention changes the trajectory, with delay, with rcd referral")



})


