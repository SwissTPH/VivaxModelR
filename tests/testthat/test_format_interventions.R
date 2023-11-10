

test_that("test formating data for intervention simulation", {

  df0=data.frame(id=c("regionA","regionB"))
  df=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), rho=c(0.1,0.2))

  interv0=list(intervention_name="0", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA , "rho.new"=NA)
  intervA=list(intervention_name="A", "alpha.new"=0.7, "beta.new"=0.8, "omega.new"=0.6, "rho.new"=0.7 )

  expectdata0=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), rho=c(0.1,0.2),
                         rho.old=c(0.1, 0.2), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.1,0.2), alpha.new=c(0.17, 0.12), beta.new=c(0.43,0.42), omega.new=c(0.18, 0.13),
                         intervention="0")

  expectdata0_0=data.frame(id=c("regionA","regionB"), rho=c(1,1),alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                           rho.old=c(1, 1), alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                           rho.new=c(1,1), alpha.new=c(0,0), beta.new=c(1,1), omega.new=c(1, 1),
                           intervention="0")

  expectdataA=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), rho=c(0.1,0.2),
                         rho.old=c(0.1, 0.2),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.7, 0.7),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                         intervention="A")

  expectdataA_0=data.frame(id=c("regionA","regionB"), rho=c(1, 1), alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                           rho.old=c(1, 1),alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                           rho.new=c(0.7, 0.7),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                           intervention="A")


  expect_equal(format_data_simulation(df, interv0), expectdata0, label = " no intervention")
  expect_equal(format_data_simulation(df, intervA), expectdataA, label = " intervention A")
  expect_equal(format_data_simulation(df0, interv0), expectdata0_0, label = " no intervention default")
  expect_equal(format_data_simulation(df0, intervA), expectdataA_0, label = " intervention A default")
})



test_that("test formating data for intervention simulation, with delay", {

  df0=data.frame(id=c("regionA","regionB"))
  df0d=data.frame(id=c("regionA","regionB"), sigma=c(1/15, 1/15))
  df=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.3,0.4))

  interv0=list(intervention_name="0", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "rho.new"=NA )
  intervA=list(intervention_name="A", "alpha.new"=0.7, "beta.new"=0.8, "omega.new"=0.6 , "sigma.new"=1/5, "rho.new"=0.7)

  expectdata0=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.3,0.4),
                         rho.old=c(0.3,0.4), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.3,0.4), alpha.new=c(0.17, 0.12), beta.new=c(0.43,0.42), omega.new=c(0.18, 0.13),
                         intervention="0")

  expectdata0_d=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.3,0.4),
                           rho.old=c(0.3,0.4), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                           rho.new=c(0.3,0.4), alpha.new=c(0.17, 0.12), beta.new=c(0.43,0.42), omega.new=c(0.18, 0.13), sigma.old=c(1/15, 1/15), sigma.new=c(1/15, 1/15),
                           intervention="0")

  expectdata0_0=data.frame(id=c("regionA","regionB"), rho=c(1,1), alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                           rho.old=c(1,1), alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                           rho.new=c(1,1), alpha.new=c(0,0), beta.new=c(1,1), omega.new=c(1, 1),
                           intervention="0")

  expectdata0_0d=data.frame(id=c("regionA","regionB"),sigma=c(1/15, 1/15),rho=c(1,1),  alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                            rho.old=c(1,1), alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                            rho.new=c(1,1), alpha.new=c(0,0), beta.new=c(1,1), omega.new=c(1, 1), sigma.old=c(1/15, 1/15), sigma.new=c(1/15, 1/15),
                            intervention="0")

  expectdataA=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.3,0.4),
                         rho.old=c(0.3,0.4),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.7,0.7),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                         intervention="A")

  expectdataAd=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.3,0.4),
                          rho.old=c(0.3,0.4),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                          rho.new=c(0.7,0.7),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6), sigma.old=c(1/15, 1/15), sigma.new=c(1/5, 1/5),
                          intervention="A")

  expectdataA_0=data.frame(id=c("regionA","regionB"), rho=c(1,1),alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                           rho.old=c(1,1), alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                           rho.new=c(0.7,0.7),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                           intervention="A")

  expectdataA_0d=data.frame(id=c("regionA","regionB"), sigma=c(1/15, 1/15), rho=c(1,1),alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                            rho.old=c(1,1), alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                            rho.new=c(0.7,0.7),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6), sigma.old=c(1/15, 1/15), sigma.new=c(1/5, 1/5),
                            intervention="A")


  expect_equal(format_data_simulation(df, interv0, delay = F), expectdata0, label = " no intervention no delay")
  expect_equal(format_data_simulation(df, interv0, delay = T), expectdata0_d, label = " no intervention delay")

  expect_equal(format_data_simulation(df, intervA, delay=F), expectdataA, label = " intervention A no delay")
  expect_equal(format_data_simulation(df, intervA, delay=T), expectdataAd, label = " intervention A delay")

  expect_equal(format_data_simulation(df0, interv0, delay=F), expectdata0_0, label = " no intervention default no delay")
  expect_equal(format_data_simulation(df0d, interv0, delay=T), expectdata0_0d, label = " no intervention default delay")

  expect_equal(format_data_simulation(df0, intervA, delay=F), expectdataA_0, label = " intervention A default")
  expect_equal(format_data_simulation(df0d, intervA, delay=T), expectdataA_0d, label = " intervention A default")
})




test_that("test formating data for intervention simulation, with RCD", {

  df0=data.frame(id=c("regionA","regionB"))
  df=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), rho=c(0.43,0.55))

  interv0=list(intervention_name="0", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA,"rho.new"=NA,
               "iota.new"=5/7/10000, "nu.new"=5, "tau.new"=2, "eta.new"=0.5 )
  intervA=list(intervention_name="A", "alpha.new"=0.7, "beta.new"=0.8, "omega.new"=0.6,"rho.new"=0.8,
               "iota.new"=5/7/10000, "nu.new"=5, "tau.new"=2, "eta.new"=0.5)
  intervB=list(intervention_name="B", "alpha.new"=0.7, "beta.new"=0.8, "omega.new"=0.6,"rho.new"=0.8,
               "iota.new"=5/7/10000, "nu.new"=5, "tau.new"=2, "eta.new"=0.5, "rho2.new"=0.82)

  expectdata0=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), rho=c(0.43,0.55),
                         rho.old=c(0.43,0.55), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.43,0.55), alpha.new=c(0.17, 0.12), beta.new=c(0.43,0.42), omega.new=c(0.18, 0.13),
                         iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                         intervention="0")

  expectdata0_0=data.frame(id=c("regionA","regionB"), rho=c(1,1), alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                           rho.old=c(1,1), alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                           rho.new=c(1,1), alpha.new=c(0,0), beta.new=c(1,1), omega.new=c(1, 1),
                           iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                           intervention="0")

  expectdataA=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), rho=c(0.43,0.55),
                         rho.old=c(0.43,0.55),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.8,0.8),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                         iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                         intervention="A")

  expectdataA_0=data.frame(id=c("regionA","regionB"), rho=c(1,1), alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                           rho.old=c(1,1), alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                           rho.new=c(0.8,0.8),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                           iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                           intervention="A")
  expectdataB=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), rho=c(0.43,0.55),
                         rho.old=c(0.43,0.55),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.8,0.8),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                         iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(0.82,0.82),
                         intervention="B")


  expect_equal(format_data_simulation(df, interv0, rcd = T), expectdata0, label = " no intervention, rcd")
  expect_equal(format_data_simulation(df, intervA, rcd = T), expectdataA, label = " intervention A, rcd")
  expect_equal(format_data_simulation(df0, interv0, rcd = T), expectdata0_0, label = " no intervention default, rcd")
  expect_equal(format_data_simulation(df0, intervA, rcd = T), expectdataA_0, label = " intervention A default, rcd")
  expect_equal(format_data_simulation(df, intervB, rcd = T), expectdataB, label = " intervention B, rcd, rho2")
})





test_that("test formating data for intervention simulation, with delay and RCD", {

  df0=data.frame(id=c("regionA","regionB"))
  df0d=data.frame(id=c("regionA","regionB"), sigma=c(1/15, 1/15))
  df=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21))

  interv0=list(intervention_name="0", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA,"rho.new"=NA,
               "iota.new"=5/7/10000, "nu.new"=5, "tau.new"=2, "eta.new"=0.5 )
  intervA=list(intervention_name="A", "alpha.new"=0.7, "beta.new"=0.8, "omega.new"=0.6 , "sigma.new"=1/5,"rho.new"=0.77,
               "iota.new"=5/7/10000, "nu.new"=5, "tau.new"=2, "eta.new"=0.5)
  intervB=list(intervention_name="B", "alpha.new"=0.7, "beta.new"=0.8, "omega.new"=0.6 , "sigma.new"=1/5,"rho.new"=0.77,
               "iota.new"=5/7/10000, "nu.new"=5, "tau.new"=2, "eta.new"=0.5, "rho2.new"=0.82)

  expectdata0=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                         rho.old=c(0.22,0.21), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.22,0.21), alpha.new=c(0.17, 0.12), beta.new=c(0.43,0.42), omega.new=c(0.18, 0.13),
                         iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                         intervention="0")

  expectdata0_d=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                           rho.old=c(0.22,0.21), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                           rho.new=c(0.22,0.21), alpha.new=c(0.17, 0.12), beta.new=c(0.43,0.42), omega.new=c(0.18, 0.13), sigma.old=c(1/15, 1/15), sigma.new=c(1/15, 1/15),
                           iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                           intervention="0")

  expectdata0_0=data.frame(id=c("regionA","regionB"),rho=c(1,1), alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                           rho.old=c(1,1), alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                           rho.new=c(1,1), alpha.new=c(0,0), beta.new=c(1,1), omega.new=c(1, 1),
                           iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                           intervention="0")

  expectdata0_0d=data.frame(id=c("regionA","regionB"),sigma=c(1/15, 1/15), rho=c(1,1),alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                            rho.old=c(1,1), alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                            rho.new=c(1,1), alpha.new=c(0,0), beta.new=c(1,1), omega.new=c(1, 1), sigma.old=c(1/15, 1/15), sigma.new=c(1/15, 1/15),
                            iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                            intervention="0")

  expectdataA=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                         rho.old=c(0.22,0.21), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.77,0.77), alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                         iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                         intervention="A")

  expectdataAd=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                          rho.old=c(0.22,0.21),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                          rho.new=c(0.77,0.77),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6), sigma.old=c(1/15, 1/15), sigma.new=c(1/5, 1/5),
                          iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                          intervention="A")

  expectdataA_0=data.frame(id=c("regionA","regionB"),rho=c(1,1), alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                           rho.old=c(1,1),alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                           rho.new=c(0.77,0.77),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                           iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                           intervention="A")

  expectdataA_0d=data.frame(id=c("regionA","regionB"), sigma=c(1/15, 1/15), rho=c(1,1), alpha=c(0,0), beta=c(1,1), omega=c(1, 1),
                            rho.old=c(1,1),alpha.old=c(0,0), beta.old=c(1,1), omega.old=c(1, 1),
                            rho.new=c(0.77,0.77),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6), sigma.old=c(1/15, 1/15), sigma.new=c(1/5, 1/5),
                            iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                            intervention="A")

  expectdataB=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                         rho.old=c(0.22,0.21), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.77,0.77), alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                         iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(0.82,0.82),
                         intervention="B")

  expectdataBd=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                          rho.old=c(0.22,0.21),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                          rho.new=c(0.77,0.77),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6), sigma.old=c(1/15, 1/15), sigma.new=c(1/5, 1/5),
                          iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(0.82,0.82),
                          intervention="B")
  expect_equal(format_data_simulation(df, interv0, delay = F, rcd=T), expectdata0, label = " no intervention no delay, rcd")
  expect_equal(format_data_simulation(df, interv0, delay = T, rcd=T), expectdata0_d, label = " no intervention delay, rcd")

  expect_equal(format_data_simulation(df, intervA, delay=F, rcd=T), expectdataA, label = " intervention A no delay, rcd")
  expect_equal(format_data_simulation(df, intervA, delay=T, rcd=T), expectdataAd, label = " intervention A delay, rcd")

  expect_equal(format_data_simulation(df0, interv0, delay=F, rcd=T), expectdata0_0, label = " no intervention default no delay, rcd")
  expect_equal(format_data_simulation(df0d, interv0, delay=T, rcd=T), expectdata0_0d, label = " no intervention default delay, rcd")

  expect_equal(format_data_simulation(df0, intervA, delay=F, rcd=T), expectdataA_0, label = " intervention A default, rcd")
  expect_equal(format_data_simulation(df0d, intervA, delay=T, rcd=T), expectdataA_0d, label = " intervention A default, rcd")

  expect_equal(format_data_simulation(df, intervB, delay=F, rcd=T), expectdataB, label = " intervention B no delay, rcd rho2")
  expect_equal(format_data_simulation(df, intervB, delay=T, rcd=T), expectdataBd, label = " intervention B delay, rcd rho2")
})



test_that("test formating data for intervention simulation, with RCD at baseline", {

  df0=data.frame(id=c("regionA","regionB"))
  df0d=data.frame(id=c("regionA","regionB"), sigma=c(1/15, 1/15))
  df=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13),
                sigma=c(1/15, 1/15), rho=c(0.22,0.21),iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8),
                eta=c(0.5, 0.75), tau=c(12, 4)
                )
  dfB=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13),
                sigma=c(1/15, 1/15), rho=c(0.22,0.21),iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8),
                eta=c(0.5, 0.75), tau=c(12, 4), rho2=c(0.7, 1)
  )

  interv0=list(intervention_name="0", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA,"rho.new"=NA,
               "iota.new"=NA, "nu.new"=NA, "tau.new"=NA, "eta.new"=NA )
  intervA=list(intervention_name="A", "alpha.new"=0.7, "beta.new"=0.8, "omega.new"=0.6 , "sigma.new"=1/5,"rho.new"=0.77,
               "iota.new"=5/7/10000, "nu.new"=5, "tau.new"=2, "eta.new"=0.5)
  intervA2=list(intervention_name="A", "alpha.new"=0.7, "beta.new"=0.8, "omega.new"=0.6 , "sigma.new"=1/5,"rho.new"=0.77,
               "iota.new"=5/7/10000, "nu.new"=5, "tau.new"=2, "eta.new"=0.5, "rho2.new"=NA)
  intervB=list(intervention_name="B", "alpha.new"=0.7, "beta.new"=0.8, "omega.new"=0.6 , "sigma.new"=1/5,"rho.new"=0.77,
               "iota.new"=5/7/10000, "nu.new"=5, "tau.new"=2, "eta.new"=0.5, "rho2.new"=0.85)

  expectdata0=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                         iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8), eta=c(0.5, 0.75), tau=c(12,4),
                         rho.old=c(0.22,0.21), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.22,0.21), alpha.new=c(0.17, 0.12), beta.new=c(0.43,0.42), omega.new=c(0.18, 0.13),rho2=c(1,1),
                         iota.old=c(5/7/1000, 1000/7/1000), nu.old=c(5,8), tau.old=c(12,4), eta.old=c(0.5,0.75),rho2.old=c(1,1),
                         iota.new=c(5/7/1000,1000/7/1000), nu.new=c(5,8), tau.new=c(12,4), eta.new=c(0.5,0.75),rho2.new=c(1,1),
                         intervention="0")

  expectdata0_d=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                           iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8), eta=c(0.5, 0.75), tau=c(12,4),
                           rho.old=c(0.22,0.21), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                           rho.new=c(0.22,0.21), alpha.new=c(0.17, 0.12), beta.new=c(0.43,0.42), omega.new=c(0.18, 0.13), sigma.old=c(1/15, 1/15), sigma.new=c(1/15, 1/15), rho2=c(1,1),
                           iota.old=c(5/7/1000, 1000/7/1000), nu.old=c(5,8), tau.old=c(12,4), eta.old=c(0.5,0.75),rho2.old=c(1,1),
                           iota.new=c(5/7/1000,1000/7/1000), nu.new=c(5,8), tau.new=c(12,4), eta.new=c(0.5,0.75),rho2.new=c(1,1),
                           intervention="0")

  expectdataA=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                         iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8), eta=c(0.5, 0.75), tau=c(12,4),
                         rho.old=c(0.22,0.21), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.77,0.77), alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),rho2=c(1,1),
                         iota.old=c(5/7/1000, 1000/7/1000), nu.old=c(5,8), tau.old=c(12,4), eta.old=c(0.5,0.75),rho2.old=c(1,1),
                         iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                         intervention="A")

  expectdataAd=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                          iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8), eta=c(0.5, 0.75), tau=c(12,4),
                          rho.old=c(0.22,0.21),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                          rho.new=c(0.77,0.77),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6), sigma.old=c(1/15, 1/15), sigma.new=c(1/5, 1/5),rho2=c(1,1),
                          iota.old=c(5/7/1000, 1000/7/1000), nu.old=c(5,8), tau.old=c(12,4), eta.old=c(0.5,0.75),rho2.old=c(1,1),
                          iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(1,1),
                          intervention="A")


  expectdataB=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                         iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8), eta=c(0.5, 0.75), tau=c(12,4),
                         rho.old=c(0.22,0.21), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.77,0.77), alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),rho2=c(1,1),
                         iota.old=c(5/7/1000, 1000/7/1000), nu.old=c(5,8), tau.old=c(12,4), eta.old=c(0.5,0.75),rho2.old=c(1,1),
                         iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(0.85,0.85),
                         intervention="B")

  expectdataBd=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                          iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8), eta=c(0.5, 0.75), tau=c(12,4),
                          rho.old=c(0.22,0.21),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                          rho.new=c(0.77,0.77),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6), sigma.old=c(1/15, 1/15), sigma.new=c(1/5, 1/5),rho2=c(1,1),
                          iota.old=c(5/7/1000, 1000/7/1000), nu.old=c(5,8), tau.old=c(12,4), eta.old=c(0.5,0.75),rho2.old=c(1,1),
                          iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(0.85,0.85),
                          intervention="B")

  expectdataBB=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                         iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8), eta=c(0.5, 0.75), tau=c(12,4),rho2=c(0.7,1),
                         rho.old=c(0.22,0.21), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                         rho.new=c(0.77,0.77), alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                         iota.old=c(5/7/1000, 1000/7/1000), nu.old=c(5,8), tau.old=c(12,4), eta.old=c(0.5,0.75),rho2.old=c(0.7,1),
                         iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(0.85,0.85),
                         intervention="B")

  expectdataBBd=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                          iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8), eta=c(0.5, 0.75), tau=c(12,4),rho2=c(0.7,1),
                          rho.old=c(0.22,0.21),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                          rho.new=c(0.77,0.77),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6), sigma.old=c(1/15, 1/15), sigma.new=c(1/5, 1/5),
                          iota.old=c(5/7/1000, 1000/7/1000), nu.old=c(5,8), tau.old=c(12,4), eta.old=c(0.5,0.75),rho2.old=c(0.7,1),
                          iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(0.85,0.85),
                          intervention="B")

  expectdataBA=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                          iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8), eta=c(0.5, 0.75), tau=c(12,4),rho2=c(0.7,1),
                          rho.old=c(0.22,0.21), alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                          rho.new=c(0.77,0.77), alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6),
                          iota.old=c(5/7/1000, 1000/7/1000), nu.old=c(5,8), tau.old=c(12,4), eta.old=c(0.5,0.75),rho2.old=c(0.7,1),
                          iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(0.7,1),
                          intervention="A")

  expectdataBAd=data.frame(id=c("regionA","regionB"), alpha=c(0.17, 0.12), beta=c(0.43,0.42), omega=c(0.18, 0.13), sigma=c(1/15, 1/15), rho=c(0.22,0.21),
                           iota=c(5/7/1000, 1000/7/1000),iota_star=c(5/7/1000, 1000/7/1000), nu=c(5,8), eta=c(0.5, 0.75), tau=c(12,4),rho2=c(0.7,1),
                           rho.old=c(0.22,0.21),alpha.old=c(0.17, 0.12), beta.old=c(0.43,0.42), omega.old=c(0.18, 0.13),
                           rho.new=c(0.77,0.77),alpha.new=c(0.7, 0.7), beta.new=c(0.8,0.8), omega.new=c(0.6, 0.6), sigma.old=c(1/15, 1/15), sigma.new=c(1/5, 1/5),
                           iota.old=c(5/7/1000, 1000/7/1000), nu.old=c(5,8), tau.old=c(12,4), eta.old=c(0.5,0.75),rho2.old=c(0.7,1),
                           iota.new=c(5/7/10000,5/7/10000), nu.new=c(5,5), tau.new=c(2,2), eta.new=c(0.5,0.5),rho2.new=c(0.7,1),
                           intervention="A")

  expect_equal(format_data_simulation(df, interv0, delay = F, rcd=T, rcd_at_baseline = T), expectdata0, label = " no intervention no delay, rcd")
  expect_equal(format_data_simulation(df, interv0, delay = T, rcd=T, rcd_at_baseline = T), expectdata0_d, label = " no intervention delay, rcd")

  expect_equal(format_data_simulation(df, intervA, delay=F, rcd=T, rcd_at_baseline = T), expectdataA, label = " intervention A no delay, rcd")
  expect_equal(format_data_simulation(df, intervA, delay=T, rcd=T, rcd_at_baseline = T), expectdataAd, label = " intervention A delay, rcd")
  expect_equal(format_data_simulation(df, intervA2, delay=F, rcd=T, rcd_at_baseline = T), expectdataA, label = " intervention A no delay, rcd, rho2=NA")

  expect_equal(format_data_simulation(df, intervB, delay=F, rcd=T, rcd_at_baseline = T), expectdataB, label = " intervention A delay, rcd, rho2=0.85")
  expect_equal(format_data_simulation(df, intervB, delay=T, rcd=T, rcd_at_baseline = T), expectdataBd, label = " intervention A delay, rcd, rho2=0.85")

  expect_equal(format_data_simulation(dfB, intervB, delay=F, rcd=T, rcd_at_baseline = T), expectdataBB, label = " intervention B delay, rcd+baseline, rho2=0.85")
  expect_equal(format_data_simulation(dfB, intervB, delay=T, rcd=T, rcd_at_baseline = T), expectdataBBd, label = " intervention B delay, rcd+baseline, rho2=0.85")

  expect_equal(format_data_simulation(dfB, intervA, delay=F, rcd=T, rcd_at_baseline = T), expectdataBA, label = " intervention A delay, rcd+baseline, rho2.old=0.7")
  expect_equal(format_data_simulation(dfB, intervA, delay=T, rcd=T, rcd_at_baseline = T), expectdataBAd, label = " intervention A delay, rcd+baseline, rho2.old=0.7")
})

