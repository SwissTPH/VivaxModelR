test_that("test compare MDA in non delay model with non MDA model (no RCD)", {

  parameters_mda=list("r"=1/60, "gamma"=1/383,
                      "f"=1/69,"lambda"=0.015,"delta"=3.562799e-05,
                      "alpha"=0.18*0.95, "beta"=0.431, "rho"=0.18,"omega"=1,
                      "I0"=0.0284293, "S0"=0.5629338, "Sl"=0.2010802, "Il"=0.2075568,
                      "h"=0.0008539906, "hr"=0.000524557, "hl"=0, "hh"=0, "hhl"=0,
                      "MDAcov"=0.5, "MDAp_length"=30, "MDArad_cure"=0.5)

  parameters_mda_0c=parameters_mda; parameters_mda_0c$MDAcov=0
  parameters_mda_0p=parameters_mda; parameters_mda_0p$MDAp_length=1
  parameters_mda_0r=parameters_mda; parameters_mda_0r$MDArad_cure=0

  vivax_0=simulate_vivax_ode(parameters=parameters_mda , ODEmodel=ode_vivax_cm_importation_vc, maxtime=5000, year=FALSE)
  mda_0c=simulate_vivax_mda_ode(parameters=parameters_mda_0c, ODEmodel=ode_vivax_cm_importation_vc, ODEmodel_mda =ode_vivax_mda, maxtime=5000, year=FALSE)
  mda_0p=simulate_vivax_mda_ode(parameters=parameters_mda_0p, ODEmodel=ode_vivax_cm_importation_vc, ODEmodel_mda =ode_vivax_mda, maxtime=5000, year=FALSE)
  mda_0r=simulate_vivax_mda_ode(parameters=parameters_mda_0r, ODEmodel=ode_vivax_cm_importation_vc, ODEmodel_mda =ode_vivax_mda, maxtime=5000, year=FALSE)
  mda_cp=simulate_vivax_mda_ode(parameters=parameters_mda, ODEmodel=ode_vivax_cm_importation_vc, ODEmodel_mda =ode_vivax_mda, maxtime=5000, year=FALSE)

  row.names(vivax_0)=NULL
  row.names(mda_0c)=NULL
  row.names(mda_cp)=NULL
  row.names(mda_0p)=NULL
  row.names(mda_0r)=NULL

  expect_equal(vivax_0, mda_0c, tolerance = 1e-08, label = "MDAcov=0")
  expect_equal(mda_cp[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_equal(mda_0p[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_lt(mda_cp[30,]$I, vivax_0[30,]$I, label = "MDA better than no MDA")
  expect_lt(mda_cp[40,]$I, mda_0r[40,]$I, label = "radical cure is better")
  expect_lt(mda_0p[30,]$I, vivax_0[30,]$I, label = "small prophylaxis time better than no MDA")
  expect_lt(mda_cp[30,]$I, mda_0p[30,]$I, label = "longer prophylaxis time better than shorter")
})


test_that("test compare MDA in non delay model with non MDA model (RCD model but no RCD)", {

  parameters_mda=list("r"=1/60, "gamma"=1/383,
                      "f"=1/69,"lambda"=0.015,"delta"=3.562799e-05,
                      "alpha"=0.18*0.95, "beta"=0.431, "rho"=0.18,"omega"=1,
                      "I0"=0.0284293, "S0"=0.5629338, "Sl"=0.2010802, "Il"=0.2075568,
                      "h"=0.0008539906, "hr"=0.000524557, "hl"=0, "hh"=0, "hhl"=0,
                      "tau"=5, "nu"=5, "iota"=0, "eta"=1,"rho2"=1,
                      "MDAcov"=0.5, "MDAp_length"=30, "MDArad_cure"=0.5)

  parameters_mda_0c=parameters_mda; parameters_mda_0c$MDAcov=0
  parameters_mda_0p=parameters_mda; parameters_mda_0p$MDAp_length=1
  parameters_mda_0r=parameters_mda; parameters_mda_0r$MDArad_cure=0

  vivax_0=simulate_vivax_ode(parameters=parameters_mda , ODEmodel=ode_vivax_cm_importation_vc, maxtime=5000, year=FALSE)
  mda_0c=simulate_vivax_mda_ode(parameters=parameters_mda_0c, ODEmodel=ode_vivax_rcd, ODEmodel_mda = ode_vivax_rcd_mda, maxtime=5000, year=FALSE)
  mda_0p=simulate_vivax_mda_ode(parameters=parameters_mda_0p, ODEmodel=ode_vivax_rcd, ODEmodel_mda = ode_vivax_rcd_mda, maxtime=5000, year=FALSE)
  mda_0r=simulate_vivax_mda_ode(parameters=parameters_mda_0r, ODEmodel=ode_vivax_rcd, ODEmodel_mda =ode_vivax_mda, maxtime=5000, year=FALSE)
  mda_cp=simulate_vivax_mda_ode(parameters=parameters_mda, ODEmodel=ode_vivax_rcd, ODEmodel_mda = ode_vivax_rcd_mda, maxtime=5000, year=FALSE)

  row.names(vivax_0)=NULL
  row.names(mda_0c)=NULL
  row.names(mda_cp)=NULL
  row.names(mda_0p)=NULL
  row.names(mda_0r)=NULL

  expect_equal(vivax_0, mda_0c, tolerance = 5e-09, label = "MDAcov=0")
  expect_equal(mda_cp[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_equal(mda_0p[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_lt(mda_cp[30,]$I, vivax_0[30,]$I, label = "MDA better than no MDA")
  expect_lt(mda_cp[40,]$I, mda_0r[40,]$I, label = "radical cure is better")
  expect_lt(mda_0p[30,]$I, vivax_0[30,]$I, label = "small prophylaxis time better than no MDA")
  expect_lt(mda_cp[30,]$I, mda_0p[30,]$I, label = "longer prophylaxis time better than shorter")
})


test_that("test compare MDA in non delay model with non MDA model (with RCD)", {

  parameters_mda=list("r"=1/60, "gamma"=1/383,
                      "f"=1/69,"lambda"=0.018,"delta"=3.562799e-05,
                      "alpha"=0.18*0.95, "beta"=0.431, "rho"=0.18,"omega"=1,
                      "I0"=0.01206, "S0"=0.8225808, "Sl"=0.08523101, "Il"=0.08012824,
                      "h"=0.0003336124, "hr"=0.0002223417, "hl"=0, "hh"=0, "hhl"=0,
                      "tau"=5, "nu"=5, "iota"=50/7/10000, "eta"=1,"rho2"=1,
                      "MDAcov"=0.5, "MDAp_length"=30, "MDArad_cure"=0.5)

  parameters_mda_0c=parameters_mda; parameters_mda_0c$MDAcov=0
  parameters_mda_0p=parameters_mda; parameters_mda_0p$MDAp_length=1
  parameters_mda_0r=parameters_mda; parameters_mda_0r$MDArad_cure=0

  vivax_0=simulate_vivax_ode(parameters=parameters_mda , ODEmodel=ode_vivax_rcd, maxtime=5000, year=FALSE)
  mda_0c=simulate_vivax_mda_ode(parameters=parameters_mda_0c, ODEmodel=ode_vivax_rcd, ODEmodel_mda = ode_vivax_rcd_mda, maxtime=5000, year=FALSE)
  mda_0p=simulate_vivax_mda_ode(parameters=parameters_mda_0p, ODEmodel=ode_vivax_rcd, ODEmodel_mda = ode_vivax_rcd_mda, maxtime=5000, year=FALSE)
  mda_0r=simulate_vivax_mda_ode(parameters=parameters_mda_0r, ODEmodel=ode_vivax_cm_importation_vc, ODEmodel_mda =ode_vivax_mda, maxtime=5000, year=FALSE)
  mda_cp=simulate_vivax_mda_ode(parameters=parameters_mda, ODEmodel=ode_vivax_rcd, ODEmodel_mda = ode_vivax_rcd_mda, maxtime=5000, year=FALSE)

  row.names(vivax_0)=NULL
  row.names(mda_0c)=NULL
  row.names(mda_cp)=NULL
  row.names(mda_0p)=NULL
  row.names(mda_0r)=NULL

  expect_equal(vivax_0, mda_0c, tolerance = 9e-07, label = "MDAcov=0")
  expect_equal(mda_cp[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_equal(mda_0p[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_lt(mda_cp[30,]$I, vivax_0[30,]$I, label = "MDA better than no MDA")
  expect_lt(mda_cp[40,]$I, mda_0r[40,]$I, label = "radical cure is better")
  expect_lt(mda_0p[30,]$I, vivax_0[30,]$I, label = "small prophylaxis time better than no MDA")
  expect_lt(mda_cp[30,]$I, mda_0p[30,]$I, label = "longer prophylaxis time better than shorter")
})



test_that("test compare MDA in delay model with non MDA model (no RCD)", {

  parameters_mda=list("r"=1/60, "gamma"=1/383,
                      "f"=1/69,"lambda"=0.013,"delta"=3.562799e-05,
                      "alpha"=0.18*0.95, "beta"=0.431, "rho"=0.18,"omega"=1, "sigma"=1/15,
                      "U0"=0.03029989, "S0"=0.5172503, "Sl"=0.2172999, "Ul"=0.2246329, "Tl"=0.01020732, "T0"=0.0003098104,
                      "h"=0.0009225553, "hl"=0, "hh"=0, "hhl"=0,
                      "MDAcov"=0.5, "MDAp_length"=30, "MDArad_cure"=0.5)

  parameters_mda_0c=parameters_mda; parameters_mda_0c$MDAcov=0
  parameters_mda_0p=parameters_mda; parameters_mda_0p$MDAp_length=1
  parameters_mda_0r=parameters_mda; parameters_mda_0r$MDArad_cure=0

  vivax_0=simulate_vivax_delay_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay, maxtime=5000, year=FALSE)
  mda_0c=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0c , ODEmodel=ode_vivax_delay, ODEmodel_mda =ode_vivax_delay_mda, maxtime=5000, year=FALSE)
  mda_0p=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0p , ODEmodel=ode_vivax_delay, ODEmodel_mda =ode_vivax_delay_mda, maxtime=5000, year=FALSE)
  mda_cp=simulate_vivax_delay_mda_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay, ODEmodel_mda =ode_vivax_delay_mda, maxtime=5000, year=FALSE)
  mda_0r=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0r , ODEmodel=ode_vivax_delay, ODEmodel_mda =ode_vivax_delay_mda, maxtime=5000, year=FALSE)

  row.names(vivax_0)=NULL
  row.names(mda_0c)=NULL
  row.names(mda_cp)=NULL
  row.names(mda_0p)=NULL
  row.names(mda_0r)=NULL

  expect_equal(vivax_0, mda_0c, tolerance = 1e-08, label = "MDAcov=0")
  expect_equal(mda_cp[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_equal(mda_0p[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_lt(mda_cp[30,]$I, vivax_0[30,]$I, label = "MDA better than no MDA")
  expect_lt(mda_cp[40,]$I, mda_0r[40,]$I, label = "radical cure is better")
  expect_lt(mda_0p[30,]$I, vivax_0[30,]$I, label = "small prophylaxis time better than no MDA")
  expect_lt(mda_cp[30,]$I, mda_0p[30,]$I, label = "longer prophylaxis time better than shorter")
})


test_that("test compare MDA in delay model with non MDA model (RCD model but no RCD)", {

  parameters_mda=list("r"=1/60, "gamma"=1/383,
                      "f"=1/69,"lambda"=0.015,"delta"=3.562799e-05,
                      "alpha"=0.18*0.95, "beta"=0.431, "rho"=0.18,"omega"=1, "sigma"=1/15,
                      "U0"=0.03029989, "S0"=0.5172503, "Sl"=0.2172999, "Ul"=0.2246329, "Tl"=0.01020732, "T0"=0.0003098104,
                      "h"=0.0009225553,"hl"=0, "hh"=0, "hhl"=0,
                      "tau"=5, "nu"=5, "iota"=0, "eta"=1,"rho2"=1,
                      "MDAcov"=0.5, "MDAp_length"=30, "MDArad_cure"=0.5)

  parameters_mda_0c=parameters_mda; parameters_mda_0c$MDAcov=0
  parameters_mda_0p=parameters_mda; parameters_mda_0p$MDAp_length=1
  parameters_mda_0r=parameters_mda; parameters_mda_0r$MDArad_cure=0

  vivax_0=simulate_vivax_delay_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  mda_0c=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0c , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=5000, year=FALSE)
  mda_0p=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0p , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=5000, year=FALSE)
  mda_cp=simulate_vivax_delay_mda_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=5000, year=FALSE)
  mda_0r=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0r , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=5000, year=FALSE)

  mda_cp_r=simulate_vivax_delay_mda_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay, ODEmodel_mda=ode_vivax_delay_rcd_referral_mda, maxtime=5000, year=FALSE)

  row.names(vivax_0)=NULL
  row.names(mda_0c)=NULL
  row.names(mda_cp)=NULL
  row.names(mda_0p)=NULL
  row.names(mda_0r)=NULL
  row.names(mda_cp_r)=NULL

  expect_equal(vivax_0, mda_0c, tolerance = 5e-09, label = "MDAcov=0")
  expect_equal(mda_cp[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_equal(mda_0p[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_lt(mda_cp[30,]$I, vivax_0[30,]$I, label = "MDA better than no MDA")
  expect_lt(mda_cp[40,]$I, mda_0r[40,]$I, label = "radical cure is better")
  expect_lt(mda_0p[30,]$I, vivax_0[30,]$I, label = "small prophylaxis time better than no MDA")
  expect_lt(mda_cp[30,]$I, mda_0p[30,]$I, label = "longer prophylaxis time better than shorter")

  expect_equal(mda_cp_r, mda_cp, tolerance = 1e-09, label = "referral = non referral (because iota=0)")
})

test_that("test compare MDA in delay model with non MDA model (with RCD)", {

  parameters_mda=list("r"=1/60, "gamma"=1/383,
                      "f"=1/69,"lambda"=0.015,"delta"=3.562799e-05, "sigma"=1/15,
                      "alpha"=0.18*0.95, "beta"=0.431, "rho"=0.18,"omega"=1,
                      "U0"=0.03029989, "S0"=0.5172503, "Sl"=0.2172999, "Ul"=0.2246329, "Tl"=0.01020732, "T0"=0.0003098104,
                      "h"=0.0009225553,"hl"=0, "hh"=0, "hhl"=0,
                      "tau"=5, "nu"=5, "iota"=50/7/10000, "eta"=1,"rho2"=1,
                      "MDAcov"=0.5, "MDAp_length"=30, "MDArad_cure"=0.5)

  parameters_mda_0c=parameters_mda; parameters_mda_0c$MDAcov=0
  parameters_mda_0p=parameters_mda; parameters_mda_0p$MDAp_length=1
  parameters_mda_0r=parameters_mda; parameters_mda_0r$MDArad_cure=0


  vivax_0=simulate_vivax_delay_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay_rcd_no_referral, maxtime=5000, year=FALSE)
  mda_0c=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0c , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=5000, year=FALSE)
  mda_0p=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0p , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=5000, year=FALSE)
  mda_cp=simulate_vivax_delay_mda_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=5000, year=FALSE)
  mda_0r=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0r , ODEmodel=ode_vivax_delay_rcd_no_referral, ODEmodel_mda =ode_vivax_delay_rcd_no_referral_mda, maxtime=5000, year=FALSE)

  vivax_0_r=simulate_vivax_delay_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay_rcd_referral, maxtime=5000, year=FALSE)
  mda_cp_r=simulate_vivax_delay_mda_ode(parameters=parameters_mda , ODEmodel=ode_vivax_delay_rcd_referral, ODEmodel_mda=ode_vivax_delay_rcd_referral_mda, maxtime=5000, year=FALSE)
  mda_0p_r=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0p , ODEmodel=ode_vivax_delay_rcd_referral, ODEmodel_mda=ode_vivax_delay_rcd_referral_mda, maxtime=5000, year=FALSE)
  mda_0c_r=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0c , ODEmodel=ode_vivax_delay_rcd_referral, ODEmodel_mda=ode_vivax_delay_rcd_referral_mda, maxtime=5000, year=FALSE)
  mda_0r_r=simulate_vivax_delay_mda_ode(parameters=parameters_mda_0r , ODEmodel=ode_vivax_delay_rcd_referral, ODEmodel_mda=ode_vivax_delay_rcd_referral_mda, maxtime=5000, year=FALSE)

  row.names(vivax_0)=NULL
  row.names(mda_0c)=NULL
  row.names(mda_cp)=NULL
  row.names(mda_0p)=NULL
  row.names(mda_0r)=NULL

  row.names(vivax_0_r)=NULL
  row.names(mda_0c_r)=NULL
  row.names(mda_cp_r)=NULL
  row.names(mda_0p_r)=NULL
  row.names(mda_0r_r)=NULL

  expect_equal(vivax_0, mda_0c, tolerance = 1e-07, label = "MDAcov=0")
  expect_equal(mda_cp[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_equal(mda_0p[5000,], vivax_0[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_lt(mda_cp[30,]$I, vivax_0[30,]$I, label = "MDA better than no MDA")

  expect_equal(vivax_0_r, mda_0c_r, tolerance = 1e-07, label = "MDAcov=0")
  expect_equal(mda_cp_r[5000,], vivax_0_r[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_equal(mda_0p_r[5000,], vivax_0_r[5000,], tolerance = 1e-06, label = "long term, MDA same as no MDA")
  expect_lt(mda_cp_r[30,]$I, vivax_0_r[30,]$I, label = "MDA better than no MDA")

  expect_gt(mda_cp_r[30,]$I, mda_cp[30,]$I, label = "MDA better than no MDA")
  expect_gt(mda_cp_r[2000,]$I, mda_cp[2000,]$I, label = "MDA better than no MDA")
})

