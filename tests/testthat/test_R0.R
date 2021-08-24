test_that("test R0 calculation", {

  calc_R0=get_r0_vivax(lambda=0.01, f=1/72, r=1/60, gamma=1/223)
  expect_equal(calc_R0, 1.483914 , tolerance = 1e-06, label = "case 1: R0")

  calc_R0_2=get_r0_vivax(lambda=0.012, f=1/40, r=1/65, gamma=1/338)
  expect_equal(calc_R0_2, 3.119468 , tolerance = 1e-06, label = "case 2: R0")

})

test_that("test Rc calculation", {

  calc_Rc=get_rc_vivax(lambda=0.013, f=1/72, r=1/60, gamma=1/223, alpha = 0.2, beta=0.6)
  expect_equal(calc_Rc, 1.347749 , tolerance = 1e-06, label = "case 1: Rc")

  calc_Rc_2=get_rc_vivax(lambda=0.015, f=1/40, r=1/65, gamma=1/338, alpha = 0.1, beta=0.8)
  expect_equal(calc_Rc_2, 2.856866 , tolerance = 1e-06, label = "case 2: Rc")

})



test_that("test Rc calculation. with vector control", {

  calc_Rc=get_rc_vivax(lambda=0.013, f=1/72, r=1/60, gamma=1/223, alpha = 0.2, beta=0.6, omega=1)
  expect_equal(calc_Rc, 1.347749 , tolerance = 1e-06, label = "no VC")

  calc_Rc_VC=get_rc_vivax(lambda=0.013, f=1/72, r=1/60, gamma=1/223, alpha = 0.2, beta=0.6, omega=0.7)
  expect_equal(calc_Rc_VC, 0.9434241 , tolerance = 1e-06, label = "no VC")

  calc_Rc_VC2=get_rc_vivax(lambda=0.015, f=1/40, r=1/65, gamma=1/338, alpha = 0.1, beta=0.8,omega=0.6 )
  expect_equal(calc_Rc_VC2, 1.714119 , tolerance = 1e-06, label = "case 2: Rc")

})
