test_that(
  "calc_cumn works.",
  {
    n <- list(list("T" = 10, "P" = 20, "C" = 30),
              list("T" = 40, "P" = 50, "C" = 60))
    D <- list(n = n)
    reference <- list(list("T" = 10, "P" = 20, "C" = 30),
                      list("T" = 50, "P" = 70, "C" = 90))
    expect_identical(calc_cumn(D), reference)
  })
test_that(
  "calc_cumc works.",
  {
    stagec <- list(list("T" = 1, "P" = 2, "C" = 3),
              list("T" = 4, "P" = 5, "C" = 6))
    D <- list(stagec = stagec)
    reference <- list(list("T" = 1, "P" = 2, "C" = 3),
                      list("T" = 5, "P" = 7, "C" = 9))
    expect_identical(calc_cumc(D), reference)
  })
test_that(
  "calc_n_from_c works.",
  {
    nT1 <- 10
    stagec <- list(list("T" = 1, "P" = 2, "C" = 3),
                   list("T" = 4, "P" = 5, "C" = 6))
    D <- list(stagec = stagec)
    reference <- list(list("T" = 10, "P" = 20, "C" = 30),
                      list("T" = 40, "P" = 50, "C" = 60))
    expect_identical(calc_n_from_c(nT1, D), reference)
  })
test_that(
  "calc_c works.",
  {
    n <- list(list("T" = 10, "P" = 20, "C" = 30),
              list("T" = 40, "P" = 50, "C" = 60))
    D <- list(n = n)
    reference <- list(list("T" = 1, "P" = 2, "C" = 3),
                      list("T" = 4, "P" = 5, "C" = 6))
    expect_identical(calc_c(D), reference)
  })
test_that(
  "calc_gamma works.",
  {
    var <- list("T" = 1, "P" = 2, "C" = 3)
    cumc <- list(list("T" = 1, "P" = 2, "C" = 3),
                   list("T" = 4, "P" = 5, "C" = 6))
    D <- list(var = var, cumc = cumc)
    reference <- list(list("TP" = sqrt(2), "TC" = sqrt(2)),
                      list("TP" = sqrt(1/4 + 2/5), "TC" = sqrt(1/4 + 3/6)))
    expect_identical(calc_gamma(D), reference)
  })


