
# Define test arguments
np <- 6
nbc <- 2
ntc <- 2
st <- 10
pt <- 0.5
v <- violate # TODO
idv <- "ids"
is <- unique(x$ids)
bh <- weihaz(scale = 0.5, shape = 0.5, x = 1:st)
ex <- fd_generator(bcov = , id_trt = , n = np, num_bvar = nbc, num_tvar = ntc, prop_trt = pt, stime = st)
h <- haz_func
bs <- bm_generator(stime = st, num_bvar = nbc, num_tvar = ntc, bcov_names = ,
                   tcov_names = , treat_beta = , violate = , n = np)

# test that incorrect argument types throw errors ####
test_that("Incorrect 'x' argument throws error", {
  expect_error(
    survtime_generator(x = "ex", hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'hazard' argument throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = "h", betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'betas' argument throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = "bs", ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'ncov' argument throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs, ncov = "nc", stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'stime' argument throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs, ncov = nc, stime = "st", idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'idvar' argument throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = "idv", ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'ids' argument throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = "is", b_haz = bh)
  )
})
test_that("Incorrect 'b_haz' argument throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = "bh")
  )
})

# test incorrect argument dimensions throw errors ####
test_that("Incorrect 'ex' argument dimension throws error", {
  expect_error(
    survtime_generator(x = ex[1:3,], hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'betas' argument dimension throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs[1:3,], ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'ncov' argument dimension throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs, ncov = c(nc, nc), stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'stime' argument dimension throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs, ncov = nc, stime = c(st,st), idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'idvar' argument dimension throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = c(idv, idv), ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'ids' argument dimension throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is[1:3], b_haz = bh)
  )
})
test_that("Incorrect 'b_haz' argument dimension throws error", {
  expect_error(
    survtime_generator(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh[1:3])
  )
})

# test that results of correct type and dimension are returned ####
test_that("Correc function call returns correct class", {
  expect_equal(
    class(survtime_generator(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)),
    "data.frame"
  )
})
