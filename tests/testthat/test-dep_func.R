# dat = fulldat | data.frame, 40,000 x 12, names(fulldat) ==  [1] "ids"    "arm"    "switch" "time"   "treat"  "b1" "b2" "b3" "M" "v1" "v2" "Mtime"
# window = i, 1
# base_var=bcov_names, "b1" "b2" "b3"
# time_var=tcov_names, "M"  "v1" "v2"
# covar_coef = covar_coef, list(baseline, varying). baseline 3x3 (colsums = c(1,1,1)), varying 4x3
# m_haz = m_haz, vector 1x100, all positive
# num_t = num_tvar, 3
# tcov_n = tcov_names, "M"  "v1" "v2"

# Define test arguments
len <- 3 # number of patients
tim <- 3 # number of follow-up times
df <- data.frame(ids = rep(c(1:tim), each = len), # TODO replace with a call to fd_generator().
                 arm = c(rep(0, 6), rep(1, 3)),
                 switch = c(0,0,0,0,1,0,0,0,0),
                 time= rep(c(1:tim), len),
                 treat = c(0,0,0,0,1,1,1,1,1),
                 b1 = rep(sample(0:1, 3, replace = TRUE), each=len),
                 b2 = rep(runif(3, -1, 2), each=3),
                 b3 = rep(runif(3, -1, 2), each=3),
                 M = c(0,0,1,0,1,1,1,1,1),
                 v1 = rep(runif(tim*len, -1, 2)),
                 v2 = rep(runif(tim*len, -1, 2)),
                 Mtime = c(0,0,1,0,1,2,1,2,3))
w <- 2 # set valid window
bv <- c("b1", "b2", "b3") # set valid base_var
tv <- c("M", "v1", "v2") # set valid time_var
cc <- list( # TODO replace with a call to cc_generator
  baseline = matrix(sample(1:(num_bvar*num_tvar), num_bvar*num_tvar), ncol = num_tvar),
  varying = matrix(sample(1:(num_tvar*(num_tvar+1)), num_tvar*(num_tvar+1)), ncol = num_tvar))
cc$baseline <- LICORS::normalize(cc$baseline, byrow = FALSE)
cc$varying <- LICORS::normalize(cc$varying, byrow = FALSE)
cc$varying[1,] <- -cc$varying[1,]
cc$varying <- cc$varying/100
cc$varying[2:(num_tvar+1),1] <- cc$varying[2:(num_tvar+1),1]*150
for(j in 1:num_tvar){
  cc$varying[(j+1),j] <- 1
}
mh <- weihaz(1:3, 1, 2)
nt <- 3
tn <- c("M", "v1", "v2")

# test that incorrect argument types throw errors ####
test_that("Incorrect 'dat' argument throws error", {
  expect_error(
    dep_func(dat = "TODO", window = w, base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'window' argument throws error", {
  expect_error(
    dep_func(dat = df, window = "TODO", base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'base_var' argument throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = "TODO",
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'time_var' argument throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = bv,
             time_var = "TODO", covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'covar_coef' argument throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = bv,
             time_var = tv, covar_coef = "TODO", m_haz = mh, num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'm_haz' argument throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = "TODO", num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'num_t' argument throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = "TODO", tcov_n = tn)
  )
})
test_that("Incorrect 'tcov_n' argument throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = "TODO")
  )
})


# test incorrect argument dimensions throw errors ####
df_corrupt <- df
df_corrupt[10,] <- df_corrupt[1,]
df_corrupt$time[10] <- 4
test_that("Incorrect 'dat' argument dimensions throws error", {
  expect_error(
    dep_func(dat = df_corrupt, window = w, base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'window' argument dimensions throws error", {
  expect_error(
    dep_func(dat = df, window = c(1,2), base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'base_var' argument dimensions throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = c("b1", "b2", "b3", "b4"),
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'time_var' argument dimensions throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = bv,
             time_var = c("M", "v1", "v2", "v3"), covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'covar_coef' argument dimensions throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = bv,
             time_var = tv, covar_coef = TODO, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'm_haz' argument dimensions throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = weihaz(1:4, 1, 2), num_t = nt, tcov_n = tn)
  )
})
test_that("Incorrect 'num_t' argument dimensions throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = c(1,4), tcov_n = tn)
  )
})
test_that("Incorrect 'tcov_n' argument dimensions throws error", {
  expect_error(
    dep_func(dat = df, window = w, base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = c("M", "v1"))
  )
})

# test that fulldat containing incorrect variables throws error ####
df_corrupt <- df
names(df_corrupt)[1] <- "hello" # change one name
test_that("'fulldat' argument containing incorrect variables throws error", {
  expect_error(
    dep_func(dat = df_corrupt, window = w, base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})

# test that results of correct type and dimension are returned ####
test_that("correct arguments return vector of correct dimensions", {
  expect_error( # TODO set up expect_length() to test dim or expect_equal() to test return-value class
    dep_func(dat = df, window = w, base_var = bv,
             time_var = tv, covar_coef = cc, m_haz = mh, num_t = nt, tcov_n = tn)
  )
})

