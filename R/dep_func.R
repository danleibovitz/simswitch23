
#' Generate time-dependent covariates for time *i*
#'
#' The default 'dependent' function, which takes the current dataset, pre-specified parameters and 'time' i to
#' return time-varying covariates, treatment decision and switch decision at 'time' i
#'
#' @param dat The full data-frame. Must have variables "ids", "arm", "switch", "time", "treat", "Mtime",
#' and all variables in arguments 'base_var' and 'time_var'
#' @param window An integer, representing the time to return values for
#' @param base_var A vector of names of baseline variables
#' @param time_var A vector of names of time-varying variables
#' @param covar_coef
#' @param m_haz A vector of hazards of length equal to the range of values in \code{dat$time}
#' @param num_t
#' @param tcov_n
#'
#' @return # should return a binary variable for M, and continuous for all other tvar
#' @export
#'
#' @examples
dep_func <- function(
  dat,
  window,
  base_var,
  time_var,
  covar_coef,
  m_haz,
  num_t,
  tcov_n) {

  # defend against incorrect argument types
  if(class(dat) != "data.frame") stop("'dat' must be a data.frame object")
  if(class(window) != "integer") stop("'window' must be an integer object")
  if(class(base_var) != "character") stop("'base_var' must be a character object")
  if(class(time_var) != "character") stop("'time_var' must be a character object")
  if(class(covar_coef) != "list") stop("'covar_coef' must be a list object")
  if(any(! c("baseline", "varying") %in% names(covar_coef))) stop("'covar_coef' must contain both a $baseline element and a $varying element")
  if(class(m_haz) != "integer" & class(m_haz) != "numeric") stop("'m_haz' must be an integer object or a numeric object")
  if(class(num_t) != "integer") stop("'num_t' must be an integer object")
  if(class(tcov_n) != "character") stop("'tcov_n' must be a character object")

  # TODO defend against incorrect argument ranges or absent specific values
  if(!all(c("time", "treat", "M", "Mtime", base_var, time_var) %in% names(dat))) {
    stop("'dat' data.frame must contain variables 'time', 'treat', 'M', 'Mtime', as well as all elements of
         'base_var' and 'time_var'")
  }

  # TODO defend against incorrect argument dimensions
  if(max(dat$time) != length(m_haz)) stop("The length of the 'm_haz' argument must match the 'time' variable
                                          in the 'dat' argument data.frame")
  if(length(i) != 1) stop("'window' must be length 1")
  if(length(base_var) != some_length) stop("") # TODO
  if(length(time_var) != some_length) stop("") # TODO
  if(all(dim(covar_coef$baseline) != c(dim, dim))) stop("") # TODO
  if(all(dim(covar_coef$varying) != c(dim, dim))) stop("") # TODO
  if(length(num_t) != some_length) stop("") # TODO
  if(length(tcov_names) != some_length) stop("") # TODO




  if(window == 1){# if window is 1, only a function of baseline.
    retval <- (as.matrix(dat[dat$time == window, names(dat) %in% base_var]) %*% covar_coef$baseline) + MASS::mvrnorm(n = sum(dat$time == window), rep.int(0,num_t), diag(0.1,nrow = num_t))
    retval[,1] <- 0 # Set initial M values to 0
    retval <- cbind(retval, rep.int(0, length(retval[,1]))) # set initial Mtimes
    return(retval)
  }else{# if window is not 1, a function of baseline, previous tvar and previous treat
    retval <- (as.matrix(dat[dat$time == (window-1), names(dat) %in% c("treat", time_var)]) %*% covar_coef$varying) + MASS::mvrnorm(n = sum(dat$time==(window-1)), rep.int(0, num_t), diag(0.015, nrow = num_t)) # TODO sum(dat$time==(window-1)) can probably just be n
    # change column of metastatic disease to binary, based upon hazard function
    retval[,1] <- ifelse(dat$M[dat$time == window-1] == 1, 1, 0) # if the previous window M is 1, continue M
    retval[retval[,1] == 0, 1] <- (stats::rbinom( n = length(retval[, 1]), size = 1,
                                                  prob = 1 - exp(-exp(log(m_haz[window]) +
                                                                        as.matrix(dat[dat$time == window-1, names(dat) %in% c("treat", tcov_n)]) %*% covar_coef$varying[,1])) ))[retval[,1] == 0] # randomly assign the treat variable with probability
    retval <- cbind(retval, rep.int(0, length(retval[,1])))
    retval[,num_t + 1] <- ifelse(retval[,1] == 0, 0,  # add a column representing Mtime, and calculate Mtime.
                                 ifelse(dat$M[dat$time == window-1] == 0, 1,
                                        dat$Mtime[dat$time == window-1]+1))
    return(retval)
  }
}
