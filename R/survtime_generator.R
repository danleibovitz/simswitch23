# Discrete-time survival time generator


#' Title
#'
#' @param x A data.frame resulting from a call to fd_generator()
#' @param hazard A hazard generating function. The default is haz_func(). The function must take arguments:
#' t, x, betas, b_haz, and ncov. For details, see the documentation of haz_func()
#' @param betas A data.frame with length equal to the length of 'x', and width equal to (ncov + 3)
#' @param ncov The number of covariates (the number of baseline covariates plus the number of time-varying covariates)
#' @param stime The number of follow-up times
#' @param idvar The name of the variable in 'x' which holds patient ids.
#' @param ids Patient ids.
#' @param b_haz A vector of "baseline hazards" for all patients, with length equal to 'stime'
#' @param n TODO
#'
#' @return
#' @export
#'
#' @examples
survtime_generator <- function(
  b_haz,
  betas,
  hazard,
  ids,
  idvar,
  n,
  ncov,
  stime,
  x) {

  # sdat <- survtime_generator( x = fulldat, hazard = haz, betas = beta.mat, ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat$ids),
  #                 b_haz = b_haz)
  # Defend against incorrect argument types
  if(!is.data.frame(x)) stop("'x' must be of class 'data.frame'")
  if(!is.function(hazard)) stop("'hazard' must be of class 'function'")
  if(!is.data.frame(betas)) stop("'betas' must be of class 'data.frame'")
  if(!is.integer(ncov)) stop("'ncov' must be of class 'integer'")
  if(!is.integer(stime)) stop("'stime' must be of class 'integer'")
  if(!is.character(idvar)) stop("'idvar' must be of class 'character'")
  if(!is.integer(ids)) stop("'ids' must be of class 'integer'")
  if(!is.numeric(b_haz)) stop("'b_haz' must be of class 'numeric'")

  # Defend against incorrect argument dimensions
  if(dim(x) != stime*n) stop("'x' must be of dim [stime*n]x[ncov+6]")
  if(dim(betas) != ncov + 3) stop("'betas' must be of dim (ncov + 3)")
  if(length(stime) != 1) stop("'stime' must be of length 1")
  if(length(idvar) != 1) stop("'idvar' must be of length 1")
  if(length(ids) != stime*n) stop("'ids' must be of length (stime*n)")
  if(length(b_haz) != stime) stop("'b_haz' must be of length equal to stime")

  # Defend against hazard() function with incorrect arguments
  # TODO how do you check for the arguments a function takes?
  if(args(hazard) != c("t", "x", "betas", "b_haz", "ncov")) stop("'hazard' must accept arguments:
                                                                 't', 'x', 'betas', 'b_haz', 'ncov'")

  store <- ids # create repository of patients with unobserved events
  df <- data.frame(ids = ids, eventtime = rep.int(0, n), status = rep.int(0, n)) # create response dataframe
  i <- 1 # set iterator

  # for all time windows: probability of failure in that window is hazard produced by haz() function
  while(any(df$eventtime == 0) & i <= stime){ # TODO should this be a while loop, and stop when all events are observed?

    # randomly assign event to remaining patients
    df$eventtime[df$ids %in% store] <- ifelse(stats::rbinom(n = length(store), size = 1, prob = hazard(
      t=i, x=x[x$ids %in% store,], betas = betas[betas$ids %in% store,], b_haz = b_haz, ncov = ncov
    )) == 1, i,0) # assign current time to for patients who observe in this window

    # update store
    store <- df$ids[df$eventtime == 0]
    i <- i + 1
  }

  # update statuses
  df$status[df$eventtime != 0] <- 1
  # apply administrative censoring to any remaining patients, and leave remaining statuses as 0
  df$eventtime[df$eventtime == 0] <- stime

  return(df)
}
