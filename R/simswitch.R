
#' Simulate TS data and apply treatment-effect-estimating methods
#'
#' @param add_tvar Number of covariates to add which have no effect on overall survival, but which are nonetheless included in IPCW, TSE and RPSFTM models (non-negative integer).
#' @param b_allowance
#' @param b_haz Basic hazard function determining survival. In a simulation with no covariates and 0 treatment effect, this is the universal hazard function.
#' Default is the hazard associated with a weibull distribution with \emph{shape} = 1 and \emph{scale} = stime. A user-defined entry must be vector or array of length \emph{stime} representing the hazard up to and inculding each time point. For specifying a parametric basic hazard function, see arguments \emph{b_scale} and \emph{b_shape}.
#' @param b_mag The intial magnitude by which the hazard function \emph{b_haz} is adjusted if the desired number of events is not observed. Default is 0.5 (and the default value is usually sufficient). Must be postive and numeric.
#' @param b_scale
#' @param b_shape
#' @param bcov
#' @param beta.mat
#' @param cens_flag
#' @param covar_coef This is a new comment
#' @param dep_func
#' @param haz
#' @param hide_tvar add option for hiding relevent covariates from switching modeling
#' @param ipcw_exclude
#' @param m_allowance
#' @param m_inflation
#' @param m_fidelity
#' @param m_hard
#' @param m_haz
#' @param m_mag
#' @param m_scale
#' @param m_shape
#' @param n   Number of participants in simulated trials.
#' @param num_bvar Number of basline covariates to include.
#' @param num_tvar Number of time-varying covariates to include.
#' @param prop_cens Proportion of patients
#' @param prop_cens_allowance
#' @param prop_cont_event
#' @param prop_switch set proportion of control participants to switch
#' @param prop_trt
#' @param prop_trt_event   set proportion of participants on experimental treatment
#' @param recens
#' @param rerun_lim   set while loop limit
#' @param s_allowance
#' @param s_haz
#' @param s_mag
#' @param s_scale
#' @param s_shape
#' @param seed
#' @param stime
#' @param switch_coef
#' @param t_allowance
#' @param t_mag
#' @param treat_beta
#' @param tse_dist
#' @param unfix
#' @param verbose
#' @param violate
#'
#' @return
#' @export
#'
#' @examples
simswitch <- function(add_tvar = 0,
                      b_allowance = 0.1,
                      b_haz,
                      b_mag = 0.5,
                      b_scale,
                      b_shape,
                      bcov,
                      beta.mat,
                      censoring,
                      covar_coef,
                      dep = dep_func,
                      haz = haz_func,
                      hide_tvar = 0,
                      analyses = c("ITT", "PP", "IPCW", "RPSFTM", "TSE"),
                      m_allowance = 0.1,
                      m_inflation = 2,
                      m_fidelity = 0.2, # the proportion of stime away from first M that switch can occur
                      m_hard = TRUE, # represents whether switch can happen only after M, or if we don't care
                      m_haz,
                      m_mag = 2,
                      m_scale,
                      m_shape,
                      n = 400,
                      num_bvar = 3,
                      num_cores = 0,
                      num_tvar = 3,
                      para = FALSE,
                      prop_cens,
                      prop_cens_allowance = 0.1,
                      prop_cont_event,
                      prop_switch = 0.5,
                      prop_trt = 0.5,
                      prop_trt_event,
                      recens = TRUE, # set flag for recensoring of RPSFTM and Two-Stage Model
                      nsim = 1000,
                      rerun_lim = 200,
                      s_allowance = 0.1, # how far from the proportion of switching requested is acceptable?
                      s_haz = weihaz(1:stime, s_shape, s_scale),  # set baseline switch hazard function. lam is set as lambda of exp distribution with expected value of prop_switch by stime
                      s_mag = 0.5,
                      s_scale = 0.7*stime,
                      s_shape = 2,
                      seed,
                      stime = 100,
                      switch_coef = c(runif(1, 1.5, 2), runif(num_bvar-1, 0.1, 0.3), runif(num_tvar, 0.2, 0.5)), # baseline covariate attached to "predisposition" baseline is first, and slightly higher than other baselines. default of switch_coef log hazard ratios. baseline switch coefs are smaller.
                      t_allowance = 0.1,
                      t_mag = 0.5,
                      treat_beta = -1,
                      tse_dist = c("loglogistic", "weibull", "lognormal"), # alternatives are weibull, lognormal, etc.
                      unfix = as.character(c()),
                      verbose = 2,
                      violate = "None"){



  # Set and check default parameters ####
  checkparams()

  # Set all cox-like functions and associated parameters ####
  surv_params <- get_surv_params(...)
  checksurvparams()

  # Build dataframe skeleton ####
  skeleton <- get_skeleton_data_frame(...)

  # set empty estimate vectors
  estimates <- get_empty_estimate_frame()



  # Tune hazard functions to obtain desired rates ####
  for(r in 1:nsim){
    while(desired rates not obtained){ # rates are: 1. secondary baseline, 2. switch, 3. Event,

      # get data with switching
      fulldat <- get_switch_data(...)
      # get data without switching
      fulldat_blocked <- get_blocked_data(...)

      update_hazards(...)
    }

    # Get Models and Plots: ####

    ## Unbiased model ##
    estimates$unbiased[r] <- get_unbiased()
    estimates$unbiased_plot[r] <- get_unbiased_plot()

    for(j in analyses){
      estimates$j[r] <- get_j_method()
      estimates$j_plot[r] <- get_j_plot()
    }
  }



  # build simulation plots ####
  estimates$compar_plot <- get_compar_plot()
  estimates$bias_plot <- get_bias_plot()



  return(estimates)
}


