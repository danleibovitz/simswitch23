

# main driver function documentation ####

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
                      cens_flag = "Random",
                      covar_coef,
                      dep = dep_func,
                      haz = haz_func,
                      hide_tvar = 0,
                      ipcw_exclude = TRUE,   # set ipcw_exclude. If TRUE, we skip bootstrapping for IPCW. Bootstrapped IPCW may not be symmetrical, and may be more accurate, but takes hella long
                      m_allowance = 0.1,
                      m_inflation = 2,
                      m_fidelity = 0.2, # the proportion of stime away from first M that switch can occur
                      m_hard = TRUE, # represents weather switch can happen only after M, or if we don't care
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
                      reps = 1000,
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
                      tse_dist = "loglogistic", # alternatives are weibull, lognormal, etc.
                      unfix = as.character(c()),
                      verbose = 2,
                      violate = "None"){



  # Set and check default parameters ####
  params <- get_default_params()

  # set giant while loop flag
  rerun <- TRUE

  # Build dataframe ####
  fulldat <- fd_generator(...)

  # Set all cox-like functions and associated parameters ####
  surv_params <- get_surv_params()


  if(missing(beta.mat)){
    beta.mat <- get_beta_mat()
  }



  # set empty estimate vectors
  unbiased <- c()
  itt <- c()
  pp <- c()
  ipcw_est <- c()
  rpsftm_est <- c()
  tse <- c()
  secondary_baseline_observed <- c()
  switch_observed <- c()
  switch_iter <- 0 # how many times have we searched on the inner loop?
  Minb0 <- c()
  Minb1 <- c()


  # Tune hazard functions to obtain desired rates ####

  # giant while loop should begin here. At this point, we have all parameters set, and a fulldat dataframe with no time
  # varying covariates, no switching and no secondary baseline
  # TODO if para == TRUE, parallelize here. Can these be done all in parallel (including the 1st rep, where conditions are particular) or is it sequential?
  # TODO ... the fuck is going on in for(r in 1)?
  for(r in 1){

    if(r != 1) unfix <- c("B", "M", "S", "T") # if we have generated the first dataset, unfix all hazards

    while(rerun){ # iteratively update switching hazard function until we get the right proportion. The first term in this check is the number of pts whose final treatment indicator

      for(i in 1:stime){
        # set covars
        fulldat[fulldat$time == i, names(fulldat) %in% c(tcov_names, "Mtime")] <-
          dep(dat = fulldat, window = i, base_var=bcov_names, time_var=tcov_names, covar_coef = covar_coef, m_haz = m_haz, num_t = num_tvar, tcov_n = tcov_names)

        # set treatment indicator
        if(i != 1){ # if its not the first time window
          fulldat$treat[fulldat$time == i] <- ifelse(fulldat$treat[fulldat$time == i-1] == 1, 1, 0) # if the previous window treat is 1, continue treatment
        }
        # if treatment has not yet begun, probability of begining in the next window is a hazard function
        if(m_hard){
          fulldat$treat[fulldat$time == i & fulldat$treat == 0 & fulldat$Mtime > 0 & fulldat$Mtime <= ceiling(m_fidelity*stime)] <-
            stats::rbinom( n = length(fulldat$treat[fulldat$time == i & fulldat$treat == 0 & fulldat$Mtime > 0 & fulldat$Mtime <= ceiling(m_fidelity*stime)]), size = 1,
                           prob = 1 - exp(-exp(log(s_haz[i]) +
                                                 as.matrix(fulldat[fulldat$time == i & fulldat$treat == 0, names(fulldat) %in% c(bcov_names, tcov_names)]) %*% switch_coef))) # randomly assign the treat variable with probability
        }else{
          fulldat$treat[fulldat$time == i & fulldat$treat == 0] <- stats::rbinom( n = length(fulldat$treat[fulldat$time == i & fulldat$treat == 0]), size = 1,
                                                                                  prob = 1 - exp(-exp(log(s_haz[i]) +
                                                                                                        as.matrix(fulldat[fulldat$time == i & fulldat$treat == 0, names(fulldat) %in% c(bcov_names, tcov_names)]) %*% switch_coef))) # randomly assign the treat variable with probability
        }

      }

      # Generate fulldat_cont (the control dataset) by blocking switching
      fulldat_cont$treat <- fulldat_cont$arm
      for(i in 1:stime){
        # set covars
        fulldat_cont[fulldat_cont$time == i, names(fulldat_cont) %in% c(tcov_names, "Mtime")] <-
          dep(dat = fulldat_cont, window = i, base_var=bcov_names, time_var=tcov_names, covar_coef = covar_coef, m_haz = m_haz, num_t = num_tvar, tcov_n = tcov_names)
      }


      if(verbose > 1){
        print("Generating confounded survival times...")
      }

      sdat <- survtime_generator( x = fulldat, hazard = haz, betas = beta.mat, ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat$ids),
                                  b_haz = b_haz)

      if(verbose > 1){
        print("Generating un-confounded survival times...")
      }

      sdat_cont <- survtime_generator( x = fulldat_cont[fulldat_cont$arm == 0,], hazard = haz, betas = beta.mat[fulldat_cont$arm == 0,], ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat_cont$ids[fulldat_cont$arm==0]),
                                       b_haz = b_haz)


      if(verbose > 1){
        print("Adding censoring...")
      }
      # merge datasets:
      if(switch_iter == 0){ # if its the first iteration
        fulldat <- merge(fulldat, sdat, by = "ids")
        fulldat_cont <- merge(fulldat_cont, sdat_cont, by = "ids", all = TRUE)
        # replace experimental group in unconfounded dataset with experimental group in confounded dat. Theoretically, these are generated identically.
        fulldat_cont[fulldat_cont$arm == 1, ] <- fulldat[fulldat$arm == 1, ]
      }else{
        fulldat <- fulldat[, !names(fulldat) %in% c("eventtime", "status")] # first remove old sdat
        fulldat <- merge(fulldat, sdat, by = "ids")

        fulldat_cont <- fulldat_cont[, !names(fulldat_cont) %in% c("eventtime", "status")] # first remove old sdat
        fulldat_cont <- merge(fulldat_cont, sdat_cont, by = "ids", all = TRUE)
        # replace experimental group in unconfounded dataset with experimental group in confounded dat. Theoretically, these are generated identically.
        fulldat_cont[fulldat_cont$arm == 1, ] <- fulldat[fulldat$arm == 1, names(fulldat) %in% names(fulldat_cont)]
      }


      # censor at prespecified, non-administrative censoring time. Censoring will have been generated independently of covars, or dependent on covars
      if(cens_flag == "Random"){
        # TODO define a censoring function, e.g., rc_generator()
        fulldat$cens <- stime # set default censoring time
        fulldat_cont$cens <- stime # set default censoring time

        cens_ids <- sample(unique(fulldat$ids), ceiling(prop_cens*n)) # sample prop_cens of ids for censoring
        rand_cens <- rbeta(n = length(cens_ids), shape1 = 2, shape2 = 1.5) # censoring times are drawn from a beta dist
        # rand_cens <- ceiling(rand_cens)
        rand_cens <- rep(rand_cens, each=stime) # spread to length of dataset
        fulldat$cens[fulldat$ids %in% cens_ids] <- ceiling(rand_cens*fulldat$eventtime[fulldat$ids %in% cens_ids])
        fulldat_cont$cens[fulldat_cont$ids %in% cens_ids] <- ceiling(rand_cens*fulldat_cont$eventtime[fulldat_cont$ids %in% cens_ids])

        # censor observations. Take minimum of eventtime and cens, and change status to 0 if eventtime == cens
        fulldat$eventtime <- pmin(fulldat$eventtime, fulldat$cens)
        fulldat$status <- ifelse(fulldat$eventtime == fulldat$cens, 0, fulldat$status)
        fulldat_cont$eventtime <- pmin(fulldat_cont$eventtime, fulldat_cont$cens)
        fulldat_cont$status <- ifelse(fulldat_cont$eventtime == fulldat_cont$cens, 0, fulldat_cont$status)


      } else if(cens_flag == "Nonrandom"){
        # TODO nonrandom censoring function!
      }


      # Does the patient (control and experimental) observe secondary baseline (M)?
      # TODO i think pts who observe NO M at all are being given NA here
      fulldat$secbase_observed <- as.numeric(rep(sapply(unique(fulldat$ids), function(x)
        ifelse(sum(fulldat$M[fulldat$ids == x]) > 0 && fulldat$time[fulldat$Mtime == 1 & fulldat$ids == x] <= fulldat$eventtime[fulldat$Mtime == 1 & fulldat$ids == x], # id has at least some m, and time at first m is less than event time
               1,
               0)), each = stime))

      # Does the patient (only control) observe switch _before eventtime_?
      fulldat$switch_status <- rep(sapply(unique(fulldat$ids), function(x)
        ifelse(sum(fulldat$treat[fulldat$ids == x & fulldat$time <= fulldat$eventtime]) > 0 & fulldat$arm[fulldat$ids == x][1] == 0,
               1,
               0)), each = stime)

      if(r == 1){ # if we're in the first loop,
        switch_iter <- switch_iter + 1 # iterate search indicator
      }
      if(switch_iter > rerun_lim) stop("Your covariate model failed to converge. Try different hazards and/or coefficients")

      rerun <- FALSE # pre-empt a rerun, unless the following conditions are met

      # if not enough M occurences, adjust M hazard
      # TODO for now, M occurence is only being adjusted upward. do we want to give it a within-window type adjustment?
      current_m_prop <- length(fulldat$ids[fulldat$arm == 0 & fulldat$Mtime == 1 & fulldat$time <= fulldat$eventtime])
      if(! "M" %in% unfix){
        if(current_m_prop < min(n*prop_trt, prop_switch*prop_trt*n*m_inflation)){ # current_m_prop is proportion of control arm that observes M. if it's less than the min of (the whole control arm, the switchers*m_inflation ), increase m_haz and rerun
          m_haz <- m_haz*m_mag # change scale of m related haz, and rerun
          rerun <- TRUE
        }
      }


      # adjust switching proportion
      current_s_prop <- sum(fulldat$switch_status[fulldat$time == 1 & fulldat$arm == 0])/length(id_con) # get current switch proportion
      if(! "S" %in% unfix){
        if(abs(current_s_prop - prop_switch) > s_allowance){ # If we're outside the window
          prev_s_direc <- s_direc # save most recent s_direc value
          s_direc <- ifelse(current_s_prop > prop_switch, -1, 1) # identify the s_direc, i.e., in what direction the baseline hazard must be adjusted
          if((prev_s_direc == -1 & s_direc == 1) | (prev_s_direc == 1 & s_direc == -1)){ # if we overshot the allowance window, and must backtrack
            s_mag <- s_mag/2 # halve switch magnitude if weve crossed over the window
          }
          s_haz <- s_haz + s_direc*s_mag*s_haz # modify the switch hazard by a factor of s_mag in the direction of s_direc
          rerun <- TRUE
        }
      }


      # adjust control group event occurence
      current_b_prop <- sum(fulldat$status[fulldat$time == 1 & fulldat$arm == 0])/length(id_con) # get current switch proportion
      if(! "B" %in% unfix){
        if(abs(current_b_prop - prop_cont_event) > b_allowance){ # If we're outside the window
          # if(b_haz == "exponential"){
          #   prev_b_direc <- b_direc # save most recent b_direc value
          #   b_direc <- ifelse(current_b_prop > prop_cont_event, -1, 1) # identify the s_direc, i.e., in what direction the baseline hazard must be adjusted
          #   if((prev_b_direc == -1 & b_direc == 1) | (prev_b_direc == 1 & b_direc == -1)){ # if we overshot the allowance window, and must backtrack
          #     b_mag <- b_mag/2 # halve switch magnitude if weve crossed over the window
          #   }
          #   lambdas <- lambdas - b_direc*b_mag*lambdas # modify the switch hazard by a factor of s_mag in the direction of s_direc
          #   if(lambdas <= 0) stop("The scale parameter for the simsurv distribution is 0 or negative; it must be positive.")
          #   rerun <- TRUE
          # }else{
          prev_b_direc <- b_direc # save most recent b_direc value
          b_direc <- ifelse(current_b_prop > prop_cont_event, -1, 1) # identify the s_direc, i.e., in what direction the baseline hazard must be adjusted
          if((prev_b_direc == -1 & b_direc == 1) | (prev_b_direc == 1 & b_direc == -1)){ # if we overshot the allowance window, and must backtrack
            b_mag <- b_mag/2 # halve switch magnitude if weve crossed over the window
          }
          b_haz <- b_haz + b_direc*b_mag*b_haz # modify the switch hazard by a factor of s_mag in the direction of s_direc
          rerun <- TRUE
          # }
        }
      }

      # adjust trt group event occurence
      current_t_prop <- sum(fulldat$status[fulldat$time == 1 & fulldat$arm == 1])/length(id_trt) # get current switch proportion
      if(! "T" %in% unfix){
        if(abs(current_t_prop - prop_trt_event) > t_allowance){ # If we're outside the window
          prev_t_direc <- t_direc # save most recent b_direc value
          t_direc <- ifelse(current_t_prop > prop_trt_event, 1, -1) # identify the t_direc, i.e., in what direction the treatment effect must be adjusted
          if((prev_t_direc == -1 & t_direc == 1) | (prev_t_direc == 1 & t_direc == -1)){ # if we overshot the allowance window, and must backtrack
            t_mag <- t_mag/2 # halve switch magnitude if weve crossed over the window
          }
          beta.mat$treat <- beta.mat$treat + t_direc*t_mag*beta.mat$treat # modify the switch hazard by a factor of s_mag in the direction of s_direc
          rerun <- TRUE
        }
      }

    }

    rerun <- TRUE # reset inner loop rerun flag

    # if(verbose > 1){
    #   print("Performing naive method estimates...")
    # }


    # Run control model (i.e., get do-treat, unconfounded estimates)
    unbiased[r] <- exp(coef(survival::coxph(survival::Surv(eventtime, status) ~ arm, data = fulldat_cont[fulldat_cont$time == 1,])))
    # Old, single rep model: survival::coxph(survival::Surv(eventtime, status) ~ arm, data = fulldat_cont[fulldat_cont$time == 1,])

    unbiased_plot <- survminer::ggsurvplot(
      fit = survminer::surv_fit(survival::Surv(eventtime, status) ~ arm, data = fulldat_cont[fulldat_cont$time == 1,]),
      xlab = "Time",
      ylab = "OS",
      title = "KM Plots for Unconfounded Analysis", conf.int = TRUE) # %++%
    # geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)


    # Run Naive models:

    ## ITT ##
    # HR estimate:
    itt[r] <- exp(coef(survival::coxph(survival::Surv(eventtime, status) ~ arm, data = fulldat[fulldat$time == 1,])))
    # Old, single rep model: survival::coxph(survival::Surv(eventtime, status) ~ arm, data = fulldat[fulldat$time == 1,])

    itt_plot <- survminer::ggsurvplot(
      fit = survminer::surv_fit(survival::Surv(eventtime, status) ~ arm, data = fulldat[fulldat$time == 1,]),
      xlab = "Time",
      ylab = "OS",
      title = "KM Plots for ITT", conf.int = TRUE) # %++%
    # geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)


    ## Per Protocol ##
    # Per Protocol (censor at switch with no correction)
    PPtime <- sapply(unique(fulldat$ids), function(x) stime - sum(fulldat$treat[fulldat$ids == x])) # get total time off treatment. since there is no back-switching, tihs is Per Protocol time for control patients
    PPtime <- rep(PPtime, each=stime) # expand PP time
    fulldat$PPtime <- ifelse(fulldat$arm == 1 | fulldat$eventtime <= PPtime, fulldat$eventtime, PPtime) # if patient is in experimental group or already sees event before PPtime, use eventtime. Else censor at switch
    fulldat$PPdeath <- ifelse(fulldat$PPtime < fulldat$eventtime, 0, fulldat$status) # within control group, if death had been recorded but is now censored, censor event indicator

    # HR estimate
    pp[r] <- exp(coef(survival::coxph(survival::Surv(PPtime, PPdeath) ~ arm, data = fulldat[fulldat$time == 1,] )))
    # Old, single rep model: survival::coxph(survival::Surv(PPtime, PPdeath) ~ arm, data = fulldat)

    # KM curves:
    pp_plot <- survminer::ggsurvplot(
      fit = survminer::surv_fit(survival::Surv(PPtime, PPdeath) ~ arm, data = fulldat[fulldat$time == 1, ]),
      xlab = "Time",
      ylab = "OS",
      title = "KM Plots for PP", conf.int = TRUE) # %++%
    # geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)

    # Run adjustment models:
    if(verbose > 1){
      print("Performing complex method estimates: IPCW...")
    }

    ## IPCW ##
    if(!ipcw_exclude){
      # TODO write ipcw analysis funciton, e.g., ipcw_analyze()
      fulldat$starttime <- fulldat$time - 1 # add starttime
      fulldat$switch_time <- rep(sapply(unique(fulldat$ids), function(x) ifelse( # get first
        fulldat$arm[fulldat$ids == x][1] == 0 & fulldat$switch_status[fulldat$ids == x][1] == 1, # if id is a control pt. and observes switch,
        fulldat$time[fulldat$ids == x & fulldat$treat == 1][1], # get the first treatment time
        0)), each = stime)
      fulldat$switch <- ifelse(fulldat$time == fulldat$switch_time, 1, 0)
      # cut fulldat at censoring or event time. Recast several fulldat variables
      fulldat_cut <- fulldat[fulldat$time <= fulldat$eventtime,] # cut data past eventtime
      cens_ids <- unique(fulldat_cut$ids[fulldat_cut$switch == 1]) # if you get censored, you are added to cens_ids
      fulldat_cut$status[fulldat_cut$ids %in% cens_ids] <- 0 # if you get censored, your status becomes 0
      fulldat_cut <- fulldat_cut[!(fulldat_cut$switch_status == 1 & fulldat_cut$switch_time < fulldat_cut$time),] # cut data from control arm switcher past the point of switch
      # TODO eventstatus is, for
      fulldat_cut$eventstatus <- ifelse(fulldat_cut$time == floor(fulldat_cut$eventtime) & fulldat_cut$status == 1, 1, 0) # new eventstatus is 1 only if event observed before switch.
      fulldat_cut$time <- as.numeric(fulldat_cut$time)
      fulldat_cut$arm <- as.factor(fulldat_cut$arm)

      # on full data
      if("All" %in% violate | "IPCW" %in% violate){
        ipdat <- ipcw(data = fulldat_cut, id = "ids", tstart = starttime,
                      tstop = time, cens = switch, arm = "arm", bas.cov = bcov_names, conf = tcov_names[1:(length(tcov_names) - hide_tvar)],
                      type = "kaplan-meier", trunc = 0.05)
      }else{
        ipdat <- ipcw(data = fulldat_cut, id = "ids", tstart = starttime,
                      tstop = time, cens = switch, arm = "arm", bas.cov = bcov_names, conf = tcov_names,
                      type = "kaplan-meier", trunc = 0.05)
      }
      # HR estimate
      msm_hr[r] <- exp(coef(survival::coxph(survival::Surv(starttime, time, eventstatus) ~ arm , data = ipdat, weights = ipdat$weights.trunc)))
    }else{
      msm_hr[r] <- as.numeric(NA)
    }

    # make weighted KM estimates
    if(!ipcw_exclude){
      ipcw_plot <- survminer::ggsurvplot(
        fit = survminer::surv_fit(survival::Surv(starttime, time, eventstatus) ~ arm, data = ipdat, weights = ipdat$weights.trunc),
        xlab = "Time",
        ylab = "OS",
        title = "KM Plots for IPCW", conf.int = TRUE) # %++%
      # geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)
    } else{
      ipcw_plot <- NULL
    }

    if(verbose > 1){
      print("Performing complex method estimates: RPSFTM...")
    }
    ## RPSFTM ##

    # TODO write RPSFTM analysis function, e.g., rpsftm_analyze()
    # get proportion of treatment per patient
    rx <- sapply(unique(fulldat$ids), function(x) sum(fulldat$treat[fulldat$ids == x & fulldat$time <= fulldat$eventtime])/
                   fulldat$eventtime[fulldat$ids==x & fulldat$time==1])
    rpsft_dat <- cbind(fulldat[fulldat$time == 1, ], rx) # create condensed RPSFTM dataset
    rpsft_dat$cens <- stime # add an administrative censoring time, which is simply stime
    rpsft_dat$arm <- as.factor(rpsft_dat$arm)

    # single rep mod: Build either recensored or unrecensored model
    if(recens == TRUE){
      # build rpsftm model
      mr <- suppressWarnings(rpsftm::rpsftm(survival::Surv(eventtime, status) ~ rand(arm, rx), data = rpsft_dat, low_psi = -2, hi_psi = 2, censor_time = cens))
    } else{
      mr <- suppressWarnings(rpsftm::rpsftm(survival::Surv(eventtime, status) ~ rand(arm, rx), data = rpsft_dat, low_psi = -2, hi_psi = 2))
    }

    # set counterfactuals with rpsftm model object:
    rpsft_dat$counterfact <- rpsft_dat$eventtime # set default counterfactual
    # TODO whole simuulator SOMETIMES bugs here:
    rpsft_dat$counterfact[rpsft_dat$arm == 0] <- mr$Sstar[rpsft_dat$arm == 0, 1] # get rpsftm counterfactual times
    rpsft_dat$cf_status <- rpsft_dat$status
    rpsft_dat$cf_status[rpsft_dat$arm == 0] <- mr$Sstar[rpsft_dat$arm == 0, 2] # get rpsftm counterfactual death/censoring
    rpsft[r] <- exp(coef(survival::coxph(survival::Surv(counterfact, cf_status) ~ arm, data = rpsft_dat)))

    # make KM estimates, censored
    rpsft_plot <- survminer::ggsurvplot(
      fit = survminer::surv_fit(survival::Surv(counterfact, cf_status) ~ arm, data = rpsft_dat),
      xlab = "Time",
      ylab = "OS",
      title = "KM Plots for RPSFTM", conf.int = TRUE) # %++%
    # geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)

    if(verbose > 1){
      print("Performing complex method estimates: TSE...")
    }
    ## TSE ##

    # TODO write TSE analysis function, i.e., tse_analyze()
    tsdat <- fulldat[fulldat$time == 1,] # Holder for wide format
    tscontrol <- fulldat[fulldat$secbase_observed == 1 & fulldat$arm == 0 & fulldat$Mtime == 1,] # take subset of pts who observe M and who are in arm == 0
    #tscontrol <- tscontrol[tscontrol$time < tscontrol$eventtime,] # subset again, removing pts who observe M (secondary baseline) after eventtime
    tscontrol$TSsurv <- tscontrol$eventtime - tscontrol$time
    tscontrol <- tscontrol[tscontrol$TSsurv > 0,] # exclude last-day-switchers

    TSEform <- formula(paste("survival::Surv(TSsurv, status) ~ switch_status +", paste(c(bcov_names, tcov_names), collapse = "+"), collapse = " "))

    # tse_est <- function(data, indices, tsdat){ # function to pass to boot::boot()
    #   d <- data[indices,] # allows boot to select sample
    #
    #   # fit AF model
    #   mod <- survreg(formula = TSEform, dist = tse_dist, data = d)
    #   AF <- exp(coef(mod))[names(exp(coef(mod))) == "switch_status"] # get acceleration factor
    #   d$counterfact <- ifelse(d$switch_status == 0,
    #                           d$TSsurv + d$time, # observed survival
    #                           (d$TSsurv / AF) + d$time # counterfactual survival
    #   )
    #
    #   #
    #   df_boot <- tsdat[!tsdat$ids %in% data$ids,] # take subset of tsdat which is NOT in tscontrol
    #   df_boot$counterfact <- df_boot$eventtime # reset counterfactuals in df_boot
    #   df_boot <- rbindlist(list(df_boot, d), fill = TRUE) # fold d back into tsdat
    #   # for(i in unique(d$ids)){
    #   #   print(paste("tsdat:", length(tsdat$counterfact[tsdat$ids == i])))
    #   #   print(paste("d: ", length(d$counterfact[d$ids == i])))
    #   #   tsdat$counterfact[tsdat$ids == i] <- d$counterfact[d$ids == i]
    #   # }
    #   fit <- survival::coxph(survival::Surv(counterfact, status) ~ arm, data = df_boot)
    #   return(exp(coef(fit)))
    # }
    # tse_wrapper <- possibly(tse_est, otherwise = NA)
    # tse <- boot::boot(data = tscontrol, statistic = tse_wrapper, R = bootrep, tsdat = tsdat)

    # fit AF model
    mod <- survival::survreg(formula = TSEform, dist = tse_dist, data = tscontrol)

    AF <- exp(coef(mod))[names(exp(coef(mod))) == "switch_status"] # get acceleration factor
    tscontrol$counterfact <- ifelse(tscontrol$switch_status == 0,
                                    tscontrol$TSsurv + tscontrol$time, # observed survival
                                    (tscontrol$TSsurv / AF) + tscontrol$time # counterfactual survival
    )

    tsdat$counterfact <- tsdat$eventtime # reset counterfactuals
    for(i in unique(tscontrol$ids)){
      tsdat$counterfact[tsdat$ids == i] <- tscontrol$counterfact[tscontrol$ids == i]
    }

    tse[r] <- exp(coef(survival::coxph(survival::Surv(counterfact, status) ~ arm, data = tsdat)))

    # make KM estimates, censored
    tse_plot <- survminer::ggsurvplot(
      fit = survminer::surv_fit(survival::Surv(counterfact, status) ~ arm, data = tsdat),
      xlab = "Time",
      ylab = "OS",
      title = "KM Plots for TSE", conf.int = TRUE) # %++%
    # geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)


    # TODO Recensor!



    secondary_baseline_observed[r] <- sum(fulldat$secbase_observed[fulldat$time == 1], na.rm = TRUE) / n
    switch_observed[r] <- sum(fulldat$switch_status[fulldat$time == 1], na.rm = TRUE) / length(id_con)
    Minb0[r] <- sum(fulldat$secbase_observed[fulldat$time == 1 & fulldat$arm == 0 & fulldat$b1==0])
    Minb1[r] <- sum(fulldat$secbase_observed[fulldat$time == 1 & fulldat$arm == 0 & fulldat$b1 == 1])

    if(verbose > 1){
      print(paste("repetition: ", r))
    }

  } # end giant loop. Build plots below.












  # Bootstrap results (in parallel, if possible) with obtained hazard functions ####

  # giant while loop should begin here. At this point, we have all parameters set, and a fulldat dataframe with no time
  # varying covariates, no switching and no secondary baseline
  # TODO if para == TRUE, parallelize here. Can these be done all in parallel (including the 1st rep, where conditions are particular) or is it sequential?
  parallel::mclapply(X = 2:reps, mc.cores = 4, FUN = function(r){

    if(r != 1) unfix <- c("B", "M", "S", "T") # if we have generated the first dataset, unfix all hazards

    while(rerun){ # iteratively update switching hazard function until we get the right proportion. The first term in this check is the number of pts whose final treatment indicator

      for(i in 1:stime){
        # set covars
        fulldat[fulldat$time == i, names(fulldat) %in% c(tcov_names, "Mtime")] <-
          dep(dat = fulldat, window = i, base_var=bcov_names, time_var=tcov_names, covar_coef = covar_coef, m_haz = m_haz, num_t = num_tvar, tcov_n = tcov_names)
        # set treatment indicator
        if(i != 1){ # if its not the first time window
          fulldat$treat[fulldat$time == i] <- ifelse(fulldat$treat[fulldat$time == i-1] == 1, 1, 0) # if the previous window treat is 1, continue treatment
        }
        # if treatment has not yet begun, probability of begining in the next window is a hazard function
        if(m_hard){
          fulldat$treat[fulldat$time == i & fulldat$treat == 0 & fulldat$Mtime > 0 & fulldat$Mtime <= ceiling(m_fidelity*stime)] <-
            stats::rbinom( n = length(fulldat$treat[fulldat$time == i & fulldat$treat == 0 & fulldat$Mtime > 0 & fulldat$Mtime <= ceiling(m_fidelity*stime)]), size = 1,
                           prob = 1 - exp(-exp(log(s_haz[i]) +
                                                 as.matrix(fulldat[fulldat$time == i & fulldat$treat == 0, names(fulldat) %in% c(bcov_names, tcov_names)]) %*% switch_coef))) # randomly assign the treat variable with probability
        }else{
          fulldat$treat[fulldat$time == i & fulldat$treat == 0] <- stats::rbinom( n = length(fulldat$treat[fulldat$time == i & fulldat$treat == 0]), size = 1,
                                                                                  prob = 1 - exp(-exp(log(s_haz[i]) +
                                                                                                        as.matrix(fulldat[fulldat$time == i & fulldat$treat == 0, names(fulldat) %in% c(bcov_names, tcov_names)]) %*% switch_coef))) # randomly assign the treat variable with probability
        }

      }

      # Generate fulldat_cont (the control dataset) by blocking switching
      fulldat_cont$treat <- fulldat_cont$arm
      for(i in 1:stime){
        # set covars
        fulldat_cont[fulldat_cont$time == i, names(fulldat_cont) %in% c(tcov_names, "Mtime")] <-
          dep(dat = fulldat_cont, window = i, base_var=bcov_names, time_var=tcov_names, covar_coef = covar_coef, m_haz = m_haz, num_t = num_tvar, tcov_n = tcov_names)
      }


      if(verbose > 1){
        print("Generating confounded survival times...")
      }

      sdat <- survtime_generator( x = fulldat, hazard = haz, betas = beta.mat, ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat$ids),
                                  b_haz = b_haz)

      if(verbose > 1){
        print("Generating un-confounded survival times...")
      }

      sdat_cont <- survtime_generator( x = fulldat_cont[fulldat_cont$arm == 0,], hazard = haz, betas = beta.mat[fulldat_cont$arm == 0,], ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat_cont$ids[fulldat_cont$arm==0]),
                                       b_haz = b_haz)


      if(verbose > 1){
        print("Adding censoring...")
      }
      # merge datasets:
      if(switch_iter == 0){ # if its the first iteration
        fulldat <- merge(fulldat, sdat, by = "ids")
        fulldat_cont <- merge(fulldat_cont, sdat_cont, by = "ids", all = TRUE)
        # replace experimental group in unconfounded dataset with experimental group in confounded dat. Theoretically, these are generated identically.
        fulldat_cont[fulldat_cont$arm == 1, ] <- fulldat[fulldat$arm == 1, ]
      }else{
        fulldat <- fulldat[, !names(fulldat) %in% c("eventtime", "status")] # first remove old sdat
        fulldat <- merge(fulldat, sdat, by = "ids")

        fulldat_cont <- fulldat_cont[, !names(fulldat_cont) %in% c("eventtime", "status")] # first remove old sdat
        fulldat_cont <- merge(fulldat_cont, sdat_cont, by = "ids", all = TRUE)
        # replace experimental group in unconfounded dataset with experimental group in confounded dat. Theoretically, these are generated identically.
        fulldat_cont[fulldat_cont$arm == 1, ] <- fulldat[fulldat$arm == 1, names(fulldat) %in% names(fulldat_cont)]
      }


      # censor at prespecified, non-administrative censoring time. Censoring will have been generated independently of covars, or dependent on covars
      if(cens_flag == "Random"){
        fulldat$cens <- stime # set default censoring time
        fulldat_cont$cens <- stime # set default censoring time

        cens_ids <- sample(unique(fulldat$ids), ceiling(prop_cens*n)) # sample prop_cens of ids for censoring
        rand_cens <- rbeta(n = length(cens_ids), shape1 = 2, shape2 = 1.5) # censoring times are drawn from a beta dist
        # rand_cens <- ceiling(rand_cens)
        rand_cens <- rep(rand_cens, each=stime) # spread to length of dataset
        fulldat$cens[fulldat$ids %in% cens_ids] <- ceiling(rand_cens*fulldat$eventtime[fulldat$ids %in% cens_ids])
        fulldat_cont$cens[fulldat_cont$ids %in% cens_ids] <- ceiling(rand_cens*fulldat_cont$eventtime[fulldat_cont$ids %in% cens_ids])

        # censor observations. Take minimum of eventtime and cens, and change status to 0 if eventtime == cens
        fulldat$eventtime <- pmin(fulldat$eventtime, fulldat$cens)
        fulldat$status <- ifelse(fulldat$eventtime == fulldat$cens, 0, fulldat$status)
        fulldat_cont$eventtime <- pmin(fulldat_cont$eventtime, fulldat_cont$cens)
        fulldat_cont$status <- ifelse(fulldat_cont$eventtime == fulldat_cont$cens, 0, fulldat_cont$status)


      } else if(cens_flag == "Nonrandom"){
        # TODO nonrandom censoring function!
      }


      # Does the patient (control and experimental) observe secondary baseline (M)?
      # TODO i think pts who observe NO M at all are being given NA here
      fulldat$secbase_observed <- as.numeric(rep(sapply(unique(fulldat$ids), function(x)
        ifelse(sum(fulldat$M[fulldat$ids == x]) > 0 && fulldat$time[fulldat$Mtime == 1 & fulldat$ids == x] <= fulldat$eventtime[fulldat$Mtime == 1 & fulldat$ids == x], # id has at least some m, and time at first m is less than event time
               1,
               0)), each = stime))

      # Does the patient (only control) observe switch _before eventtime_?
      fulldat$switch_status <- rep(sapply(unique(fulldat$ids), function(x)
        ifelse(sum(fulldat$treat[fulldat$ids == x & fulldat$time <= fulldat$eventtime]) > 0 & fulldat$arm[fulldat$ids == x][1] == 0,
               1,
               0)), each = stime)

      if(r == 1){ # if we're in the first loop,
        switch_iter <- switch_iter + 1 # iterate search indicator
      }
      if(switch_iter > rerun_lim) stop("Your covariate model failed to converge. Try different hazards and/or coefficients")

      rerun <- FALSE # pre-empt a rerun, unless the following conditions are met

      # if not enough M occurences, adjust M hazard
      # TODO for now, M occurence is only being adjusted upward. do we want to give it a within-window type adjustment?
      current_m_prop <- length(fulldat$ids[fulldat$arm == 0 & fulldat$Mtime == 1 & fulldat$time <= fulldat$eventtime])
      if(! "M" %in% unfix){
        if(current_m_prop < min(n*prop_trt, prop_switch*prop_trt*n*m_inflation)){ # current_m_prop is proportion of control arm that observes M. if it's less than the min of (the whole control arm, the switchers*m_inflation ), increase m_haz and rerun
          m_haz <- m_haz*m_mag # change scale of m related haz, and rerun
          rerun <- TRUE
        }
      }


      # adjust switching proportion
      current_s_prop <- sum(fulldat$switch_status[fulldat$time == 1 & fulldat$arm == 0])/length(id_con) # get current switch proportion
      if(! "S" %in% unfix){
        if(abs(current_s_prop - prop_switch) > s_allowance){ # If we're outside the window
          prev_s_direc <- s_direc # save most recent s_direc value
          s_direc <- ifelse(current_s_prop > prop_switch, -1, 1) # identify the s_direc, i.e., in what direction the baseline hazard must be adjusted
          if((prev_s_direc == -1 & s_direc == 1) | (prev_s_direc == 1 & s_direc == -1)){ # if we overshot the allowance window, and must backtrack
            s_mag <- s_mag/2 # halve switch magnitude if weve crossed over the window
          }
          s_haz <- s_haz + s_direc*s_mag*s_haz # modify the switch hazard by a factor of s_mag in the direction of s_direc
          rerun <- TRUE
        }
      }


      # adjust control group event occurence
      current_b_prop <- sum(fulldat$status[fulldat$time == 1 & fulldat$arm == 0])/length(id_con) # get current switch proportion
      if(! "B" %in% unfix){
        if(abs(current_b_prop - prop_cont_event) > b_allowance){ # If we're outside the window
          # if(b_haz == "exponential"){
          #   prev_b_direc <- b_direc # save most recent b_direc value
          #   b_direc <- ifelse(current_b_prop > prop_cont_event, -1, 1) # identify the s_direc, i.e., in what direction the baseline hazard must be adjusted
          #   if((prev_b_direc == -1 & b_direc == 1) | (prev_b_direc == 1 & b_direc == -1)){ # if we overshot the allowance window, and must backtrack
          #     b_mag <- b_mag/2 # halve switch magnitude if weve crossed over the window
          #   }
          #   lambdas <- lambdas - b_direc*b_mag*lambdas # modify the switch hazard by a factor of s_mag in the direction of s_direc
          #   if(lambdas <= 0) stop("The scale parameter for the simsurv distribution is 0 or negative; it must be positive.")
          #   rerun <- TRUE
          # }else{
          prev_b_direc <- b_direc # save most recent b_direc value
          b_direc <- ifelse(current_b_prop > prop_cont_event, -1, 1) # identify the s_direc, i.e., in what direction the baseline hazard must be adjusted
          if((prev_b_direc == -1 & b_direc == 1) | (prev_b_direc == 1 & b_direc == -1)){ # if we overshot the allowance window, and must backtrack
            b_mag <- b_mag/2 # halve switch magnitude if weve crossed over the window
          }
          b_haz <- b_haz + b_direc*b_mag*b_haz # modify the switch hazard by a factor of s_mag in the direction of s_direc
          rerun <- TRUE
          # }
        }
      }

      # adjust trt group event occurence
      current_t_prop <- sum(fulldat$status[fulldat$time == 1 & fulldat$arm == 1])/length(id_trt) # get current switch proportion
      if(! "T" %in% unfix){
        if(abs(current_t_prop - prop_trt_event) > t_allowance){ # If we're outside the window
          prev_t_direc <- t_direc # save most recent b_direc value
          t_direc <- ifelse(current_t_prop > prop_trt_event, 1, -1) # identify the t_direc, i.e., in what direction the treatment effect must be adjusted
          if((prev_t_direc == -1 & t_direc == 1) | (prev_t_direc == 1 & t_direc == -1)){ # if we overshot the allowance window, and must backtrack
            t_mag <- t_mag/2 # halve switch magnitude if weve crossed over the window
          }
          beta.mat$treat <- beta.mat$treat + t_direc*t_mag*beta.mat$treat # modify the switch hazard by a factor of s_mag in the direction of s_direc
          rerun <- TRUE
        }
      }

    }

    rerun <- TRUE # reset inner loop rerun flag




    # Run control model:
    unbiased[r] <- get_unbiased()
    unbiased_plot[r] <- get_unbiased_plot()

    # Run Naive models: ####

    ## ITT ##
    itt[r] <- get_itt()
    itt_plot[r] <- get_itt_plot()

    ## Per Protocol ##
    pp[r] <- get_pp()
    pp_plot[r] <- get_pp_plot()

    # Run adjusted models: ####

    ## IPCW ##
    if(!ipcw_exclude){
      ipcw_est[r] <- get_ipcw()
      ipcw_plot[r] <- get_ipcw_plot()
      }else{
      ipcw_est[r] <- as.numeric(NA)
      ipcw_plot[r] <- as.numeric(NA)
    }

    ## RPSFTM ##
    rpsftm_est[r] <- get_rpsftm()
    rpsftm_plot[r] <- get_rpsftm_plot()

    ## TSE ##
    tse[r] <- get_tse()
    tse_plot[r] <- get_tse_plot()





  } # end giant loop. Build plots below.
  )



  # build simulation plots ####
  compar_plot <- get_compar_plot()
  bias_plot <- get_bias_plot()


  if(verbose > 1){
    print("Done!")
  }


  return(sim)


}


