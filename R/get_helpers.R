get_empty_estimate_frame <- function(){
  data.frame(
    unbiased = NA,
    itt = NA,
    pp = NA,
    ipcw_est = NA,
    rpsftm_est = NA,
    tse = NA,
    secondary_baseline_observed = NA,
    switch_observed = NA,
    switch_iter = NA, # how many times have we searched on the inner loop?
    Minb0 = NA,
    Minb1 = NA
    )
}


get_surv_params <- function(){

}

checkparams <- function(){

}

checksurvparams <- function(){

}
