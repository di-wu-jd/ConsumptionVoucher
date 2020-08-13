
# Function to calculate average cost and outcome based on trimmed data by outcome
costs <- function(s, data, trim=0) {
  x <- data[s]
  n <- x[, .N]
  lo <- floor(n * trim) + 1
  hi <- n + 1 - lo
  x <- x[order(outcome)][lo:hi]
  c(x[, sum(outcome)/.N], x[, sum(cost)/.N])
}


# Function to calculate confidence interval on treatment incremental ROI.
# Parameters:
#   c: vector of outcome from control group.
#   t: data.table of outcome and treatment cost from treatment group. It is a data.table with two column (outcome and cost)
#   alpha: significant level.
#   tr: the fraction of observations to be trimmed from the calculation. Same as trim in mean().
#   nboot: number of bootstrap resamples required.
#   ni: batch size of bootstraps to be processed. (Used to reduce memory usage.
#   pr: TRUE to print out progressional information.
#   SEED: TRUE to set random seed for replication.
boot_return <- function(c, t, alpha=.05, tr=NA, nboot=NA, ni=nboot, pr=TRUE, SEED=TRUE) {
  # set ni if not specified in the inputs
  if (is.na(ni)) {
    ni <- nboot
  }

  # calculate overall statistics
  c_stat <- costs(seq(nrow(c)), data=c, trim=tr)
  t_stat <- costs(seq(nrow(t)), data=t, trim=tr)

  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.

  # We calculate two statistics for the treatment group (trimmed mean, and trimmed average discount) and 1 for the control group (trimmed mean)
  # vector of trimmed means for treatment group
  tmean <- numeric()
  # vector of average discount (based on the trimmed data)
  tcost <- numeric()
  # vector of trimmed means for control group
  cmean <- numeric()
  ccost <- numeric()

  # process control group
  if (pr) {
    print("Processing control group.")
  }
  # To avoid memory issues, we only sample n_i replications at a time
  nbootr <- nboot # total number of remaining replications
  nbootc <- min(ni, nbootr) # number of replications to be processed in this iteration.
  while (nbootc > 0) {
    if (pr) {
      print(paste0("Processing bootstrap replications. ", nbootr, " replications remaining."))
    }
    s<-matrix(sample(nrow(c),size=nrow(c)*nbootc,replace=TRUE),nrow=nbootc)
    tmp <- apply(s,1,costs,data=c,trim=tr)
    cmean <- c(cmean, tmp[1, ])
    ccost <- c(ccost, tmp[2, ])
    nbootr <- nbootr - nbootc
    nbootc <- min(ni, nbootr)
  }

  # process treatment group
  if (pr) {
    print("Processing treatment group.")
  }
  # To avoid memory issues, we only sample n_i replications at a time
  nbootr <- nboot # total number of remaining replications
  nbootc <- min(ni, nbootr) # number of replications to be processed in this iteration.
  while (nbootc > 0) {
    if (pr) {
      print(paste0("Processing bootstrap replications. ", nbootr, " replications remaining."))
    }
    s<-matrix(sample(nrow(c),size=nrow(c)*nbootc,replace=TRUE),nrow=nbootc)
    tmp <- apply(s,1,costs,data=t,trim=tr)
    tmean <- c(tmean, tmp[1, ]) # mean
    tcost <- c(tcost, tmp[2, ]) # cost
    nbootr <- nbootr - nbootc
    nbootc <- min(ni, nbootr)
  }

  ret <- (tmean - cmean)/(tcost - ccost)

  lo <- floor(nboot * alpha/2) + 1
  hi <- nboot + 1 - lo

  ret <- sort(ret)
  list(lower=ret[lo], higher=ret[hi], mean=(t_stat[1] - c_stat[1])/(t_stat[2] - c_stat[2]))
}

