# Revised version of linconpb function from WRS package to add a parameter to enable memory efficient calculation for large dataset.
# Code of original linconpb function as well as installation instruction of 'WRS' can be found at https://github.com/nicebread/WRS

linconpba<-function(x,alpha=.05,nboot=NA,grp=NA,est=mean,con=0,bhop=FALSE,
                   ni=nboot,pr=TRUE,SEED=TRUE,...){
  #
  #   Multiple comparisons for  J independent groups.
  #
  #   The data are assumed to be stored in x
  #   which either has list mode or is a matrix.  In the first case
  #   x[[1]] contains the data for the first group, x[[2]] the data
  #   for the second group, etc. Length(x)=the number of groups = J.
  #   If stored in a matrix, the columns of the matrix correspond
  #   to groups.
  #
  #   est is the measure of location and defaults to a 20% trimmed mean
  #   ... can be used to set optional arguments associated with est
  #
  #   The argument grp can be used to analyze a subset of the groups
  #   Example: grp=c(1,3,5) would compare groups 1, 3 and 5.
  #
  #   Missing values are allowed.
  #
  # ni specifies the maximum number of replications to be processed at one time.
  #     specifying ni in the parameter allows the function to use bounded
  #     memory during the bootstraping step to enable processing of large data
  #     set without running out of memory.
  # IF pr=TRUE, some information will be printed along with the execution.
  # If bhop=FALSE
  # Probability of one or more Type I errors controlled using Hochberg's method.
  # If bhop=TRUE, the Benjamini--Hochberg method is used. 
  #
  # When using a onestep M-estimator or mom, pbmcp is a good choice, which uses method SR
  # in section 7.6.2 of Wilcox, 2012, Introduction to Robust Estimation and Hypothesis Testing
  #
  con<-as.matrix(con)
  if(is.matrix(x))x<-listm(x)
  if(!is.list(x))stop('Data must be stored in list mode or in matrix mode.')
  if(!is.na(sum(grp))){  # Only analyze specified groups.
    xx<-list()
    for(i in 1:length(grp))xx[[i]]<-x[[grp[i]]]
    x<-xx
  }
  J<-length(x)

  # sample esitmations by group
  tempn<-0
  mvec<-NA
  for(j in 1:J){
    temp<-x[[j]]
    temp<-temp[!is.na(temp)] # Remove missing values.
    tempn[j]<-length(temp)
    x[[j]]<-temp
    mvec[j]<-est(temp,...)
  }
  nmax <- max(tempn)
  Jm<-J-1
  #
  # Determine contrast matrix
  #
  if(sum(con^2)==0){
    ncon<-(J^2-J)/2
    con<-matrix(0,J,ncon)
    id<-0
    for (j in 1:Jm){
      jp<-j+1
      for (k in jp:J){
        id<-id+1
        con[j,id]<-1
        con[k,id]<-0-1
      }}}
  ncon<-ncol(con)
  if(nrow(con)!=J){
    stop('Something is wrong with con; the number of rows does not match the number of groups.')
  }
  #  Determine nboot if a value was not specified
  if(is.na(nboot)){
    nboot<-5000
    if(J <= 8)nboot<-4000
    if(J <= 3)nboot<-2000
  }
  # Determine critical values
  if(bhop)dvec<-(ncon-c(1:ncon)+1)*alpha/ncon
  if(!bhop)dvec<-alpha/c(1:ncon)
  bvec<-matrix(NA,nrow=J,ncol=nboot)
  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  
  # set ni if not specified in the inputs
  if (is.na(ni)) {
    ni <- nboot
  }
  est.vec<-NA
  for(j in 1:J){
    est.vec[j]<-est(x[[j]],...)
    # To avoid memory issues, we only sample n_i replications at a time
    nbootr <- nboot # total number of remaining replications
    nbootc <- min(ni, nbootr) # number of replications to be processed in this iteration.
    bvecc <- numeric()
    while (nbootc > 0) {
      if (pr) {
        print(paste0("Processing bootstrap replications. ", nbootr, " replications remaining."))
      }
      data<-matrix(sample(x[[j]],size=length(x[[j]])*nbootc,replace=TRUE),nrow=nbootc)
      bvecc <- c(bvecc, apply(data,1,est,...)) # bootstrapped values in the current iteration
      nbootr <- nbootr - nbootc
      nbootc <- min(ni, nbootr)
    }
    bvec[j,]<-bvecc # Bootstrapped values for jth group
    rm(data) # free up memory
  }
  if (pr) {
    print("All samples generated, aggregating results.")
  }

  # calculate p-values
  test<-NA
  bcon<-t(con)%*%bvec #ncon by nboot matrix
  tvec<-t(con)%*%mvec
  for (d in 1:ncon){
    test[d]<-(sum(bcon[d,]>0)+.5*sum(bcon[d,]==0))/nboot
    if(test[d]> .5)test[d]<-1-test[d]
  }
  test<-2*test
  output<-matrix(0,ncon,6)
  dimnames(output)<-list(NULL,c('con.num','psihat','p.value','p.crit',
                                'ci.lower','ci.upper'))
  temp2<-order(0-test)
  zvec<-dvec[1:ncon]
  sigvec<-(test[temp2]>=zvec)
  output[temp2,4]<-zvec
  icl<-round(dvec[ncon]*nboot/2)+1
  icu<-nboot-icl-1
  for (ic in 1:ncol(con)){
    output[ic,2]<-tvec[ic,]
    output[ic,1]<-ic
    output[ic,3]<-test[ic]
    temp<-sort(bcon[ic,])
    output[ic,5]<-temp[icl]
    output[ic,6]<-temp[icu]
  }
  num.sig<-sum(output[,3]<=output[,4])
  list(output=output,con=con,est.loc=est.vec,num.sig=num.sig)
}