run_replicate <- function(replicate){
  require(glmnet)
  require(randomForest)
  require(mvtnorm)
  require(assertthat)
  
  Z.fixed = list(
    # factors here that are not going to vary, but still must be set
    n.test = 10000
  )
  
  # notation change below makes this consistent with Bingham and Chipman
  Z.levels = list(
    # replicate runs
    replicate = replicate,  # to add replicates later, rerun with bigger values (e.g., 1:3, then 4:6) 
    # external settings
    n = c(250, 1000),   # sample size of training set
    # data generating process
    q = c(20, 50),      # number of measured variables
    ENE = c(10, 20),    # number of terms in the model (including interactions)
    beta.mu = c(1, 3),  # coefficient values for nonzero terms
    sigma = c(0.5,2),   # noise level 
    x.cor = c(0,0.8),   # collinearity of predictors
    # method 
    model = c("lasso","rf")
  )
  
  # make a full factorial in the levels...
  Z <- expand.grid(Z.levels)
  
  ## Seeds are built so that:
  # 1. Every run gets a different seed.  i.e. no randomization restrictions.
  # 2. The seed numbers begin with 10,000*(replicate index), 
  #    eg for 3 replicates, they'd be {10000, 20000, 30000} + base seed.
  #    This should make it possible to add replicates without re-running previous runs.
  #    For example in run 1, you might use replicate = 1:3, and then later use replicate=4:6
  Z.collapsed <- apply(Z[,-1],1,paste,collapse='.')  # CHANGED from other case
  Z.collapsed.u <- unique(Z.collapsed)
  Z.seed <- rep(NA,length(Z.collapsed))
  
  for(i in 1:length(Z.collapsed.u)){
    whch <- Z.collapsed == Z.collapsed.u[i] 
    Z.seed[whch] <- i + 10000*Z[whch,"replicate"]
  }
  
  results <- matrix(0,nrow(Z),10)
  colnames(results) <- c("train.oracle.rsq", "model.rsq.test",
                         "model.rsq.Ey.test","model.mse.test",
                         "meanSStot","nme","n2fi1par","n2fi2par",
                         "execution.time", "seed")
  rfopt <- matrix(0,nrow(Z),2,dimnames = list(NULL,c("mtry.pct","ntree")))
  
  # code could be modified to run different rows of the Z matrix in parallel
  for (run in 1:nrow(Z)){
    
    cat(replicate,":",run, "\n")
    gc() # garbage collection before starting gives reliable timings
    start.time <- Sys.time()
    
    set.seed(Z.seed[run])
    n <- Z$n[run]
    q <- Z$q[run]
    ENE <- Z$ENE[run]
    beta.mu <- Z$beta.mu[run]
    sigma <- Z$sigma[run]
    x.cor <- Z$x.cor[run]
    model <- Z$model[run]
    
    # Generate X matrix
    xcormat <- matrix(0,q,q)
    for (i in 1:q) {
      for (j in i:q){
        if (x.cor){
          xcormat[i,j] <- xcormat[j,i] <- 
            exp(log(0.8)*abs(i-j))  # 0.8 will be correlation of adjacent x vars.
          
        } else xcormat[i,j] <- (i==j)
      }
    }
    
    X <- rmvnorm(n, sigma = xcormat)
    X.test <- rmvnorm(Z.fixed$n.test, sigma = xcormat)
    colnames(X) <- colnames(X.test) <- paste("X",1:ncol(xcormat),sep='')
    
    # from Bingham and Chipman
    # with a 2fi prior on activity given as 
    # 
    # P(A) = P(main effect A active) = p
    # P(AB | A, B) = c1*p if both A, B inactive  [ we choose c1=0 ]
    #              = c2*p if exactly 1 of A, B active
    #              = c3*p if both A, B active
    # 
    # We specify ENE = expected number of effects as a number (level in Z).
    # Assume we want to divide ENE among expected number of main effects, 
    # expected number of 2fi /w 1 parent and expected number of 2fi /w 2 parents
    # according to proportions a1+a3+a2=1
    a1 = 0.5
    a3 = 0.25
    a2 = 1 - a1 - a3
    
    # See "expected heredity calculations.pdf". 
    p = a1*ENE/q
    c3 =  2*a3*ENE/(p^3*q*(q-1))
    c2 = a2*ENE/(p^2*q*(1-p)*(q-1)) 
    
    # check our settings are OK
    assert_that(p >= 0 & p <= 1)
    assert_that(c2*p >= 0 & c2*p <= 1)
    assert_that(c3*p >= 0 & c3*p <= 1)
    
    # sample the active main effects
    me.active <- runif(q) < p
    while(sum(me.active)==0){ # if there are 0 active ME's roll the dice again.
      me.active <- runif(q) < p
    }
    indices.of.active.me <- which(me.active)
    nme = length(indices.of.active.me)
    indices.of.inactive.me <- which(!me.active)
    nime = length(indices.of.inactive.me)
    
    # set up null cases, in case there are no ints of either kind
    indices.of.twofi.2active.parent <- 
      indices.of.twofi.1active.parent <- matrix(0,nrow=0,ncol=2)
    
    if (nme > 0){
      # conditional on active main effects, use weak heredity to sample ints
      
      # terms with 1 active parent
      if (nime > 0) {
        n1par = rbinom(n=1, size = nime*nme, prob= c2*p) # sample number of terms
        tmp1 = as.matrix(expand.grid(indices.of.active.me,indices.of.inactive.me))
        # given number of terms, sample pairs uniformly
        tmp1ind <- sample(nrow(tmp1),n1par) 
        indices.of.twofi.1active.parent <- tmp1[tmp1ind,,drop=FALSE]
        rm(tmp1)
      } 
      
      # similar for 2 active parents
      if (nme > 1){
        n2par = rbinom(n=1, size = choose(nme,2), prob = c3*p)
        tmp2 = matrix(0,choose(nme,2),2)
        counter <- 1
        for (i in 1:(nme-1)){
          for (j in (i+1):nme){
            tmp2[counter,]=c(indices.of.active.me[i],indices.of.active.me[j])
            counter <- counter + 1
          }
        }
        
        tmp2ind <- sample(nrow(tmp2),n2par)
        indices.of.twofi.2active.parent <- tmp2[tmp2ind,,drop=FALSE]
        rm(tmp2)
      }
    } 
    
    # Note: interaction columns are formed by multiplying 2 columns of X. These
    # product columns do not appear to need scaling.  I empirically verified
    # that the mean and sd of the product terms are still reasonably close to 0
    # and 1, when the correlation between Xs is decaying from 0.8.  The sample
    # sd's did vary more from product to product, but were still centered around
    # 1. If we were to worry about scaling, I think this paper asserts a
    # distributional result for the product of 2 normal RVs:
    # https://arxiv.org/pdf/1807.03981.pdf
    
    # generate the beta values and calculate E(Y|X) for train and test sets. Also
    # make an "oracle" regression matrix containing columns for only the active
    # effects.
    
    # Setting sd for the beta draws to 0 yields all nonzero betas = beta.mu
    betasd = 0 
    
    muvec <- rep(0,n)
    muvec.test <- rep(0,Z.fixed$n.test)
    truebeta1 <- rnorm(nme,beta.mu,betasd)
    for(i in 1:nme){
      muvec <- muvec + truebeta1[i]*X[,indices.of.active.me[i]]
      muvec.test <- muvec.test + truebeta1[i]*X.test[,indices.of.active.me[i]]
    }
    X.oracle <- X[,indices.of.active.me]
    
    n2fi1par <- nrow(indices.of.twofi.1active.parent)
    if (n2fi1par){ # if there are any active 2fi /w 1 parent terms...
      truebeta2 <- rnorm(n2fi1par,beta.mu,betasd)  
      for(i in 1:n2fi1par){
        muvec <- muvec + 
          truebeta2[i]*X[,indices.of.twofi.1active.parent[i,1]]*
          X[,indices.of.twofi.1active.parent[i,2]]
        muvec.test <- muvec.test + 
          truebeta2[i]*X.test[,indices.of.twofi.1active.parent[i,1]]*
          X.test[,indices.of.twofi.1active.parent[i,2]]
        X.oracle <- cbind(X.oracle,X[,indices.of.twofi.1active.parent[i,1]]*
                            X[,indices.of.twofi.1active.parent[i,2]])
        colnames(X.oracle)[ncol(X.oracle)] <- 
          paste(colnames(X)[indices.of.twofi.1active.parent[i,]],collapse=':')
      }
    }
    n2fi2par <- nrow(indices.of.twofi.2active.parent)
    if (n2fi2par){ # if there are any active 2fi /w 2 parent terms...
      truebeta3 <- rnorm(n2fi2par,beta.mu,betasd)
      for(i in 1:n2fi2par){
        muvec <- muvec + 
          truebeta3[i]*X[,indices.of.twofi.2active.parent[i,1]]*
          X[,indices.of.twofi.2active.parent[i,2]]
        muvec.test <- muvec.test + 
          truebeta3[i]*X.test[,indices.of.twofi.2active.parent[i,1]]*
          X.test[,indices.of.twofi.2active.parent[i,2]]
        X.oracle <- cbind(X.oracle,X[,indices.of.twofi.2active.parent[i,1]]*
                            X[,indices.of.twofi.2active.parent[i,2]])
        colnames(X.oracle)[ncol(X.oracle)] <- 
          paste(colnames(X)[indices.of.twofi.2active.parent[i,]],collapse=':')
      }
    }
    
    # add noise to response...
    y <- muvec + rnorm(n,0,sigma)
    y.test <- muvec.test + rnorm(Z.fixed$n.test,0,sigma)
    
    if (model=="lasso"){
      fullx <- model.matrix(~.^2, data = data.frame(X))[,-1]
      fullx.test <- model.matrix(~.^2, data = data.frame(X.test))[,-1]
      junk <- cv.glmnet(x = fullx, y = y, nfolds = 10)
      yhat.test <- predict(junk, newx = fullx.test, s = "lambda.min")
    }
    
    if(model=="rf"){
      
      # Note: in 128 runs, all but 6 chose largest "maxnodes" and
      # all but 3 chose ntree=500.  So I am collapsing the search to one over mtry
      rfsettings <- expand.grid( # grid over which to search for best settings
        mtry.pct = c(0.05,0.10,0.20,0.30,0.40,0.50,0.60,0.70,1),
        ntree = 500
      )
      
      rfres <- rep(0,nrow(rfsettings))
      
      for (ii in 1:nrow(rfsettings)){ # find best settings using OOB preds
        junk <- randomForest(x = X, y = y, 
                             mtry = floor(q*rfsettings$mtry.pct[ii]),
                             ntree = rfsettings$ntree[ii])
        rfres[ii] <- cor(y,predict(junk))
      }
      best.ii <- which.max(rfres) # best OOB performance
      junk <- randomForest(x = X, y = y, 
                           mtry = floor(q*rfsettings$mtry.pct[best.ii]),
                           ntree = rfsettings$ntree[best.ii])  # refit best
      rfopt[run,"mtry.pct"] <- rfsettings$mtry.pct[best.ii]
      rfopt[run,"ntree"] <- rfsettings$ntree[best.ii]
      
      yhat.test <- predict(junk, newdata = X.test)
    }
    end.time <- Sys.time()
    
    results[run,] <- c(
      summary(lm(y~X.oracle))$r.sq,  # train Rsq from estimated oracle model
      cor(yhat.test, y.test)^2,    # test Rsq from best model
      cor(yhat.test, muvec.test)^2, # test Rsq using E(Y|x) instead of Y.test
      mean((yhat.test - y.test)^2), # test MSE from best lasso model
      mean((y.test - mean(y.test))^2), # mean SS from normal error model
      nme,n2fi1par,n2fi2par, 
      as.numeric(end.time - start.time, units = "secs"),
      Z.seed[run]
    )
    
  }
  
  results <- cbind(Z,results)
  save(results,file = paste("replicate_",replicate,"_results.RData",sep=""))
  
}
