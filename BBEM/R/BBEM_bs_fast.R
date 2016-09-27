BBEM_bs_fast <-
function(y, X, n.var, v0=0.1, v1=100, maxIter=20, T=100, threshold=3){

  # some constants
  n = dim(X)[1]; old.p = dim(X)[2]
  c1 = 1/v0 - 1/v1
  c2 = log(v1 / v0)
  old.X = X;

  # input hyperparameters
  nu0=1; lambda0=1; a0=1.1; b0=1.1;

  gamma_all = matrix(NA, T, p); theta_all=rep(0,T); sigmaSq_all=rep(0,T); selected=rep(0,p)
  var.weights = abs(t(X) %*% y /diag(t(X) %*% X))
  var.weights = var.weights/sum(var.weights)

  for(t in 1:T){
    var.id=sample(1:(old.p), n.var, prob=var.weights); 
    X = old.X[, var.id];
    p=dim(X)[2];
    
    # generate bs weights
    w=rgamma(n,1,1); 
    w=w/sum(w)*n;
    Xy = t(X) %*%(w*y);
  
    # initializing
    theta = min(sqrt(n)/p,a0/(a0+b0));
    sigmaSq=lambda0;
    gamma = rep(0, p); 
    alpha = v0 * (1 - gamma) + v1 * gamma

    if(n > p) M = solve(t(X) %*% diag(w) %*% X + diag(1/alpha)) else  {
      D = diag(alpha); 
      M = D - D %*% t(X) %*% solve(diag(1/w) + X %*% D %*% t(X)) %*% X %*% D
    }

    active.var=NULL; stop.flag=0

    for(k in 1:maxIter){
                                        # E step
        # find pdf of beta given last iter parameters and data
      if (length(active.var)>1) {
        M.active = M[, active.var];
        tmp=solve(diag((1-2*gamma[active.var])/c1) + M.active[active.var,]); 
        M = M - M.active %*% tmp %*% t(M.active)
      }

      if (length(active.var)==1) {
        M.active = M[, active.var];
        tmp=1/((1-2*gamma[active.var])/c1+ M.active[active.var]);
        M = M - tmp*M.active %*% t(M.active)
      }

      meanBeta = M %*% Xy
      varBeta = sigmaSq * M
        
      postBetaSq = diag(varBeta) + meanBeta^2
      postRSS = sum(w*diag(X %*% varBeta %*% t(X))) + sum((y - X %*% meanBeta)^2*w)

      gamma.tmp= as.numeric(postBetaSq > sigmaSq / c1 * (c2 -2 * log(theta / (1-theta))))                
     
      gamma.diff = gamma.tmp - gamma;
      gamma = gamma.tmp
      active.var = (1:p)[abs(gamma.diff)>0];
      alpha = v0 * (1 - gamma) + v1 * gamma 
      theta = (sum(gamma) + a0-1) / (p + a0 + b0-2)
      sigmaSq = (sum(postBetaSq / alpha) + postRSS + nu0*lambda0) / (n + p + nu0)

        # stopping criteria
      if(k > 1){
        if (sum(gamma.diff^2)==0) stop.flag=stop.flag+1; 
        if (stop.flag > threshold) break
      }   
    }
    gamma_all[t,var.id] = gamma;
    selected[var.id] = selected[var.id]+1
    theta_all[t] = theta
    sigmaSq_all[t] = sigmaSq
  }
  return(list('gamma' = rowMeans(t(gamma_all), na.rm=TRUE),'theta'=mean(theta_all),'sigmaSq'=mean(sigmaSq_all),'selvar'=selected))
}
