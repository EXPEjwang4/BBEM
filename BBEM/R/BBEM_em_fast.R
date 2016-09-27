BBEM_em_fast <-
function(y, X, v0=0.1, v1=100,  maxIter=20, threshold=3){

    # some constants
    n = dim(X)[1]; p = dim(X)[2]
    c1 = 1/v0 - 1/v1
    c2 = log(v1 / v0)
    Xy = t(X) %*%y; 

    # initializing
    nu0=1; lambda0=1; a0=1.1; b0=1.1; 
    theta = min(sqrt(n)/p,a0/(a0+b0));
    sigmaSq=lambda0;
    gamma = rep(0, p); 
    gamma_all = gamma;
    alpha = v0 * (1 - gamma) + v1 * gamma

    if(n > p) M = solve(t(X) %*% X + diag(1/alpha)) else  {
        D = diag(alpha); 
        M = D - D %*% t(X) %*% solve(diag(n) + X %*% D %*% t(X)) %*% X %*% D
    }

    active.var = NULL; stop.flag=0; 

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

        # take expectations over the pdf above					
        postBetaSq = diag(varBeta) + meanBeta^2
        postRSS = sum(diag(X %*% varBeta %*% t(X))) + sum((y - X %*% meanBeta)^2)
       
        D = diag(alpha)

        # M step          
        gamma.tmp= as.numeric(postBetaSq > sigmaSq / c1 * (c2 -2 * log(theta / (1-theta))))                     
        gamma.diff = gamma.tmp - gamma;
        gamma = gamma.tmp
        active.var = (1:p)[abs(gamma.diff)>0];
        alpha = v0 * (1 - gamma) + v1 * gamma
        gamma_all = cbind(gamma_all, gamma);

        theta = (sum(gamma) + a0 -1) / (p + a0 + b0-2)
        sigmaSq = (sum(postBetaSq / alpha) + postRSS + nu0*lambda0) / (n + p + nu0)
        
        # stopping criteria
        if(k > 1){
            if (sum(gamma.diff^2)==0) stop.flag=stop.flag+1; 
            if (stop.flag > threshold) break
        }                   
    }    	
    return(list('gamma'=gamma,'gamma_all'=gamma_all,'theta'=theta,'sigmaSq'= sigmaSq,'numOfIter'=k))
}
