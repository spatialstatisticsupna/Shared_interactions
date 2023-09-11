inla.rgeneric.SCM = function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL){
  
  envir = parent.env(environment())
  
  interpret.theta = function(){
    return(
      list(prec = exp(theta[1L]),
           a = sapply(theta[as.integer(2:(k+1))], function(x){exp(x)})
      )
    )
    
  }
  
  graph = function(){
    return(Q())
  }
  
  Q = function(){
    param = interpret.theta()
    
    C0 <- Diagonal(x = rep(1,nrow(W)))
    ##RW1
    D1 <- diff(diag(k), differences=1)
    R <- t(D1)%*%D1
    C1 <- kronecker(R,C0)
    
    delta <- kronecker(solve(Diagonal(x = param$a)),Diagonal(x = rep(1,nrow(W))))
    delta2 <- delta %*% delta
    prec_z <- delta %*% C1 %*% delta
    
    C <- bdiag(param$prec * prec_z + 
                 exp(15) * delta2 %*% C1 %*% delta2,
               exp(15) * C1)
    C[1:(k*nrow(W)),(k*nrow(W)+1):(nrow(W)*k*2)] <- - exp(15) * delta2 %*% C1
    C[(k*nrow(W)+1):(nrow(W)*k*2), 1:(k*nrow(W))] <- - exp(15) * C1 %*% delta2
    
    Q <- inla.as.sparse(C)
    return(Q)
  }
  
  mu = function(){
    return(numeric(0))
  }
  
  log.norm.const = function(){
    return(numeric(0))
  }
  
  log.prior = function(){
    param = interpret.theta()
    
    res <- dgamma(param$prec, 1, 5e-05, log = TRUE) + theta[1L] +
      sum(dgamma(param$a, 10, 10, log = TRUE)) + sum(theta[2:(k+1)])
    
    return(res)
  }
  
  initial = function(){
    return(as.vector(initial.values))
  }
  
  quit = function(){
    return(invisible())
  }
  
  if (!length(theta)){
    theta <- initial()
  }
  val <- do.call(match.arg(cmd), args = list())
  return(val)
}