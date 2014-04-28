#'MODEL FIT MEASURES
#' 
modelfit_mg<-function(Z0, W, A, nvar, nlv, ng, case_index){
#MODEL FIT MEASURES
gfi_1 <- 0
gfi_2 <- 0  
srmr_1 <- 0
kk <- 0
ss <- 0
COR_RES <- list() # correlation residual matrix 
  for (g in 1:ng){
    k <- kk + 1
    kk <- kk + nvar
    s <- ss + 1
    ss  <- ss + nlv
   
    zz <- Z0[case_index[g,1]:case_index[g,2], k:kk]#HELA: case_index, what are the realistic dimensions?
    w <- W[k:kk, s:ss]
    v <- cbind(diag(nvar), w)
    t <- A[s:ss, ]
    omega <- v - w%*%t
    ee <- zz%*%omega
   
    samcov <- cov(zz)  #sample covariances for each group
    samcorr <- cor(zz)
    precov <- solve(omega%*%t(omega))%*%omega%*%diag(apply(ee,2,var))%*%t(omega)%*%solve(omega%*%t(omega))
   
   browser()
   
  }
res<-list()     
return(res)
}