#' A function to estimate parameters in GCSA.
#' 
#' N = number of observations
#' Ngen = number of observed Genotypes
#' Nphen = number of observed Phenotypes
#' J = number of observed variables (Genotypes+Phenotypes)=Ngen+Nphen
#' P = number of latent variables  (Genes+Clinical pathways)
#' T = total number of variables (= J + P)
#' f0 = 100000000;
#' imp = 100000000;
#' V = [eye(J), W]
#' A = [C,B]
#' No Kronecker products used in estimating all parameters (W,C,& B)
#'
#' Ridge penalties are added (no regularization)
#' lambda_w = ridge parameter for W   default 0
#' lamdbda_b = ridge parameter for B  default 0
#' if lambda_w = lambda_b = 0 (no regularization)
#'    
#' @param  N = number of observations
#' @param Z0 = N x J matrix of observed variables: Genotypes first then Phenotypes
#' @param W0 = J x P matrix with 0's and 1's. W0[j,p]=1 indicates an arrow from the observed variable j to the latent variable p.
#' @param B0 = P x P matrix with 0's and 1's. B0[p1,p2]=1 indicates an arrow from the latent variable p1 to the latent variable p2.
#' @param lambda_w = ridge parameter for W   default 0
#' @param lamdbda_b = ridge parameter for B  default 0
#'
#' @return \item{Westim}{matrix with weight coefficients estimates}
#' @return \item{B.estim}{matrix with path coefficients estimates}
#' @return \item{vecW.estim}{vector with the weight coefficients estimates}
#' @return \item{vecB.estim}{vector with the path coefficients estimates}
#' @return \item{FIT,FIT_M,FIT_S,AFIT,GFI,SRMR}{??????????}
#'
#' @seealso \code{\link{????}}, \code{\link{???????}}, \code{\link{?????}}
#'
#' @references
#'
#' GRAMMAR-Raw:
#' Aulchenko YS, de Koning DJ, Haley C.
#' Genomewide rapid association using mixed model and regression: a fast and
#' simple method for genomewide pedigree-based quantitative trait loci
#' association analysis. Genetics. 2007 Sep;177(1):577-85.
#'
#' GRAMMAR-GC:
#' Amin N, van Duijn CM, Aulchenko YS.
#' A genomic background based method for association analysis in related individuals.
#' PLoS One. 2007 Dec 5;2(12):e1274.
#'
#' GRAMMAR-Gamma:
#' Svischeva G, Axenovich TI, Belonogova NM, van Duijn CM, Aulchenko YS.
#' Rapid variance components-based method for whole-genome association analysis.
#' Nature Genetics. 2012 44:1166-1170. doi:10.1038/ng.2410
#'
#' @examples
#' # Using clean ge03d2 data
#' require(GenABEL.data)
#' data(ge03d2.clean)
#' # take only a small piece for speed
#' ge03d2.clean <- ge03d2.clean[1:200,]
#'
#' @author Hela, Stepan, .....
#'
#' @keywords GSCA

GSCAestim<-function(Z0,W0,B0,lambda_w = 0,lambda_b = 0){

  N <- nrow(Z0)
  
  Z0 <- as.matrix(Z0);W0 <- as.matrix(W0);B0 <- as.matrix(B0);
  
  J <- nrow(W0); P <- ncol(W0)
  
  JP <- J + P
  
  W <- W0
  B <- B0
  C <- matrix(0, P, J)
  
  for (p in c(1:P)){
    WINDp <- which(W0[ , p] != 0)  #find the non null elements
    print(WINDp)
    if (length(WINDp) > 1){    
      W0[WINDp, p] <- 99}
    }
  
  B0[B0 == 1] <- 99
  
  WIND <- which(W0 == 99, arr.ind = T)
  
  BIND <- which(B0 == 99, arr.ind = T)
  
  W[WIND] <- runif(nrow(WIND))         
  B[BIND] <- runif(nrow(BIND))
  
  V <- cbind(diag(J), W)
  
  Z <- scale(Z0)*sqrt(N)/sqrt(N-1) 
  if (rankMatrix(t(Z)%*%Z) == J){
    Z <- chol(t(Z)%*%Z)}
  sizeZ <- dim(Z)[1]
  Gamma <- Z%*%W
  Psi <- Z%*%V
  
  f0 <- 100000000
  imp <- 100000000
  
  it <- 0             
  
    while(imp > 0.000001){print(imp)
      it <- it+1
      #Step 1: Update B
      tr_b <- 0
        for (p in 1:P){ print(p)
          ee <- matrix(0, 1, P)
          ee[p] <- 1
          LL <- diag(P)
          LL[p, p] <- 0
          b0 <- B0[, p]
          bindex_p <- which(b0 == 99)
          YY <- Gamma - Gamma%*%B%*%LL
          if (length(bindex_p)!=0){
            B[bindex_p, p] <- solve(t(Gamma[ , bindex_p])%*%Gamma[ , bindex_p] + lambda_b*diag(length(bindex_p)))%*%((t(Gamma[ , bindex_p])%*%YY)%*%t(ee));
            tr_b <- tr_b + t(B[bindex_p, p])%*%B[bindex_p, p]
          }
        }
      browser()
    }

}


