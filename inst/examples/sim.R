\dontrun{
# directional distribution
rbetaiso <- function(kappa) {
   phi <- runif(1,0,1)*2*pi
   q <- runif(1,0,1)
   theta=acos((1-2*q)/sqrt(kappa*kappa-(1-2*q)*(1-2*q)*(kappa*kappa-1)))
   list("u"=c(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)),
		   "theta"=theta,"phi"=phi)					
}
   
# multivariate size distribution and orientation distribution 
rmulti <- function(m,s,kappa) {	
   dir <- rbetaiso(kappa)
   M <- chol(s, pivot = TRUE)
   M <- M[, order(attr(M, "pivot"))]
   x <- exp(matrix(m,nrow=1) +
          matrix(rnorm(ncol(s)), nrow = 1, byrow = TRUE) %*%M)
   a <- min(x)
   b <- max(x)
   
   list("a"=a,"b"=b,"u"=dir$u,"shape"=a/b,
        "theta"=dir$theta, "phi"=dir$phi)

}

set.seed(1234)
sigma <- matrix(c(0.1,0.1,0.1,0.25), ncol=2)
theta <- list("lam"=500,"rmulti"=list("m"=c(-3.0,-2.0),"s"=sigma,"kappa"=0.5))
S <- simSpheroidSystem(theta,rjoint="rmulti",box=list(c(0,5)),pl=101)

# Spheroids with lognormal distributed length of major axis
theta <- list("lam"=1000,"size"=list("meanlog"=-2.5,"sdlog"=0.5),
              "shape"=list("s"=0.5),
              "orientation"=list("kappa"=1.5))
S <- simSpheroidSystem(theta,size="rlnorm",orientation="rbetaiso",
                       box=list(c(0,5)),pl=101)

# Spheroids of constant sizes
theta <- list("lam"=1000,"size"=list(0.25),"shape"=list("s"=0.5),
              "orientation"=list("kappa"=1))
S <- simSpheroidSystem(theta,size="const",orientation="rbetaiso",
                       box=list(c(0,5)),pl=101)
}
