\dontrun{
	
lam <- 10
box <- list("xrange"=c(0,3),"yrange"=c(0,3),"zrange"=c(0,9))

# Spheroids of constant sizes
theta <- list("size"=list(0.25),
			  "shape"=list("radius"=0.15),
			  "orientation"=list("kappa"=1))

S <- simCylinderSystem(theta,lam,size="const",
						orientation="rbetaiso",box=box,pl=101)

#require("rgl")
#open3d()
#cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
#cylinders3d(S, box, col=cols)
  
}
