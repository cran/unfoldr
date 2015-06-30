data(spheroids)
## set number of cpu cores (optional)
# options(par.unfoldr=8)

## Setup spheroid system 
setupSpheroidSystem(spheroids,box=list(c(0,5)))

## unfolding
sp <- verticalSection(spheroids,2.5)
ret <- unfold(sp,c(10,9,8),kap=1.3)
