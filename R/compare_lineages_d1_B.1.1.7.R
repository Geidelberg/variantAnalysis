# Maximum likeihood estimation of Ne(t) from resampled ML trees
# Designed to be run on CLIMB server to retrieve sample time information

library( ape )
library( lubridate )
library( glue )
library( mlesky )
library( treedater ) 

library( sarscov2 ) 
require(ggplot2)
require(grid)
require(gridExtra)
require(ggtree)
library(alakazam)
require(stringi)



#95% HPD interval	[5.2065E-4, 6.7144E-4]
mr = 5.9158E-4
mrci = 	c( 5.2065E-4, 6.7144E-4)
mrsd = diff( mrci ) / 4 / 1.96





# metadata 
civetfn =  list.files(  '/cephfs/covid/bham/climb-covid19-volze/phylolatest/civet/' , patt = 'cog_global_[0-9\\-]+_metadata.csv', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
civmd = read.csv( civetfn , stringsAs=FALSE , header=TRUE )
civmd$central_sample_id <-  sapply( strsplit( civmd$sequence_name , split='/' ) , '[', 2 ) # for linkage 
civmd$sample_date <- as.Date( civmd$sample_date )
civmd$sample_time <- decimal_date( civmd$sample_date ) 



datetree <- function(mltr, civmd, meanrate)
{
  sts <- setNames( civmd$sample_time[  match( mltr$tip.label , civmd$central_sample_id ) ], mltr$tip.label )
  tr <- di2multi(mltr, tol = 1e-05)
  tr = unroot(multi2di(tr))
  tr$edge.length <- pmax(1/29000/5, tr$edge.length)
  dater(unroot(tr), sts[tr$tip.label], s = 29000, omega0 = meanrate, numStartConditions = 0, meanRateLimits = c(meanrate, meanrate + 1e-6), ncpu = 6)
}



# bootrep <- function(mltr, civmd = civmd, meanrate = meanrate, taxis = taxis, ...)
# {
#   td = datetree( mltr , civmd = civmd, meanrate = meanrate ) 
#   res = diff( range( epiweek( date_decimal( td$sts )  )  )  ) + 1
#   res <- res * 2
#   sg = mlskygrid( td, tau = NULL, tau_lower=.001, tau_upper = 10 , sampleTimes = td$sts , res = res, ncpu = 6, NeStartTimeBeforePresent = 0.25)
#   res = with( sg, approx( time, ne, rule = 1, xout = taxis )$y )
#   print( res ) 
#   res 
# }




# datetree same function as in variantAnalysis
date_trees <- function(mltr_fn, ofn, n_tree_dating = 10, civmd, meanrate, meanratesd, ncpu = 4, ...)
{
  mltr = read.tree(mltr_fn)
  
  # checking all samples have metadata attached... removing tips that aren't able to be matched
  mltr = lapply(mltr, function(tr) ape::drop.tip(tr, tr$tip.label[!tr$tip.label %in% civmd$central_sample_id]))
  
  
  tds = parallel::mclapply(mltr, function(tr) {
    tmp = lapply(1:n_tree_dating, function(x) {
      td = datetree( tr , civmd = civmd, meanrate =  max( 0.0001, rnorm( 1, meanrate, sd = meanratesd ) ) )
      td$tip.label =  paste0(td$tip.label, '|', as.Date(date_decimal(td$sts)), '|', td$sts )
      td
    })
    tmp
  }, mc.cores = ncpu)
  
  saveRDS( tds , file=paste0(ofn, "_dated_trees", '.rds' ))
  
  tds
}



n_tree_dating = 10
# meanrate = 0.0005
# matchedfn = "/cephfs/covid/bham/climb-covid19-volze/subsampling/matchSample_notB.1.1.7_2021-02-13.nwk"
# ncpu = 8


# make treedater trees for each ML tree.
tds_list = date_trees(mltr_fn = "/cephfs/covid/bham/climb-covid19-volze/subsampling/sampler1_B.1.1.7_2021-02-13_n=3000.nwk",
                      ofn = paste0('Sample_England_', 'sampler1_B.1.1.7_2021-02-13_n=3000', '_n_tree_dating_', n_tree_dating), civmd = civmd, 
                      meanrate = mr,n_tree_dating = 10,
                      meanratesd = mrsd, ncpu = 4)
# error likely due to not defining civmd$sample_time

###############################



