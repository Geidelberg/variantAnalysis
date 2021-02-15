
#######
# tds_list is a list of lists. First level = number of alignments/ML trees (nreps). Second level = number of times each ML is dated with treedater (n_tree_dating).
# res_mlesky = run_mlesky(tds_list = "/cephfs/covid/bham/climb-covid19-geidelbergl/lockdown/2021-02-02/Sample_England_2021-02-02_n_2000_nreps_20_test_dated_trees.rds",
#                         ofn = paste0('Sample_England_', Sys.Date(), '_n_', n, '_nreps_', nreps, '_test'), taxis = taxis)

run_mlesky <- function(tds_list, ofn, taxis = taxis, NeStartTimeBeforePresent = 0.25) {
  
  if(!inherits(tds_list, what = c("list")))
    tds_list = readRDS(tds_list)
  
  
  res_mlesky_list = lapply(tds_list, function(tds) {
    
    # times = diff( range( epiweek( date_decimal( tds[[1]]$sts )  )  )  ) + 1 
    weeks = round(as.numeric((date_decimal(max(tds[[1]]$sts))-date_decimal(min(tds[[1]]$sts)))/7))
    res <- weeks * 2
    class( tds ) <- 'multiPhylo' 
    
    tds = lapply( tds , 
                  function(x) {x$tip.label = unlist(lapply(strsplit(x$tip.label, '[|]'), function(y) paste0(y[1])))
                  return(x)}
    )
    
    sgs = parallel::mclapply( tds, function(td) {
      mlskygrid(td, tau = NULL, tau_lower=.001, tau_upper = 10 , sampleTimes = td$sts , res = res, ncpu = 6, NeStartTimeBeforePresent = NeStartTimeBeforePresent)
    }, mc.cores = 5 )
    
    
    out = lapply(sgs, function(sg) {
      with( sg, approx( time, ne, rule = 1, xout = taxis )$y )
    })
    
    out
    
  })
  
  # I am collapsing the results from all alignments and all dated trees together as one.
  res_mlesky <- 
    do.call( cbind, lapply(res_mlesky_list, function(x) do.call( cbind, x ) ))
  
  saveRDS( list( time = taxis, ne = res_mlesky ) , file=paste0(ofn, "_mlesky", '.rds' ))
  
  res_mlesky
}
