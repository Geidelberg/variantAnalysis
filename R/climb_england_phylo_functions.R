

#' Select a random sample from aligment stratified through time 
#'
#' Select samples based on quantile of sample time distribution. Requires date to be at and of sequence label
#'
#" @param path_to_align A DNAbin alignment *or* a system path (type character) where original alignment can be found, such as /gisaid/gisaid_cov2020_sequences_March14_aligned.fas
#' @param path_to_save Where to store (as fasta) the filtered alignment
#' @param q_threshold Clock outlier threshold 
#' @param minEdge minimum branch length (substitutions per site) to stabilize clock inference 
#' 
#' @return A DNAbin alignment. Will also save to path_to_save 
#' @export 
time_stratified_sample_cog <- function(n, nreps=1, mindate, maxdate, md, d
) {
  
  
  library( ape ) 
  library( treedater )
  library( lubridate )
  
  
  sts <-    decimal_date( ymd( md$sample_date))
  
  
  names( sts ) <-  md$central_sample_id 
  
  ssts = sort( sts )
  N <- length( sts )
  
  keep_list = list()
  for(rep in 1:nreps) {
    ibins = floor( seq( N/n , N , length = n ) )
    i0 = 1 
    keep <- c() 
    for ( k in 1:n){
      keep <- c( keep, sample( names(ssts)[i0:ibins[k]], size = 1 )  )
      i0 <- ibins[ k ] + 1
    }
    keep <- unique( keep )
    
    keep_list[[rep]] = keep
    
    d1 = d[rownames(d) %in% keep,]
    
    
    outfn <- paste0('Sample_England_', Sys.Date(), '_n_', n, '_rep_', rep, '.fasta')
    write.dna( d1, file = outfn, format = 'fasta' )
  }
  
  
  list( 
    # alignment = d1 
    # , sids = keep 
    # , sts= sts[keep ]
    # , outfn = outfn
    sids = keep_list
    
  ) 
}



# tds_list = date_trees(mltr_fn = "/cephfs/covid/bham/climb-covid19-geidelbergl/lockdown/2021-01-28/Sample_England_2021-01-28_n_30_nreps_3.nwk",
#            ofn = paste0('Sample_England_', Sys.Date(), '_n_', n, '_nreps_', nreps, '_test'), civmd = md, meanrate = meanrate)



# datetree same function as in variantAnalysis
date_trees <- function(mltr_fn, ofn, n_tree_dating = 2, civmd, meanrate, ncpu = 4, ...)
{
  mltr = read.tree(mltr_fn)
  
  tds = parallel::mclapply(mltr, function(tr) {
    tmp = lapply(1:n_tree_dating, function(x) {
      td = datetree( tr , civmd = civmd, meanrate = meanrate )
      td$tip.label =  paste0(td$tip.label, '|', as.Date(date_decimal(td$sts)), '|', td$sts )
      td
    })
    tmp
  }, mc.cores = ncpu)
  
  saveRDS( tds , file=paste0(ofn, "_dated_trees", '.rds' ))
  
  tds
}


#######
# tds_list is a list of lists. First level = number of alignments/ML trees (nreps). Second level = number of times each ML is dated with treedater (n_tree_dating).
# res_mlesky = run_mlesky(tds_list = "/cephfs/covid/bham/climb-covid19-geidelbergl/lockdown/2021-02-02/Sample_England_2021-02-02_n_2000_nreps_20_test_dated_trees.rds",
#                         ofn = paste0('Sample_England_', Sys.Date(), '_n_', n, '_nreps_', nreps, '_test'), taxis = taxis)

run_mlesky <- function(tds_list, ofn, taxis = taxis, NeStartTimeBeforePresent = 0.25) {
  
  if(!inherits(tds_list, what = c("list")))
    tds_list = readRDS(tds_list)
  
  
  res_mlesky_list = lapply(tds_list, function(tds) {
    
    times = diff( range( epiweek( date_decimal( tds[[1]]$sts )  )  )  ) + 1 
    res <- times * 2
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

############
run_skygrowth <- function(tds_list, ofn, taxis = taxis, covar_data = NULL, covar_formula = ~ transit_stations_percent_change_from_baseline, ncpu=6) {
  
  if(!inherits(tds_list, what = c("list")))
    tds_list = readRDS(tds_list)
  
  if(is.null(covar_data)) {
    covar_info = "no_covar"
  }  else  { 
    covar_info = paste0("covar_", as.character(covar_formula)[2])
    covar_data = covar_data[, c("time", as.character(covar_formula)[2])]
  }
  
  res_skygrowth_list = lapply(tds_list, function(tds) {
    
    times = diff( range( epiweek( date_decimal( tds[[1]]$sts )  )  )  ) + 1 
    res <- times * 2
    class( tds ) <- 'multiPhylo' 
    
    # tds = list(td)
    if(is.null(covar_data)){
      out = skygrowth0_cog(tds =  tds, tstart = min(tds[[1]]$sts), logRmean = log(0.70), logRsd = .5, res = res, ncpu = ncpu)
      out = out[c("taxis", "date", "Ne", "growth", "R", "sgs")]
      
    } else    {
      out = skygrowth1_cog( tds = tds, X=covar_data, formula = covar_formula, tstart = min(tds[[1]]$sts), logRmean = log(0.70), logRsd = .5, res = res, ncpu = ncpu) 
      out = out[c("taxis", "date", "Ne", "growth", "R", "beta", "data")]
      
    }
    
    out
    
    
  })
  
  saveRDS( res_skygrowth_list, file=paste0(ofn, "_skygrowth", "_covar_", covar_info, '.rds' ))
  
  list(res_skygrowth = res_skygrowth_list, covar_info = covar_info)
  
}
# res_skygrowth = run_skygrowth(tds_list = "/cephfs/covid/bham/climb-covid19-geidelbergl/lockdown/2021-02-02/Sample_England_2021-02-02_n_2000_nreps_20_test_dated_trees.rds",
#                         ofn = paste0('Sample_England_', Sys.Date(), '_n_', n, '_nreps_', nreps, '_test'), taxis = taxis, 
#                         covar_data = google_mob,
#                         covar_formula = ~ transit_stations_percent_change_from_baseline, ncpu = 8)


