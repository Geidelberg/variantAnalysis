


# library( sarscov2 ) 
library( ape ) 
library( lubridate )
library( treedater ) 
require(ggplot2)
require(grid)
require(gridExtra)
require(ggtree)
library(alakazam)
require(stringi)
require( skygrowth )
require( mlesky )

source("/cephfs/covid/bham/climb-covid19-geidelbergl/R/erik_skygrowth_functions_cog.R")
source("/cephfs/covid/bham/climb-covid19-geidelbergl/R/functions.R")
source("/cephfs/covid/bham/climb-covid19-geidelbergl/R/sampling.R")
source("/cephfs/covid/bham/climb-covid19-geidelbergl/R/dating_and_mlesky_functions.R")
source("/cephfs/covid/bham/climb-covid19-geidelbergl/R/climb_england_phylo_functions.R")




# metadata 
civetfn =  list.files(  '/cephfs/covid/bham/results/msa/20210604/alignments/' , patt = 'cog_[0-9\\-]+_metadata.csv', full.names=TRUE) #'../phylolatest/civet/cog_global_2020-12-01_metadata.csv'
civmd = read.csv( civetfn , stringsAs=FALSE , header=TRUE )
civmd$central_sample_id=civmd$sequence_name

# civmd$central_sample_id <-  sapply( strsplit( civmd$sequence_name , split='/' ) , '[', 2 ) # for linkage 
civmd$sample_date <- as.Date( civmd$sample_date )
civmd$sample_time <- decimal_date( civmd$sample_date ) 


mltr_fn = 'B.1.617.2'
mltr_list = list.files(  '/cephfs/covid/bham/climb-covid19-geidelbergl/617.2/f0-trees' , patt = mltr_fn, full.names=TRUE)
mltr = lapply(mltr_list, read.tree)

sts <- lapply(mltr, function(x) {
  civmd$sample_time[  match( x$tip.label , civmd$central_sample_id ) ]
})

# taxis = decimal_date( seq( as.Date( '2020-10-15') , as.Date('2021-01-24'), by = 1) )
taxis = decimal_date( seq( as.Date(date_decimal(min(unlist(sts)))) , as.Date(date_decimal(max(unlist(sts)))), by = 1) )

res_mlesky = run_mlesky(tds_list = "/cephfs/covid/bham/climb-covid19-geidelbergl/617.2/Sample_England_B.1.617.2_n_tree_dating_5_dated_trees.rds",
                        ofn = "/cephfs/covid/bham/climb-covid19-geidelbergl/617.2/Sample_England_B.1.617.2_n_tree_dating_5_dated_trees_mlesky.rds", taxis = taxis)

