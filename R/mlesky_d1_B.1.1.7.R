


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



taxis = decimal_date( seq( as.Date( '2020-10-15') , as.Date('2021-01-24'), by = 1) )


res_mlesky = run_mlesky(tds_list = "/cephfs/covid/bham/climb-covid19-geidelbergl/variantAnalysis/Sample_England_sampler1_B.1.1.7_2021-02-13_n=3000_n_tree_dating_10_dated_trees.rds",
                        ofn = "/cephfs/covid/bham/climb-covid19-geidelbergl/variantAnalysis/B.1.1.7_2021-02-13_n_tree_dating_10_dated_trees_mlesky.rds", taxis = taxis)

