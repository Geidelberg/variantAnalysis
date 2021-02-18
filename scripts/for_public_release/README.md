# Estimating effective population size (Ne) through time of lineage B.1.1.7 using mlesky

Steps to perform [phylodynamic inference](https://www.biorxiv.org/content/10.1101/2021.01.18.427056v1) of effective population size for time-scaled phylogenies. 

These scripts require installation of the `mlesky` package, which can be [downloaded here](https://github.com/emvolz-phylodynamics/mlesky). To download the package directly in R, you can run the following:


```r
devtools::install_github('emvolz-phylodynamics/mlesky') 
```

We will also need the `treedater` package [as described here](https://academic.oup.com/ve/article/3/2/vex025/4100592). To install run the following:

```r
install.packages("treedater")
```



Other packages required:


```r
library( ape )
library( lubridate )
library( glue )
library( sarscov2 ) 
library( ggplot2 )
library( grid )
library( gridExtra )
library( ggtree )
library( alakazam )
library( stringi )
```





## Lineage B.1.1.7

To perform this analysis, you will need 
1) a list of maximum-likelihood trees in the `.nwk` format built from sequences from the B.1.1.7 lineage (and any other lineage you wish to analyse)
2) a `.csv` metadata file with a column of sequence names (that match `tip.label` of the trees in your `.nwk` file) and their sample dates


Loading tree:
```r
mltr <- ape::read.tree("sampler1_B.1.1.7_2021-02-13_n=3000.nwk")
mltr

```


```
## 100 phylogenetic trees
```

Loading metadata:
```r
metadata = read.csv( "metadata.csv" , stringsAs = FALSE , header=TRUE )
metadata$sample_date <- as.Date( metadata$sample_date )
metadata$sample_time <- decimal_date( metadata$sample_date ) # converting to decimal date for use in treedater
```

## Dating the trees


To date the trees we have to assume a clock rate. To incorporate uncertainty around this parameter in our results we sample this from a distribution. Here we define this as a normal distribution with mean and standard deviation derived from the posterior from model-based phylodynamic analyses.

Mean clock rate and 95% HPD interval: 5.9158 x 10<sup>-4</sup>	(5.2065 x 10<sup>-4</sup> - 6.7144 x 10<sup>-4</sup>)
```r
mr = 5.9158E-4
mrci = 	c( 5.2065E-4, 6.7144E-4)
mrsd = diff( mrci ) / 4 / 1.96

```

We are using `treedater` to date every ML tree ten times; each time will assume a clock rate drawn from the above distribution. The code to perform this is below:


```r
# make treedater trees for each lineage
tds_list_B.1.1.7 = date_trees(
  mltr_fn = "sampler1_B.1.1.7_2021-02-13_n=3000.nwk",
  ofn = paste0('Sample_England_sampler1_B.1.1.7_2021-02-13_n=3000_n_tree_dating_10'), 
  metadata = metadata, 
  meanrate = mr,
  n_tree_dating = 10,
  meanratesd = mrsd, 
  ncpu = 4)

```
This will save your dated trees as an `.RDS` file. The full code (and functions) to do this are provided in `d1_date_trees.R`. Also given is the code to perform the same methods on trees from lineage B.1.177 and also control sequences from across England matched by time and place to the B.1.1.7 alignments.



## Performing mlesky on the dated trees


From `d1_date_trees.R`, we have produced a list of dated trees. Next we define a time period over which we would like to estimate Ne(t), and run mlesky:


```r
# time period over which to estimate Ne
taxis = decimal_date( seq( as.Date( '2020-10-15') , as.Date('2021-01-24'), by = 1) )

# run mlesky
res_mlesky_B.1.1.7 = run_mlesky(
  tds_list = "Sample_England_sampler1_B.1.1.7_2021-02-13_n=3000_n_tree_dating_10_dated_trees.rds",
  ofn = "Sample_England_sampler1_B.1.1.7_2021-02-13_n=3000_n_tree_dating_10", 
  taxis = taxis)
```
This will save your mlesky outputs as an `.RDS` file. The functions are provided in `d2_mlesky.R`, as well as the code to perform the same analysis for other lineages.



## Comparing Ne(t) with TPP-adjusted SGTF
