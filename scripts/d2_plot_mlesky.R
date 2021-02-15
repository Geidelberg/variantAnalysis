
require(ggplot2)
require(lubridate)
library(cowplot)
library(grid)
library(gridExtra)


plot_mlesky <- function(ofn, ofn2, Lineage_main, Lineage_matched, dedup, meanrate) {
  tN = readRDS( ofn )
  tN2 = readRDS( ofn2 )
  
  # attach( tN )
  q_ne = t(apply( tN$ne, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ))
  q_mane = t(apply( tN2$ne, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ))
  
  gr =  apply( tN$ne, 2, function(x) c(exp( diff(log(x)) )^6.5, NA)  ) 
  magr =  apply( tN2$ne, 2, function(x) c(exp( diff(log(x)) )^6.5, NA)  ) 
  
  q_gr = t( apply( gr, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ) )
  q_magr = t( apply( magr, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ) )
  
  colnames( q_ne ) = colnames( q_mane ) = c( 'y', 'ylb', 'yub' )
  pldf0 = as.data.frame( q_ne ) ; pldf0$Lineage = Lineage_main; pldf0$time = tN$time
  pldf1 = as.data.frame( q_mane ); pldf1$Lineage = Lineage_matched; pldf1$time = tN2$time
  pldf = rbind( pldf0, pldf1 )
  
  p0 = ggplot( aes(x = as.Date( date_decimal( time)), y = y, colour = Lineage, fill = Lineage , ymin = ylb, ymax = yub ) , data = pldf ) +
    geom_path(size=1) + geom_ribbon( alpha = .25, col = NA ) + xlab('') + ylab('Effective population size' ) +
    theme_minimal() + theme(legend.position='none',panel.grid.minor = element_blank())+
    scale_x_date(date_breaks = "1 month", date_labels = '%b')+theme(axis.text=element_text(size=12),
                                                                    axis.title=element_text(size=14)) + scale_y_log10()+
    annotation_logticks(colour = 'grey', short = unit(.05, "cm"), mid = unit(.05, "cm"), long = unit(.05, "cm"))
  
  colnames( q_gr ) = colnames( q_magr ) = c( 'y', 'ylb', 'yub' )
  gpldf0 = as.data.frame( q_gr ) ; gpldf0$Lineage = Lineage_main; gpldf0$time = tN$time
  gpldf1 = as.data.frame( q_magr ); gpldf1$Lineage = Lineage_matched; gpldf1$time = tN2$time
  gpldf = rbind( gpldf0, gpldf1 )
  
  Rratio = gr / magr
  qRr =  t( apply( Rratio, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ) )
  colnames( qRr ) =  c( 'y', 'ylb', 'yub' )
  Rrpldf = as.data.frame( qRr ); Rrpldf$time = tN$time
  
  
  r_range = range(c(gpldf$yub, Rrpldf$yub, gpldf$ylb, Rrpldf$ylb), na.rm = T)
  
  p1 =  ggplot( aes(x = as.Date( date_decimal( time)), y = y, colour = Lineage, fill = Lineage , ymin = ylb, ymax = yub ) , data = gpldf ) +
    geom_path(size=1) + geom_ribbon( alpha = .25 , col = NA) + xlab('') + ylab('Reproduction number' ) + theme_minimal() +
    scale_y_log10(limits = r_range) + annotation_logticks(colour = 'grey', short = unit(.05, "cm"), mid = unit(.05, "cm"), long = unit(.05, "cm"))   +
    geom_hline( aes( yintercept=1), lty = 2 )  + theme(legend.position='', panel.grid.minor = element_blank())+
    scale_x_date(date_breaks = "1 month", date_labels = '%b')+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  
  p2 = ggplot( aes(x = as.Date( date_decimal( time)), y = y, ymin = ylb, ymax = yub ) , data = Rrpldf ) + geom_path(size=1) +
    
    geom_ribbon( alpha = .25, col = NA ) + xlab('') + ylab('Ratio of reproduction numbers' ) + theme_minimal() +
    scale_y_log10(limits = r_range) + annotation_logticks(colour = 'grey', short = unit(.05, "cm"), mid = unit(.05, "cm"), long = unit(.05, "cm")) + theme( panel.grid.minor = element_blank())+
    geom_hline( aes( yintercept=1), lty = 2 ) +scale_x_date(date_breaks = "1 month", date_labels = '%b')+
    theme(axis.text=element_text(size=12),  axis.title=element_text(size=14))
  
  legend <- cowplot::get_legend(
    p1 + theme(legend.box.margin = margin(0, 0, 0, 12), legend.position = "top", legend.title = element_text(size=14), legend.text =element_text(size=12) )
  )
  
  
  
  P0 = cowplot::plot_grid( plotlist = list( p0, p1, p2 ), nrow = 1, align = "v" )
  
  P0 = cowplot::plot_grid(legend, P0, ncol = 1, rel_heights =  c(0.1, 1))
  
  ggsave( plot = P0, file = paste0('results/d1_', Lineage_main, '_', Lineage_matched, dedup, '_', "meanrate", meanrate, '.pdf'), width = 12, height = 4.5 )
  
  print(P0)
  
  Time=time
  
  # detach(tN)
  return(list(ofn = ofn, 
              time = Time,
              pldf0 = pldf0,
              gpldf0 = gpldf0
  ))
}



# plot_mlesky(ofn ="C:/Users/lilyl/OneDrive/Documents/variantAnalysis/results/notB.1.1.7_2021-02-13_n_tree_dating_10_dated_trees_mlesky.rds_mlesky.rds", dedup = "", meanrate = 0.0005 )
# plot_mlesky(ofn ="C:/Users/lilyl/OneDrive/Documents/variantAnalysis/results/B.1.177_2021-02-13_n_tree_dating_10_dated_trees_mlesky.rds_mlesky.rds", dedup = "", meanrate = 0.0005 )
plot_mlesky(ofn ="C:/Users/lilyl/OneDrive/Documents/variantAnalysis/results/B.1.177_2021-02-13_n_tree_dating_10_dated_trees_mlesky.rds_mlesky.rds",
            ofn2 ="C:/Users/lilyl/OneDrive/Documents/variantAnalysis/results/notB.1.1.7_2021-02-13_n_tree_dating_10_dated_trees_mlesky.rds_mlesky.rds", 
            Lineage_main = "B.1.177",
            Lineage_matched = "Control",
            dedup = "", meanrate = 0.0005 )
