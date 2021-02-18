# NOTE THAT THIS IS d4 in PUBLIC RELEASE
# 
#  _   _                    ____   ____ _____ _____ 
# | \ | | ___  __   _____  / ___| / ___|_   _|  ___|
# |  \| |/ _ \ \ \ / / __| \___ \| |  _  | | | |_   
# | |\  |  __/  \ V /\__ \  ___) | |_| | | | |  _|  
# |_| \_|\___|   \_/ |___/ |____/ \____| |_| |_|    
# 

require(hrbrthemes)
require(scales)



# plot_mlesky function in d2_plot_mlesky.R

mlesky_B.1.1.7 = plot_mlesky(ofn ="C:/Users/lilyl/OneDrive/Documents/variantAnalysis/results/B.1.1.7_2021-02-13_n_tree_dating_10_dated_trees_mlesky.rds_mlesky.rds",
            ofn2 ="C:/Users/lilyl/OneDrive/Documents/variantAnalysis/results/notB.1.1.7_2021-02-13_n_tree_dating_10_dated_trees_mlesky.rds_mlesky.rds", 
            Lineage_main = "B.1.1.7",
            Lineage_matched = "Control",
            dedup = "", meanrate = "sampled" )





sgss_stp_new_43_52_weeks <- readRDS("~/N501Y_analyses/dropouts/sgss_stp_new_43_52_weeks.rds")
sgss_stp_new_43_56_weeks <- readRDS("~/N501Y_analyses/dropouts/sgss_stp_new_43_56_weeks.rds")


sgss_stp_new_43_52_weeks$true_negative = sgss_stp_new_43_52_weeks$sgss_s_negative * sgss_stp_new_43_52_weeks$total_cases /
  (sgss_stp_new_43_52_weeks$sgss_s_negative + sgss_stp_new_43_52_weeks$sgss_s_positive)


sgss_stp_new_43_52_weeks$true_negative_corrected = sgss_stp_new_43_52_weeks$sgss_s_negative_corrected * sgss_stp_new_43_52_weeks$total_cases /
  (sgss_stp_new_43_52_weeks$sgss_s_negative_corrected + sgss_stp_new_43_52_weeks$sgss_s_positive_corrected)



sgss_stp_new_43_56_weeks$true_negative = sgss_stp_new_43_56_weeks$sgss_s_negative * sgss_stp_new_43_56_weeks$total_cases /
  (sgss_stp_new_43_56_weeks$sgss_s_negative + sgss_stp_new_43_56_weeks$sgss_s_positive)


sgss_stp_new_43_56_weeks$true_negative_corrected = sgss_stp_new_43_56_weeks$sgss_s_negative_corrected * sgss_stp_new_43_56_weeks$total_cases /
  (sgss_stp_new_43_56_weeks$sgss_s_negative_corrected + sgss_stp_new_43_56_weeks$sgss_s_positive_corrected)

#  _ _ _ 
# | | | |
# | | | |
# |_|_|_|
# (_|_|_)
## !!!
#  _ _ _ 
# | | | |
# | | | |
# |_|_|_|
# (_|_|_)

sgss_stp_new_43_52_weeks = sgss_stp_new_43_56_weeks

mlesky_B.1.1.7$epiweek = lubridate::epiweek(lubridate::date_decimal(mlesky_B.1.1.7$time))

mlesky_B.1.1.7$epiweek = ifelse(mlesky_B.1.1.7$epiweek < 4, mlesky_B.1.1.7$epiweek+53, mlesky_B.1.1.7$epiweek)

#  _ _ _ 
# | | | |
# | | | |
# |_|_|_|
# (_|_|_)
## !!!
#  _ _ _ 
# | | | |
# | | | |
# |_|_|_|
# (_|_|_)
## !!!


pldf <- as.data.frame(do.call(rbind, lapply(unique(sgss_stp_new_43_52_weeks$epiweek), function(week) {
  if(nrow(sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$epiweek == week, ]) == 42) {# checking there are no duplicates; there are 42 STPs
    total_S_neg = sum(sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$epiweek == week, "sgss_s_negative_corrected"])
    
    ne = tail( mlesky_B.1.1.7$q_ne[which(mlesky_B.1.1.7$epiweek == week), ][,"y"], 1)
    neub = tail( mlesky_B.1.1.7$q_ne[which(mlesky_B.1.1.7$epiweek == week), ][,"yub"], 1)
    nelb = tail( mlesky_B.1.1.7$q_ne[which(mlesky_B.1.1.7$epiweek == week), ][,"ylb"], 1)
    
    
    return(c(week = week, total_S_neg = total_S_neg, ne = ne, neub = neub, nelb = nelb))
  }
}
)
)
)
pldf = pldf[!is.na(pldf$ne),]


saveRDS(pldf, file = "pldf_Ne_SGTF.rds")

require(ggplot2)
library(ggrepel)
require(hrbrthemes)

pl = ggplot(pldf, aes(x = ne, y = total_S_neg)) + geom_point( shape = 15) + 
  geom_errorbarh(aes(xmin = nelb, xmax = neub, y = total_S_neg)) +
  # scale_x_log10() + scale_y_log10() +
  scale_y_continuous(
    trans = "log10",
    breaks = function(x) {
      brks <- extended_breaks(Q = c(1, 5))(log10(x))
      10^(brks[brks %% 1 == 0])
    },
    labels = math_format(format = log10)
  ) +
  scale_x_continuous(
    trans = "log10",
    breaks = function(x) {
      brks <- extended_breaks(Q = c(1, 5))(log10(x))
      10^(brks[brks %% 1 == 0])
    },
    labels = math_format(format = log10)
  )+
  theme_bw() + labs(x = "Effective population size B.1.1.7", y = "TPP-adjusted SGTF case numbers") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14,face="bold"), panel.grid.minor = element_blank())+ annotation_logticks() +
  geom_label_repel(aes(x = ne, y = total_S_neg, label = week),alpha = 0.8)
pl

ggsave( plot = pl, file = "results/TPP-adjusted_SGTF_vs_Ne_mlesky_to_week_56.pdf", width = 8, height = 8 )


fit = lm(log(pldf$total_S_neg)~ log(pldf$ne))
coef(fit)

# my.formula <- log(pldf$total_S_neg) ~ (log(pldf$ne) - 1)
# 
# 
# my.formula <- log(pldf$total_S_neg) ~ (log(pldf$ne) - 1)

# nmod <- (lm(I(log(total_S_neg))~I(log(ne)-1) , pldf))

# pl + stat_smooth(formula = nmod,  method = "lm", alpha = 0.2)

# pl + geom_abline(intercept = exp(3.3979), slope = exp(1.4371) )



# pl + geom_abline(intercept = log(3.3979), slope = 1.437075 )
# 
# nmod <- (lm(I(log(total_S_neg))~I(log(ne)-1) , pldf))
# 
# 
# nmod = I(log(pldf$total_S_neg))~I(log(pldf$ne)-1)
# 
# 
# pl + geom_smooth(method = "lm", formula = nmod)

pl + geom_smooth(method = "lm", formula = y ~ x-10)



# nmod <- (lm(I(y-50)~I(x-10) +0, pldf))
nmod <- (lm(I(log(total_S_neg))~I(log(ne))  , pldf))
# predict(nmod, newdata = list(x=0))


abline(predict(nmod, newdata = list(ne=1)), coef(nmod), col='red')


### just to see what it looks like with raw S-

pldf <- as.data.frame(do.call(rbind, lapply(unique(sgss_stp_new_43_52_weeks$epiweek), function(week) {
  if(nrow(sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$epiweek == week, ]) == 42) {# checking there are no duplicates; there are 42 STPs
    total_S_neg = sum(sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$epiweek == week, "sgss_s_negative_corrected"])
    total_S_neg_raw = sum(sgss_stp_new_43_52_weeks[sgss_stp_new_43_52_weeks$epiweek == week, "sgss_s_negative"])
    
    ne = tail( mlesky_B.1.1.7$q_ne[which(mlesky_B.1.1.7$epiweek == week), ][,"y"], 1)
    neub = tail( mlesky_B.1.1.7$q_ne[which(mlesky_B.1.1.7$epiweek == week), ][,"yub"], 1)
    nelb = tail( mlesky_B.1.1.7$q_ne[which(mlesky_B.1.1.7$epiweek == week), ][,"ylb"], 1)
    
    
    return(c(week = week, TPR_adjusted = total_S_neg, Raw = total_S_neg_raw, ne = ne, neub = neub, nelb = nelb))
  }
}
)
)
)
names(pldf) = c("Week", "TPP-adjusted", "Raw", "ne", "neub" ,     "nelb")
pldf = reshape2::melt(pldf, id.vars = names(pldf)[!names(pldf) %in% c("Raw" , "TPP-adjusted"  )])


pldf = pldf[!is.na(pldf$ne),]

tmp_unadj = pldf[pldf$variable == "Raw", ]
tmp_adj = pldf[pldf$variable == "TPP-adjusted", ]

my.formula_unadj <- log(tmp_unadj$value) ~ log(tmp_unadj$ne)
my.formula_adj <- log(tmp_adj$value) ~ log(tmp_adj$ne)




pl_2 = ggplot(pldf,aes(x = ne, y = value, col = variable, fill = variable)) + geom_point( shape = 15, size = 2)  + 
  geom_errorbarh(aes(xmin = nelb, xmax = neub, y = value, col = variable)) + 
  # scale_x_log10() + 
  scale_y_continuous(
    trans = "log10",
    breaks = function(x) {
      brks <- extended_breaks(Q = c(1, 5))(log10(x))
      10^(brks[brks %% 1 == 0])
    },
    labels = math_format(format = log10)
  ) + 
  scale_x_continuous(
    trans = "log10",
    breaks = function(x) {
      brks <- extended_breaks(Q = c(1, 5))(log10(x))
      10^(brks[brks %% 1 == 0])
    },
    labels = math_format(format = log10)
  )+    
  stat_smooth(method = "lm", alpha = 0.2)+
  
  theme_bw() + labs(x = "Effective population size B.1.1.7", y = "SGTF case numbers", col = "SGTF", fill = "SGTF") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14,face="bold"), legend.position = c(0.7,0.3), legend.title = element_text(size=14), 
        legend.text =element_text(size=12) ,legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        plot.margin = margin(10, 10, 10, 10), panel.grid.minor = element_blank())+ 
  annotation_logticks() +
  ggpmisc::stat_poly_eq(data = tmp_unadj, formula = my.formula_unadj, size = 5,
                        mapping = aes(log(ne), log(value),  col = variable, fill = variable,
                                      label = paste( ..adj.rr.label.., sep = "~~~")),  parse = TRUE)+  
  ggpmisc::stat_poly_eq(data = tmp_adj, label.y = 0.85, size = 5,
                        formula = my.formula_adj, mapping = aes(log(ne), log(value), col = variable, fill = variable,
                                                                label = paste( ..adj.rr.label.., sep = "~~~")),  parse = TRUE)

pl_2


ggsave( plot = pl_2, file = "results/TPP-adjusted_SGTF_vs_Ne_mlesky_raw_vs_adjusted_to_week_56.pdf", width = 8, height = 8 )


