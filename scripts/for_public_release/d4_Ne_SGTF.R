# Compares the Ne of B.1.1.7 estimated from mlesky against TPP-adjusted SGTF
# Required: 
# - RDS files of mlesky outputs from d2_mlesky.R

require( hrbrthemes )
require( scales )
require( ggplot2 )
require( lubridate )
require( cowplot )
require( grid )
require( gridExtra )
require( epitrix )
library( ggrepel )



# Compares B.1.1.7 to control; function from d3_plot_mlesky.R
# This functino outputs Ne(t) of B.1.1.7 which we use in this script
Ne_t_B.1.1.7 = plot_mlesky(ofn ="Sample_England_sampler1_B.1.1.7_2021-02-13_n=3000_n_tree_dating_10_mlesky.rds",
                           ofn2 ="Sample_England_matchSample_control_2021-02-13_n_tree_dating_10_mlesky.rds", 
                           lineage_main = "B.1.1.7",
                           lineage_matched = "Control" , si = si)

# Converts decimal date to epiweek; the first few weeks of 2021 are interpreted as weeks 54-56 of 2020
Ne_t_B.1.1.7$epiweek = lubridate::epiweek(lubridate::date_decimal(Ne_t_B.1.1.7$time))
Ne_t_B.1.1.7$epiweek = ifelse(Ne_t_B.1.1.7$epiweek < 4, Ne_t_B.1.1.7$epiweek+53, Ne_t_B.1.1.7$epiweek)

# Read SGSS data
sgss_stp_new_43_56_weeks <- readRDS("data/sgss_stp_new_43_56_weeks.rds")


# Creates data frame extracting the total number of TPP-adjusted SGTF across all STPs for each epiweek; Ne (and HPD) at the end of each epi week estimated from mlesky.
pldf <- as.data.frame(do.call(rbind, lapply(unique(sgss_stp_new_43_56_weeks$epiweek), function(week) {
  if(nrow(sgss_stp_new_43_56_weeks[sgss_stp_new_43_56_weeks$epiweek == week, ]) == 42) {# checking there are no duplicates; there are 42 STPs
    total_S_neg = sum(sgss_stp_new_43_56_weeks[sgss_stp_new_43_56_weeks$epiweek == week, "sgss_s_negative_corrected"])
    ne = tail( Ne_t_B.1.1.7$q_ne[which(Ne_t_B.1.1.7$epiweek == week), ][,"y"], 1)
    neub = tail( Ne_t_B.1.1.7$q_ne[which(Ne_t_B.1.1.7$epiweek == week), ][,"yub"], 1)
    nelb = tail( Ne_t_B.1.1.7$q_ne[which(Ne_t_B.1.1.7$epiweek == week), ][,"ylb"], 1)
    return(c(week = week, total_S_neg = total_S_neg, ne = ne, neub = neub, nelb = nelb))
  }
}
)
)
)

# Removes any rows with NAs
pldf = pldf[!is.na(pldf$ne),] 

# Plot
pl = ggplot(pldf, aes(x = ne, y = total_S_neg)) + geom_point( shape = 15) + 
  geom_errorbarh(aes(xmin = nelb, xmax = neub, y = total_S_neg)) +
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

ggsave( plot = pl, file = "TPP-adjusted_SGTF_vs_Ne_mlesky_to_week_56.pdf", width = 8, height = 8 )

