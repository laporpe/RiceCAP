############################################################################
#
#                 visualize rice GP
#
############################################################################
library(dplyr)
library(ggplot2)

gp_results_LA <- read.table("/Users/maryfrancis/Documents/Diepenbrock_Lab/rice-cgm-wgp/figures/GP_explore_LA_alltraits.txt", header = TRUE)
gp_results_MS <- read.table("/Users/maryfrancis/Documents/Diepenbrock_Lab/rice-cgm-wgp/figures/GP_explore_MS_alltraits.txt", header = TRUE)
gp_results_AR <- read.table("/Users/maryfrancis/Documents/Diepenbrock_Lab/rice-cgm-wgp/figures/GP_explore_AR_alltraits.txt", header = TRUE)

gp_results_LA["enviro"] <- "LA"
gp_results_MS["enviro"] <- "MS"
gp_results_AR["enviro"] <- "AR"

gp_results <- rbind(gp_results_LA,gp_results_MS,gp_results_AR)


traits <- unique(gp_results$Trait)

day_traits <- c("Avg_MATDATE","Avg_HeadDate","Avg_EmergeDate","Avg_DaysEMERG","Avg_DaysHEAD","Avg_DaysMATURITY","Avg_DaysGrainfill")
  
non_day_traits <- c("Avg_SEEDNUMperPANMEAN","Avg_PANLENGTHPLANTMEAN","Avg_TILLNUMMEAN","Avg_MeanHT","Avg_SEEDSAMPLEwtPANmean","Avg_TOTSEEDWTmean","Avg_TOTBIOMASS")


gp_results %>% ggplot(aes(x=enviro, y=Accuracy, fill = enviro)) + facet_wrap("Trait") + geom_violin() + scale_fill_brewer(palette="YlGnBu") +
  theme_bw()  + ylab("GP accuracy") + xlab("Trait") + scale_x_discrete(guide = guide_axis(angle = 45)) + ggtitle("All Environments")

data_summary <- function(x) {
  m <- mean(x)
  ymin <- unname(quantile(x, probs = 0.25))
  ymax <- unname(quantile(x, probs = 0.75))
  return(c(y=m,ymin=ymin,ymax=ymax))
}
gp_results <- gp_results %>% filter(Trait == c("Avg_SEEDNUMperPANMEAN","Avg_PANLENGTHPLANTMEAN","Avg_TILLNUMMEAN","Avg_MeanHT","Avg_DaysHEAD","Avg_DaysMATURITY",
                      "Avg_DaysGrainfill","Avg_SEEDSAMPLEwtPANmean","Avg_TOTSEEDWTmean","Avg_TOTBIOMASS"))


gp_results$Trait <- factor(gp_results$Trait, levels = c("Avg_SEEDNUMperPANMEAN","Avg_PANLENGTHPLANTMEAN","Avg_TILLNUMMEAN","Avg_MeanHT","Avg_DaysHEAD","Avg_DaysMATURITY",
                                                        "Avg_DaysGrainfill","Avg_SEEDSAMPLEwtPANmean","Avg_TOTSEEDWTmean","Avg_TOTBIOMASS"), 
                           labels = c("Seed Number per Panicle","Panicle Length","Tiller Number","Plant Height","Days to Heading","Days to Maturity","Days to Grainfill",
                                      "Seed Weight per Panicle","Total Seed Weight","Total Biomass"))


gp_results %>% group_by(enviro, Trait, Iteration)  %>%
  dplyr::summarize(median = median(Accuracy, na.rm=TRUE)) %>% ggplot(aes(x=enviro, y=median, fill = enviro)) +
  geom_violin() + facet_wrap("Trait")  + stat_summary(fun.data=data_summary,color = "black", shape = 20, position = position_dodge(0.9))  +
  theme_bw(base_size = 14)  + ylab("predictive ability") + xlab("environment") + 
  labs(fill = "GP method") + theme(strip.text.x = element_text(size = 14))





