# ======================================================================================================= #
# R script for reproduction of the analyses in Piatscheck et al 2018 Acta Oecologica
# ======================================================================================================= #

# Paper title: Ecological factors associated with pre-dispersal predation of fig seeds and wasps by fig-specialist 
# lepidopteran larvae
# doi.org/10.1016/j.actao.2018.03.001

# The following script allows the reproduction of figure 3 and 4 and the test of association between pre-dispersal 
# fig predation by Omiodes stigmosalis on Ficus petiolaris and ecological variables using generalized linear mixed 
# models.

# The phylogenetic analysis and spatial aggregation of fig predation analysis (as well as reproduction of Figures 1,
# 2 and 5) are not presented here as they were not realized in R.

# ------------------------------------------------------------------------------ #
# Import the data and formatting
# ------------------------------------------------------------------------------ #

lep_data <- read.csv("lep_data_Piatscheck.et.al.2018.csv")

lep_data$lep_resp_mat <- cbind(lep_data$damaged_syconia, lep_data$counted_syconia-lep_data$damaged_syconia)
lep_data$site <- factor(lep_data$site, levels = c("158", "172","112","113","95","179","201","96","70"))
lep_data$season <- factor(lep_data$season, levels = c("F2012", "S2013","F2013","S2014")) 
lep_data$treeID <- as.factor(paste(lep_data$site, lep_data$tree, sep = "-"))

# ------------------------------------------------------------------------------ #
# Damage proportions plots
# ------------------------------------------------------------------------------ #

library(ggplot2)
lep_data$damage_prop <- lep_data$damaged_syconia/lep_data$counted_syconia

# ______________________________________________________________________________ #
# Figure 3: Damage proportion density plot
ggplot(lep_data, aes(damage_prop)) + 
  geom_density() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        text = element_text(size=22)) + 
  labs(x = "Damage proportions", y = "Density")

# ______________________________________________________________________________ #
# Figure 4: Damage proportion by sites and seasons plot
ggplot(lep_data, aes(site, damage_prop)) + 
  geom_boxplot(aes(fill=season)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        text = element_text(size=20),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.9, 0.85),
        legend.background = element_rect(color = "black", 
                                         fill = "grey90", 
                                         size = 1, 
                                         linetype = "solid")) +
  labs(x = "Sites", y = "Larval Damage Proportion")  +
  scale_fill_grey()

# ------------------------------------------------------------------------------ #
# Ecological correlates of lepidopteran damage
# ------------------------------------------------------------------------------ #
  
library(lme4)
lep_data$interphases_prop <- (lep_data$early_interph + lep_data$late_interph)/lep_data$total_syconia

# Because the glmer function works with variables approximatively the same scale, some varibles need to be rescaled.
lep_data$volume_tr  <- lep_data$volume/1000
lep_data$est_n_fruits_tr <- lep_data$est_n_fruits/1000

# ______________________________________________________________________________ #
# Generalized Linear Mixed Models (logistic regression)

# Definitive results
summary(glmer(lep_resp_mat ~ site + 
                season + 
                distance_to_NN_20km +
                volume_tr + 
                reproduction_prop + 
                interphases_prop + 
                (1|treeID) + 
                (1|season:treeID), 
              data=lep_data, 
              family=binomial, 
              control=glmerControl(optimizer="nlminbwrap")))
# Note: the positive association (marginally significant) of damaged syconia relative to non-damaged with tree density
# needs to be interpreted as negative since we used neighbor joining. 

# GLMM with "asynchrony" instead of "proportion of interphase"
summary(glmer(lep_resp_mat ~ site + 
                season + 
                distance_to_NN_20km + 
                est_n_fruits_tr + 
                interphases_prop + 
                (1|treeID) + 
                (1|season:treeID), 
              data=lep_data, 
              family=binomial, 
              control=glmerControl(optimizer="nlminbwrap")))

# with "crop size" instead of "volume" and "reproductive effort"
summary(glmer(lep_resp_mat ~ site + 
                season + 
                distance_to_NN_20km + 
                est_n_fruits_tr + 
                synchrony + 
                (1|treeID) + 
                (1|season:treeID), 
              data=lep_data, 
              family=binomial, 
              control=glmerControl(optimizer="nlminbwrap")))

