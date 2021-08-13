# ---------------------------------------------------------------------------------------------------------------------------#
#   Sensitivity Analysis - checking sensitivity to breaks used in ordinal regression
# ---------------------------------------------------------------------------------------------------------------------------#
library(tidyverse)
library(dplyr)
library(brms)
library(RColorBrewer) 
library(ggmcmc)
library(ggthemes)
library(bayesplot)
library(tidybayes)
library(ggridges)
library(bayestestR)
library(ggplot2)
library(sf)
library(patchwork)
mytime <- format(Sys.time(), "%Y-%m-%d")

dat2000 <- read.csv("Data/Data_seagrass_trends_pressures_2000-2010.csv")

# create different breaks to apply to models
categories <- list(c(-Inf, -0.025, -0.005, 0.005, 0.025, Inf),
              c(-Inf, -0.05, -0.005, 0.005, 0.05, Inf),
              c(-Inf, -0.025, -0.01, 0.01, 0.025, Inf),
              c(-Inf, -0.05, -0.01, 0.01, 0.05, Inf),
              c(-Inf, -0.1, -0.02, 0.02, 0.1, Inf))

model_list <- vector(mode = "list", length = length(categories))    


for(i in 1:length(categories)){
  
  dat2000$category <- cut(dat2000$IRtrend_X50., 
                          breaks= categories[[i]], 
                          labels=c("rapid_decline","slow_decline","stable", "slow_increase", "rapid_increase"),ordered_result = TRUE)
  model <- brm(
    formula = bf(category ~ mean_Dd_100k + mean_Np_100k +  mean_Oa_100k + mean_Ocp_100k+  
                   mean_Ship_100k +  mean_Sst_100k +    
                   Turbidity.Mean_100k + Turbidity.Coefficient.of.Variance_100k +
                   (1 | bioregion)) +
      lf(disc ~ 0+lh, cmc = FALSE), 
    prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
              prior(normal(0, 2), class = b),          # weakly informative
              prior(exponential(1), class = sd)),         # weakly informative
    iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
    data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

  model_list[[i]] <- model  
  
}

# ----------------------------------------------------------------------------------------------------
# FIGURE S5
# ----------------------------------------------------------------------------------------------------
plotlist <- vector(mode = "list", length = length(categories))

for(i in 1:5){
  
  probs <- p_direction(
    model_list[[i]],
    effects = c("all"),
    component = c("all"),
    parameters = NULL,
    method = "direct"
  )
  
  # Get probs to add to parameter plot 
  addtoplot <- data.frame(probs$pd[5:15])
  
  par_plot <-  mcmc_plot(model_list[[i]])
  plotme <- par_plot$data[5:15,] %>% droplevels()
  plotme <- cbind(plotme, addtoplot)
  colnames(plotme)[10] <- c("Probability")
  plotme$Probability <- sprintf("%.2f", plotme$Probability)
  levels(plotme$parameter) <- c("Destructive demersal fishing",
                                "Nutrient pollution",          
                                "Ocean acidification",         
                                "Organic chemical pollution",  
                                "Shipping",                    
                                "Extreme sea surface temperature events",
                                "Turbidity (mean)", 
                                "Turbidity (CV)",
                                "Opportunistic", 
                                "Mixed",
                                "Persistent")
  g1 <- plotme[1:8,] %>%
    arrange(desc(m)) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
    mutate(parameter=factor(parameter, levels=parameter)) %>%   # This trick update the factor levels
    ggplot(aes(x=m, y=parameter)) +
    geom_pointrange(aes(xmin = ll, xmax = hh), color='grey27', shape=21, fatten = 2, size = 0.5) +
    geom_pointrange(aes(xmin = l, xmax = h), fill='dodgerblue1', color='grey52', shape=21, fatten = 2, size = 2) +
    geom_point(size=3, color="dodgerblue1") +
    theme_bw() +
    xlab("")+
    ylab("")+
    ggtitle(label = "Pressures \n \n Trajectory category breaks =", subtitle = paste0(categories[[i]], collapse = " "))+
    xlim(-5,5)+
    geom_text(aes(x = 4, label = Probability))+
    geom_vline(aes(xintercept=0), linetype = 'dotted')
  
  
  g2 <- plotme[9:11,] %>%
    arrange(desc(m)) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
    mutate(parameter=factor(parameter, levels=parameter)) %>%   # This trick update the factor levels
    ggplot(aes(x=m, y=parameter)) +
    geom_pointrange(aes(xmin = ll, xmax = hh), color='grey27', shape=21, fatten = 2, size = 0.5) +
    geom_pointrange(aes(xmin = l, xmax = h), fill='blue', color='grey52', shape=21, fatten = 2, size = 2) +
    geom_point(size=3, color="dodgerblue1") +
    theme_bw() +
    xlab("Parameter estimate")+
    ylab("")+
    ggtitle("Life history")+
    xlim(-5,5)+
    geom_text(aes(x = 4, label = Probability))+
    geom_vline(aes(xintercept=0), linetype = 'dotted')+
    coord_fixed(ratio = 1.3)
  plot <- g1/g2
  
plotlist[[i]] <- plot
}

for(i in 1:5){
  ggsave(plot = plotlist[[i]], file = paste("Figures/", "model",paste0(categories[[i]], collapse = " "),".png"))
}

