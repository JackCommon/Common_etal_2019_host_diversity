#### host_density_figs.R by Jack Common 
#### Code to generate figures that display the raw host density data

rm(list=ls())
options(show.signif.stars = F)

#### ---- Dependencies ---- ####
library(tidyverse)
library(magrittr)
library(scales)
library(lme4)
library(cowplot)
library(lmerTest)
library(ggdark)

#### ---- Functions ---- ####
## Compare AIC values for model fit
# This function extracts the AIC for each GLM, and then compares the absolute relative 
# differences for each AIC to the model with the lowest AIC. This acts as a measure of 
# model fit. More can be found at: http://faculty.washington.edu/skalski/classes/QERM597/papers_xtra/Burnham%20and%20Anderson.pdf

compare_AICs = function(df){          # df is a dataframe of AIC values 
  print(df)                           # prints the origina AIC values 
  col_len = length(df[,2])            # extracts the number of number of models
  AIC_min = abs(min(df[,2]))          # finds the minimum AIC value
  for (i in seq(1, col_len, 1)){      # loop through the AIC values and prints the absolute differences from AIC_min
    print( (abs(df[i,2])) - AIC_min)
  }
}

## Convert log-odds to probabilities from binomial GLM output
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
treatment_names <- list(
  "1-clone" = "1-clone",
  "3-clone" = "3-clone",
  "6-clone" = "6-clone",
  "12-clone" = "12-clone",
  "24-clone" = "24-clone",
  "1-clone_control" = "1-clone (ancestral phage)",
  "24-clone_control" = "24-clone (ancestral phage)"
)

treatment_labeller <- function(variable, value) {
  return(treatment_names[value])
}

treatment_names <- c("1-clone", "3-clone", "6-clone", "12-clone", "24-clone",
                     "1-clone (ancestral phage)", "24-clone (ancestral phage)")

pfu_labels <- c(
  expression("0"), 
  expression('10'^2*''), 
  expression('10'^4*''),
  expression('10'^6*''),
  expression('10'^8*''),
  expression('10'^10*''),
  expression("10"^12*"")
)



#### ---- Original data ---- ####
data_original <- read.csv("./original_data/dynamics_master.csv", header=T, skip=1)  %>% 
  select(Treatment, Timepoint, Replicate, BIM_mix, Escape_phage, Tracked_BIM,
         PFU, SM_CFU, CRISPR_CFU, BIM_CFU, white_CFU, total_CFU,
         w_SM, w_CRISPR, w_BIM, w_white)

data_original$Timepoint %<>% as.factor()
data_original$Replicate %<>% as.factor
data_original$Tracked_BIM %<>% as.factor()

# Host titre dataset
host <- data_original %>% 
  select(-PFU, -w_SM, -w_CRISPR, -w_BIM, -w_white, -white_CFU, -total_CFU, 
         SM=SM_CFU, CRISPR=CRISPR_CFU, BIM=BIM_CFU) %>% 
   gather("SM", "CRISPR", "BIM", key="Strain", value="CFU", factor_key = T)

host$Treatment %<>% relevel(ref="24-clone_control")
host$Treatment %<>% relevel(ref="1-clone_control")
host$Treatment %<>% relevel(ref="24-clone")
host$Treatment %<>% relevel(ref="12-clone")
host$Treatment %<>% relevel(ref="6-clone")
host$Treatment %<>% relevel(ref="3-clone")
host$Treatment %<>% relevel(ref="1-clone")

# make a new replicate column for graphic purposes
host$Replicate_2 <- rep("0", nrow(host))
for (i in 1:nrow(host)){
  for(n in 1:24){
    if(host$Strain[i]=="SM" & host$Replicate[i]==n)
      host$Replicate_2[i] <- paste(n, "SM", sep=".")
    if(host$Strain[i]=="CRISPR" & host$Replicate[i]==n)
      host$Replicate_2[i] <- paste(n, "CRISPR", sep=".")
    if(host$Strain[i]=="BIM" & host$Replicate[i]==n)
      host$Replicate_2[i] <- paste(n, "BIM", sep=".")
  }
}
host$Replicate_2 %<>% as.factor()


#### --- Figures ---- ####
library(wesanderson)
pal <- wes_palette("IsleofDogs1", 3)

host_density_fig <- ggplot(aes(y=CFU+1, x=Timepoint, group=Replicate_2), 
                         data=host)+
  
  #geom_point(stat='identity')+
  geom_line(stat='identity', aes(colour=Strain))+
  facet_wrap(~Treatment, labeller = treatment_labeller)+
  
  labs(x="Days post-infection", y=expression(bold("Host density (Cfu ml"*{}^{-1}*")")))+
  #ggtitle("Density of CRISPR and SM clones ")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        strip.text = element_text(face='bold', size=14),
        axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(face="bold", size=16),
        legend.title.align = 0.5,
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key = element_rect(fill="grey95"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  #coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  scale_colour_manual(values=pal,
                      labels=c("Surface mutant\n(PA14 ∆pilA)", "Resistant\nCRISPR clones",
                               "Susceptible\nCRISPR clones"),
                      name=c("Host immunity"))+
   NULL
# quartz()
# last_plot()

ggsave("Figure_S4.png", host_density_fig, path="./figs/", 
              device="png", dpi=600, width=28, height = 20, units = c("cm"))
