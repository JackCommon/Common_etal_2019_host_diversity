#### phage_analysis.R by Jack Common 
#### Code to generate figures that display the raw phage data, and also conduct
#### statistical analysis (using GLMMs) of how titre changes over time and between
#### CRISPR diversity treatments

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
    "1-clone_control" = "1-clone\n(ancestral phage)",
    "24-clone_control" = "24-clone\n(ancestral phage)"
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

# Phage titre dataset
phage <- data_original %>% 
  select(-SM_CFU, -CRISPR_CFU, -BIM_CFU, -white_CFU, -w_SM, -w_CRISPR, -w_BIM, -w_white) %>% 
  gather("PFU", key="Strain", value="Titre", factor_key = T)


# Relevel the datasets so things make sense
data_original$Treatment %<>% relevel(ref="24-clone_control")
data_original$Treatment %<>% relevel(ref="1-clone_control")
data_original$Treatment %<>% relevel(ref="24-clone")
data_original$Treatment %<>% relevel(ref="12-clone")
data_original$Treatment %<>% relevel(ref="6-clone")
data_original$Treatment %<>% relevel(ref="3-clone")
data_original$Treatment %<>% relevel(ref="1-clone")

phage$Treatment %<>% relevel(ref="24-clone_control")
phage$Treatment %<>% relevel(ref="1-clone_control")
phage$Treatment %<>% relevel(ref="24-clone")
phage$Treatment %<>% relevel(ref="12-clone")
phage$Treatment %<>% relevel(ref="6-clone")
phage$Treatment %<>% relevel(ref="3-clone")
phage$Treatment %<>% relevel(ref="1-clone")

write.csv(phage, "./original_data/phage_data.csv", row.names = F)

for(i in 1:length(data_original$w_BIM)){
  if(data_original$w_BIM[i]=="#DIV/0!"){
    data_original$w_BIM[i] <- NA
  } 
}

#### ---- Visualise phage titre data ---- ####
phage_plot <- ggplot(aes(y=log10(Titre+1), x=Timepoint, group=Replicate), data=phage)+
  
  #geom_point(stat='identity')+
  geom_path(stat='identity', size=.3)+
  geom_hline(yintercept = 100, linetype=2, size=.3)+
  facet_wrap(~Treatment, labeller = treatment_labeller, scales = "free_x")+
  
  labs(x="Days post-infection", y=expression(bold("Phage density (Pfu ml"*{}^{-1}*")")))+
  #ggtitle("Density of phage")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 10),
        axis.title = element_text(face="bold", size=10),
        axis.text = element_text(size=8),
        axis.text.x = element_text(margin=margin(0,0,0,0,"pt")),
        legend.text = element_text(size=8),
        strip.text.x = element_text(face='bold', size=8, margin=margin(1.6,0,1.6,0, "pt")),
        strip.background = element_rect(fill="transparent", colour=NA))+
  
  coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(breaks=c(seq(0,12,2)),
                     labels=pfu_labels)+
  NULL

#

last_plot()

#ggsave("Figure_1.png", phage_plot, path="./figs/", 
#       device="png", dpi=600,width=28, height = 20, units = c("cm"))

#### ---- Original data: continuous time---- ####
data_original <- read.csv("./original_data/dynamics_master.csv", header=T, skip=1)  %>% 
  select(Treatment, Timepoint, Replicate, BIM_mix, Escape_phage, Tracked_BIM,
         PFU, SM_CFU, CRISPR_CFU, BIM_CFU, white_CFU, total_CFU,
         w_SM, w_CRISPR, w_BIM, w_white)

data_original$Replicate %<>% as.factor
data_original$Tracked_BIM %<>% as.factor()

# Phage titre dataset
phage <- data_original %>% 
  select(-SM_CFU, -CRISPR_CFU, -BIM_CFU, -white_CFU, -w_SM, -w_CRISPR, -w_BIM, -w_white) %>% 
  gather("PFU", key="Strain", value="Titre", factor_key = T)


# Relevel the datasets so things make sense
phage$Treatment %<>% relevel(ref="24-clone_control")
phage$Treatment %<>% relevel(ref="1-clone_control")
phage$Treatment %<>% relevel(ref="24-clone")
phage$Treatment %<>% relevel(ref="12-clone")
phage$Treatment %<>% relevel(ref="6-clone")
phage$Treatment %<>% relevel(ref="3-clone")
phage$Treatment %<>% relevel(ref="1-clone")


#### --- Continuous time models ---- #####
# Use this dataframe to analyse just the polyclonal treatments
# d.cont.phage <- phage %>%
#   na.exclude %>%
#   filter(Timepoint!="0",
#          Treatment %in% c("3-clone", "6-clone", "12-clone", "24-clone", "24-clone_control"))

m1 <- lmer(log(Titre+1)~Treatment+(1|Replicate), data=phage)

m2 <- lmer(log(Titre+1)~Timepoint+(1|Replicate), data=phage)

m3 <- lmer(log(Titre+1)~Treatment + Timepoint+(1|Replicate), data=phage)

m4 <- lmer(log(Titre+1)~Treatment * Timepoint+(1|Replicate), data=phage)

AIC(m1, m2, m3, m4) %>% compare_AICs()

plot(m1)
plot(m2)
plot(m3)
plot(m4)

summary(m3)
anova(m4, type="marginal")
confint(m3, parm="beta_", method="boot")

#### ---- Figures ----####

# Set up target dataframe to store means and CIs
m3 <- lmer(log(Titre+1)~Treatment + Timepoint+(1|Replicate), data=phage)

coefs <- data.frame(term =  factor(8), 
                    beta = numeric(8), 
                    l.95 = numeric(8), h.95 = numeric(8), 
                    l.89 = numeric(8), h.89 = numeric(8), 
                    l.67 = numeric(8), h.67 = numeric(8))

coefs$term <- c("Intercept", "3-clone", "6-clone", "12-clone", "24-clone", "1-clone (ancestral phage)",
                "24-clone (ancestral phage)","Days post-infection") %>% 
  as.factor

# Get coefficients and confidence intervals
coefs$beta <- fixef(m3)
coefs$l.95 <- confint(m3, parm="beta_", level=0.95) %>% 
  as.data.frame() %>%    
  slice(1:8) %>% 
  select("2.5 %") %>% 
  .[1:8,1]
coefs$h.95 <- confint(m3, parm="beta_", level=0.95) %>% 
  as.data.frame() %>%    
  slice(1:8) %>% 
  select("97.5 %") %>% 
  .[1:8,1]

coefs$l.89 <- confint.merMod(m3, parm="beta_", level=0.89) %>% 
  as.data.frame() %>%    
  slice(1:8) %>% 
  select("5.5 %") %>% 
  .[1:8,1]
coefs$h.89 <- confint(m3, parm="beta_", level=0.89) %>% 
  as.data.frame() %>%    
  slice(1:8) %>% 
  select("94.5 %") %>% 
  .[1:8,1]

coefs$l.67 <- confint.merMod(m3, parm="beta_", level=0.67) %>% 
  as.data.frame() %>%    
  slice(1:8) %>% 
  select("16.5 %") %>% 
  .[1:8,1]
coefs$h.67 <- confint(m3, parm="beta_", level=0.67) %>% 
  as.data.frame() %>%    
  slice(1:8) %>% 
  select("83.5 %") %>% 
  .[1:8,1]

# write.csv(coefs, "./summary_data/phage_model_coefs.csv", row.names = F)
coefs <- read.csv("./summary_data/phage_model_coefs.csv")

coefs$term %<>% relevel(ref="Days post-infection")
coefs$term %<>% relevel(ref="24-clone (ancestral phage)")
coefs$term %<>% relevel(ref="24-clone")
coefs$term %<>% relevel(ref="12-clone")
coefs$term %<>% relevel(ref="6-clone")
coefs$term %<>% relevel(ref="3-clone")
coefs$term %<>% relevel(ref="1-clone (ancestral phage)")
coefs$term %<>% relevel(ref="Intercept")

p1 <- ggplot(aes(y=beta, x=term), data=coefs)+
  geom_errorbar(aes(ymin=l.67, ymax=h.67), width=0, size=3, alpha=0.5)+
  geom_errorbar(aes(ymin=l.89, ymax=h.89), width=0, size=1.5, alpha=0.5)+
  geom_errorbar(aes(ymin=l.95, ymax=h.95), width=0)+
  geom_point(fill="white", pch=21, colour="black", size=1)+
  geom_hline(yintercept=0, linetype=2, size=.4)+
  coord_flip()+
  cowplot::theme_cowplot()+
  labs(y=expression(bold("Phage density (ln[pfu ml"*{}^{-1}*"])")), x="Fixed effect level")+ 
  scale_y_continuous(breaks=seq(-14, 20, 2))+
  theme(#panel.grid.major = element_line(colour="lightgrey"),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=10))+
  NULL
p1

#ggsave("Figure_S2.png", p1, path="./figs/", 
#       device="png", dpi=600,width=20, height=12, units=c("cm"))

#### --- Save fig ---- ####
Figure1 <- plot_grid(phage_plot, p1, labels = c("A", "B"), label_size = 10,
                     ncol=1, nrow=2, rel_heights = c(1.7, 1))
# last_plot()

ggsave("Figure_1.tiff", Figure1, path="./figs/", compression="lzw",
       device="tiff", dpi=300, width=11, height=16, units=c("cm"))
