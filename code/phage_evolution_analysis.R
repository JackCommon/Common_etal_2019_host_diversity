##### phage_evolution_analysis.R by Jack Common
#### This performs the analysis of the binary phage evolution data. Bayesian GLMMs are used
#### to generate analyses and figures used in the manuscript, as conversion problems (due to all zeros
#### in the 24-clone treatments) in lme4 models mean they are pretty useless. I show these models here
#### so that hopefully you can see they don't work usefully!

rm(list=ls())
options(show.signif.stars = F)

#### ---- Dependencies ---- ####
library(tidyverse)
library(magrittr)
library(scales)
library(ggstatsplot)
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
  "3-clone" = "3-clone",
  "6-clone" = "6-clone",
  "12-clone" = "12-clone",
  "24-clone" = "24-clone"
)

treatment_labeller <- function(variable, value) {
  return(treatment_names[value])
}

treatment_names <- c("3-clone", "6-clone", "12-clone", "24-clone")

timepoint_names <- list(
  "1" = "1 dpi",
  "2" = "2 dpi",
  "3" = "3 dpi"
)

timepoint_labeller <- function(variable, value) {
  return(timepoint_names[value])
}

#### ---- Original data ---- ####
data <- read.csv("./original_data/phage_evolution.csv", header=T)
data$Timepoint %<>% as.factor
data$Replicate %<>% as.factor
data$Phage.isolate %<>% as.factor

data$Treatment %<>% relevel(ref="24-clone")
data$Treatment %<>% relevel(ref="12-clone")
data$Treatment %<>% relevel(ref="6-clone")
data$Treatment %<>% relevel(ref="3-clone")

data <- data %>% 
  select(-Other.ID.1, -Other.ID.2)

data$expanded <- rep(0,length(data$Treatment))
data$shift <- rep(0,length(data$Treatment))

data$expanded <- ifelse (data$Targetted == 1 & data$Other==1, 1, 0)
data$shift <- ifelse (data$Targetted == 0 & data$Other==1, 1, 0)

oth <- data %>% 
  #select(-Other.ID.1, -Other.ID.2) %>% 
  na.exclude

summarised <- data %>% 
  group_by(Treatment, Timepoint, Replicate) %>% 
  summarise(targ= mean(Targetted), non.targ = mean(Other), exp=mean(expanded), shifted=mean(shift))
summarised$Timepoint %<>% as.factor


#### ---- Models ---- ####

m.null <- glmer(expanded~1+(1|Replicate), data=data, family=binomial)
m1 <- glmer(expanded~Timepoint+(1|Replicate), data=data, family=binomial)
m2 <- glmer(expanded~Treatment+(1|Replicate), data=data, family=binomial)
m3 <- glmer(expanded~Timepoint+Treatment+(1|Replicate), data=filter(data, Treatment!="24-clone"),
            family=binomial)
m4 <- glm(expanded~Timepoint*Treatment, data=data, family=binomial)

anova(m.null, m1, m2, m3, m4, test="LRT")
AIC(m.null, m1, m2, m3, m4) %>% compare_AICs()

MuMIn::r.squaredGLMM(m.null)
MuMIn::r.squaredGLMM(m1)
MuMIn::r.squaredGLMM(m2)
MuMIn::r.squaredGLMM(m3)
MuMIn::r.squaredGLMM(m4)

anova(m3, test="LRT")
summary(m3)
summary(m4)

ggcoefstats(m3, exclude.intercept = F) 

# Lots of conversion problems here, so let's use a Bayesian approach

#### Bayesian models ####
library(MCMCglmm)
library(VCVglmm)
library(data.table)
library(plyr)
library(ape)
library(coda)
library(lattice)
library(aod)

d <- data %>% 
  na.exclude

d$Treatment %<>% relevel(ref="24-clone")
d$Treatment %<>% relevel(ref="12-clone")
d$Treatment %<>% relevel(ref="12-clone")
d$Treatment %<>% relevel(ref="6-clone")
d$Treatment %<>% relevel(ref="3-clone")

prior1 <- list(R = list(V = 1, fix = 1), G = list(G1 = list(V = 1, nu = 0.002)))

m1 <- MCMCglmm(expanded~Timepoint, random=~Replicate,
               family="threshold", prior=prior1,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE, data=d)

m2 <- MCMCglmm(expanded~Treatment, random=~Replicate,
               family="threshold", prior=prior1,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE, data=d)

m3 <- MCMCglmm(expanded~Timepoint+Treatment, random=~Replicate, 
               family="threshold", prior=prior1,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE, 
               data=d)

m4 <- MCMCglmm(expanded~Timepoint*Treatment, random=~Replicate, 
               family="threshold", prior=prior1,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE, 
               data=d)

Wald.test.auto(m1)
Wald.test.auto(m2)
Wald.test.auto(m3)

### Plotting predictions ####
d$m4.means <- predict(m4, d, type="response")
qplot(y=m4.means, x=Timepoint, data=d)+
  geom_point(aes(colour=Treatment))

d$m4.95.lower <- predict(m4, d, type="response",
                                interval="confidence", level=0.95)[,"lwr"]
d$m4.95.upper <- predict(m4, d, type="response",
                                interval="confidence", level=0.95)[,"upr"]
d$m4.89.lower <- predict(m4, d, type="response",
                                interval="confidence", level=0.89)[,"lwr"]
d$m4.89.upper <- predict(m4, d, type="response",
                                interval="confidence", level=0.89)[,"upr"]
d$m4.67.lower <- predict(m4, d, type="response",
                                interval="confidence", level=0.67)[,"lwr"]
d$m4.67.upper <- predict(m4, d, type="response",
                                interval="confidence", level=0.67)[,"upr"]

preds <- d %>% 
  select(Treatment, Timepoint, m4.means, 
         m4.95.lower, m4.95.upper,
         m4.89.lower, m4.89.upper,
         m4.67.lower, m4.67.upper) %>% 
  unique

pd <- position_dodge(0.5)

write.csv(x = preds, file = "./summary_data/phage_evo_preds.csv", row.names = F)

p <- ggplot(aes(x=Treatment, y=m4.means), data=preds)+
  labs(x="CRISPR allele diversity", y="Proportion of phage with\n expanded host range")+
  geom_errorbar(aes(ymin=m4.67.lower, ymax=m4.67.upper), 
                width=0, size=4, alpha=0.5, position=pd)+
  geom_errorbar(aes(ymin=m4.89.lower, ymax=m4.89.upper), 
                width=0, size=2, alpha=0.5, position=pd)+
  geom_errorbar(aes(ymin=m4.95.lower, ymax=m4.95.upper), 
                width=0, position=pd)+
  geom_point(aes(position=Treatment), pch=21, fill="white", 
             colour="black", position=pd, size=2)+
  facet_wrap(~Timepoint, labeller = timepoint_labeller)+
  cowplot::theme_cowplot()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=16),
        panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"))+
  NULL

p

ggsave("Figure_2.png", p, path="./figs/", 
       device='png', dpi=600, width=25, height=10, units=c("cm"))
