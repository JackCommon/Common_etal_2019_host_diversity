##### BIM dynamics experiment - escape phage analysis
#### Created: 11/6/19

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
data <- read.csv("./dynamics_exp/original_data/escape_phage_assay_moretimepoints.csv", header=T)
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

#### ---- Graphs ---- ####
p1 <- ggplot(aes(y=exp, x=Timepoint), data=summarised)+
  #geom_violin(aes(colour=Timepoint))+
  geom_boxplot()+
  geom_point()+
  #stat_summary(fun.y=mean, geom="point", size=3, pch=21,
  #             aes(fill=Timepoint), position=position_dodge(.9))+
  labs(y="Proportion", x="Diversity")+
  cowplot::theme_cowplot()+
  #scale_colour_discrete(name="Days\npost-infection")+
  facet_wrap(~Treatment, labeller = treatment_labeller)+
  cowplot::theme_cowplot()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=16),
        panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"))+
  NULL
p1

ggsave("Fig. 2.png", p1, device="png", path="./docs/",
       width=18, height=12, units=c("cm"), dpi=300)

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

m5 <- glm(expanded~Targetted, data=oth, family=binomial)
summary(m5)

ggcoefstats(m5)

exp.only <- filter(oth, expanded==1)
nrow(exp.only)/nrow(oth)*100

shift.only <- filter(oth, shift==1)
nrow(shift.only)/nrow(oth)*100

#### Bayesian attempts ####
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

prior1 <- list(G = list(G1 = list(V = diag(1), nu = 0.002)),
               R = list(V = 0.5, nu = 0.002, fix = TRUE))

prior2 <- list(R = list(V = 1, fix = 1), G = list(G1 = list(V = 1, nu = 0.002)))

m1 <- MCMCglmm(expanded~Timepoint, random=~Replicate,
               family="threshold", prior=prior2,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE, data=d)

m2 <- MCMCglmm(expanded~Treatment, random=~Replicate,
               family="threshold", prior=prior2,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE, data=d)

m3 <- MCMCglmm(expanded~Timepoint+Treatment, random=~Replicate, 
               family="threshold", prior=prior2,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE, 
               data=d)

m4 <- MCMCglmm(expanded~Timepoint*Treatment, random=~Replicate, 
               family="threshold", prior=prior2,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE, 
               data=d)

Wald.test.auto(m1); Wald.test.auto(m2); Wald.test.auto(m3)
m6$DIC; m7$DIC; m8$DIC

Solapply(m8)
MCMCRepnorm2(m8)
plot(m8)
MCMCfixplot(m8)

(sum.m6 <- summary(m6))
plot(m6)

MCMCfixplot(m6)
MCMCranef(m6)

means <- apply(m6$Sol, 2, mean)
posterior.mode(m6$Sol)
HPDinterval(m6$Sol, prob=0.95)
posterior.mode(m6$VCV)
HPDinterval(m6$Sol, prob=.89)

posterior.mode(m8$VCV[, "Replicate"]/(m8$VCV[, "Replicate"] + m8$VCV[, "units"]))

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

p2 <- ggplot(aes(x=Treatment, y=m4.means), data=preds)+
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

p2

ggsave("Figure_2.png", p2, path="~/Documents/OneDrive - University of Exeter/Papers/Common, Walker-Sunderhauf and Westra 2019/", 
       device='png', dpi=600, width=25, height=10, units=c("cm"))


### ---- BLOCK OUT ---- ####
# pred.m8.95 <- predict(m8, d, type="response", 
#                    interval="confidence", level=0.95) %>% 
#   unique %>% 
#   as.data.frame() %>% 
#   rename(lwr.95=lwr, upr.95=upr)
# 
# pred.m8.89 <- predict(m8, d, type="response", 
#                       interval="confidence", level=0.89) %>% 
#   unique %>% 
#   as.data.frame() %>% 
#   select(lwr, upr) %>% 
#   rename(lwr.89=lwr, upr.89=upr)
# 
# pred.m8.67 <- predict(m8, d, type="response", 
#                       interval="confidence", level=0.67) %>% 
#   unique %>% 
#   as.data.frame() %>% 
#   select(lwr, upr) %>% 
#   rename(lwr.67=lwr, upr.67=upr)
# 
# pred.m8 <- bind_cols(pred.m8.95, pred.m8.89, pred.m8.67)
# 
# pred.treats <- c(rep(c("3-clone", "6-clone", "12-clone", "24-clone"),3))
# pred.times <- c(rep(1, 4), rep(2, 4), rep(3,4))
# pred.vars <- data.frame(pred.treats, pred.times)
# 
# pred.m8.sum <- bind_cols(pred.vars, pred.m8)
# pred.m8.sum$pred.times %<>% as.factor
# pred.m8.sum$pred.treats %<>% relevel(ref="24-clone")
# pred.m8.sum$pred.treats %<>% relevel(ref="12-clone")
# pred.m8.sum$pred.treats %<>% relevel(ref="6-clone")
# pred.m8.sum$pred.treats %<>% relevel(ref="3-clone")
# 
# pd <- position_dodge(0.5)
# 
# p1 <- ggplot(aes(x=pred.treats, y=fit), data=pred.m8.sum)+
#   labs(x="Diversity", y="Proportion")+
#   geom_errorbar(aes(ymin=lwr.67, ymax=upr.67, colour=pred.times), 
#                 width=0, size=3, alpha=0.4, position=pd)+
#   geom_errorbar(aes(ymin=lwr.89, ymax=upr.89, colour=pred.times), 
#                 width=0, size=3, alpha=0.4, position=pd)+
#   geom_errorbar(aes(ymin=lwr.95, ymax=upr.95, colour=pred.times), 
#                 width=0, size=3, alpha=0.4, position=pd)+
#   geom_point(aes(position=pred.times), pch=21, colour="white", position=pd, size=2)+
#   coord_cartesian(ylim=c(0,1))+
#   scale_colour_discrete(name="Days\npost-infection")+
#   dark_theme_bw()+
#   NULL
# 
# p1
# 
# p1 <- ggplot(aes(x=Treatment, y=Posterior.Mean), data=est.data)+
#   #geom_hline(yintercept=0, linetype=2)+
#   labs(x="Diversity", y="Proportion")+
#   geom_errorbar(aes(ymin=X95.l, ymax=X95.h, colour=Timepoint), 
#                 width=0, size=3, alpha=0.4, position=pd)+
#   geom_errorbar(aes(ymin=X89.l, ymax=X89.h, colour=Timepoint), 
#                 width=0, size=3, alpha=0.4, position=pd)+
#   geom_errorbar(aes(ymin=X67.l, ymax=X67.h, colour=Timepoint), 
#                 width=0, size=3, alpha=0.4, position=pd)+
#   geom_point(aes(position=Timepoint), pch=21, colour="white", position=pd, size=2)+
#   scale_colour_discrete(name="Days\npost-infection")+
#   dark_theme_bw()+
#   theme(
#     legend.title.align = 0.5)+
#   NULL
# p1
# 
# last_plot()
# 
# means <- apply(m6$Sol, 2, mean)
# HPD.95 <- HPDinterval(m6$Sol, prob=0.95)
# HPD.89 <- HPDinterval(m6$Sol, prob=0.89)
# HPD.67 <- HPDinterval(m6$Sol, prob=0.67)
# 
# mus<- c()
# HPD.95.lower <- c()
# HPD.89.lower <- c()
# HPD.67.lower <- c()
# HPD.95.upper <- c()
# HPD.89.upper <- c()
# HPD.67.upper <- c()
# 
# for (i in 1:6){
#   if(i==1) { 
#     mus[i]      <- means[i]
#     HPD.95.lower[i] <- HPD.95[i,"lower"]
#     HPD.95.upper[i] <- HPD.95[i,"upper"]
#     HPD.89.lower[i] <- HPD.89[i,"lower"]
#     HPD.89.upper[i] <- HPD.89[i,"upper"]
#     HPD.67.lower[i] <- HPD.67[i,"lower"]
#     HPD.67.upper[i] <- HPD.67[i,"upper"]
#   } else{
#   mus[i] <- means[1]+means[i]
#   HPD.95.lower[i] <- HPD.95[1, "lower"]+HPD.95[i,1]
#   HPD.95.upper[i] <- HPD.95[1, "upper"]+HPD.95[i,2]
#   HPD.89.lower[i] <- HPD.89.lower[1]+HPD.89[i,1]
#   HPD.89.upper[i] <- HPD.89.upper[1]+HPD.89[i,2]
#   HPD.67.lower[i] <- HPD.67.lower[1]+HPD.67[i,1]
#   HPD.67.upper[i] <- HPD.67.upper[1]+HPD.67[i,2]
#   }
# }
# 
# ests <- data.frame(mus, HPD.95.lower, HPD.95.upper,
#                         HPD.89.lower, HPD.89.upper,
#                         HPD.67.lower, HPD.67.upper
# )
# 
# ests
# clip <- pipe("pbcopy", "w")                       
# write.table(ests, sep = "\t", file=clip, row.names = F)                               
# close(clip)
# 
# 
# est.data <- read.csv("./Bayesian stuff/escape phage MCMCglmm estimates 2.csv")
# est.data$Timepoint %<>% as.factor
# est.data$Treatment %<>% relevel(ref="24-clone")
# est.data$Treatment %<>% relevel(ref="12-clone")
# est.data$Treatment %<>% relevel(ref="6-clone")
# est.data$Treatment %<>% relevel(ref="3-clone")
# 
# pd <- position_dodge(.5)
# 
# p1 <- ggplot(aes(x=Treatment, y=Posterior.Mean), data=est.data)+
#   #geom_hline(yintercept=0, linetype=2)+
#   labs(x="Diversity", y="Proportion")+
#   geom_errorbar(aes(ymin=X95.l, ymax=X95.h, colour=Timepoint), 
#                 width=0, size=3, alpha=0.4, position=pd)+
#   geom_errorbar(aes(ymin=X89.l, ymax=X89.h, colour=Timepoint), 
#                 width=0, size=3, alpha=0.4, position=pd)+
#   geom_errorbar(aes(ymin=X67.l, ymax=X67.h, colour=Timepoint), 
#                 width=0, size=3, alpha=0.4, position=pd)+
#   geom_point(aes(position=Timepoint), pch=21, colour="white", position=pd, size=2)+
#   scale_colour_discrete(name="Days\npost-infection")+
#   dark_theme_bw()+
#   theme(
#     legend.title.align = 0.5)+
#   NULL
# p1
# 
# cowplot::ggsave("phage expanded HPD.png", p1, device="png", path="./docs/",
#                 dpi=300, width=15, height=8, units=c("cm"))
