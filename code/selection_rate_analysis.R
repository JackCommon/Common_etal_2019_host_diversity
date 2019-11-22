#### selection_rate_analysis.R by Jack Common
#### Code to generate figures that display the raw selection rate data, conduct
#### statistical analysis (using GLMMs) of how selection rate changes over time and between
#### CRISPR diversity treatments, and generate figures that display the means predicted from
#### those models

rm(list=ls())
options(show.signif.stars = F)

#### ---- Dependencies ---- ####
library(scales)
library(cowplot)
library(ggplot2)
library(plyr)
library(lme4)
library(magrittr)
library(tidyverse)
library(lmerTest)

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

timepoint_names <- list(
  "0" = "0 dpi",
  "1" = "1 dpi",
  "2" = "2 dpi",
  "3" = "3 dpi"
)

timepoint_labeller <- function(variable, value) {
  return(timepoint_names[value])
}

treatment_list <- c("1-clone", "1-clone\n(ancestral phage)", "3-clone", "6-clone", "12-clone",
                    "24-clone", "24-clone\n(ancestral phage)")

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
d <- read.csv("./original_data/selection rates.csv", header=T)  %>% 
  select(Treatment, Timepoint, Replicate, BIM_mix, Escape_phage, Tracked_BIM,
         rCRISPR, rBIM) %>% 
  gather("rCRISPR", "rBIM", key="Strain", value="r", factor_key = T)

d$Timepoint %<>% as.factor()
d$Replicate %<>% as.factor
d$Tracked_BIM %<>% as.factor()
d$Escape_phage %<>% as.factor

# Relevel the datasets so things make sense
d$Treatment %<>% relevel(ref="24-clone_control")
d$Treatment %<>% relevel(ref="24-clone")
d$Treatment %<>% relevel(ref="12-clone")
d$Treatment %<>% relevel(ref="6-clone")
d$Treatment %<>% relevel(ref="3-clone")
d$Treatment %<>% relevel(ref="1-clone_control")
d$Treatment %<>% relevel(ref="1-clone")

#### ---- CRISPR boxplot ---- ####
pd <- position_dodge(0.5)
library(wesanderson)
pal <- wes_palette("IsleofDogs1", 3)

CRISPR_fitness_plots <- ggplot(aes(x=Treatment, y=r), 
                               data=filter(d, Strain=="rCRISPR", Timepoint!="0"))+
  geom_boxplot(width=0.5, outlier.shape = NA, aes(colour=Timepoint), position = pd)+
  geom_point(aes(colour=Timepoint), alpha=0.5, pch=21, position = pd)+
  geom_hline(yintercept=0, linetype=2, size=1)+
  #stat_summary(fun.y="mean", geom="point", colour="red")+
  #facet_wrap(~Timepoint, scales="free", labeller = timepoint_labeller, ncol=2)+
  labs(y="Selection rate of\nall CRISPR clones", x="CRISPR diversity")+
  scale_x_discrete(labels = treatment_list)+
  scale_colour_manual(values=pal,
                      name=c("Days\npost-infection"))+
  theme_cowplot()+
  theme(axis.title = element_text(face="bold", size=10),
        #panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"),
        axis.text = element_text(size=8),
        legend.text = element_text(size=10),
        legend.title = element_text(face="bold", size=10),
        legend.title.align = 0.5,
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key = element_rect(fill="grey95"))+
  NULL

last_plot()

#### ---- BIM boxplot ---- ####
BIM_fitness_plots <- ggplot(aes(x=Treatment, y=r), 
                            data=filter(d, Strain=="rBIM", Timepoint!="0",
                                        Treatment%in%c("3-clone", "6-clone", "12-clone",
                                                       "24-clone", "24-clone_control")))+
  geom_boxplot(aes(colour=Timepoint), width=0.5, outlier.shape = NA, position = pd)+
  geom_point(aes(colour=Timepoint), alpha=0.5, pch=21, position=pd)+
  geom_hline(yintercept=0, linetype=2, size=1)+
  #stat_summary(fun.y="mean", geom="point", colour="red")+
  labs(y="Selection rate of\nsusceptible CRISPR clone", x="CRISPR diversity")+
  scale_x_discrete(labels= treatment_list[3:7])+
  scale_colour_manual(values=pal,
                      name=c("Days\npost-infection"))+
  theme_cowplot()+
  theme(axis.title = element_text(face="bold", size=10),
        #panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"),
        axis.text = element_text(size=8),
        legend.text = element_text(size=10),
        legend.title = element_text(face="bold", size=10),
        legend.title.align = 0.5,
        legend.key.height = unit(.4, "cm"),
        legend.key.width = unit(.4, "cm"),
        legend.key = element_rect(fill="grey95"))+
  NULL
last_plot()

FigS3 <- plot_grid(CRISPR_fitness_plots+
                     theme(legend.position = "none"), 
                   BIM_fitness_plots, ncol=1, 
                   labels = c("A", "B"), label_size = 10)


last_plot()
ggsave("Figure_S3.tif", FigS3, path="./figs/",
       device="tiff",dpi=300, width=12, height = 13, units=c("cm"), compression="lzw")


#### ---- lme4 CRISPR models ---- ####
d.CRISPR <- d %>% 
  na.exclude %>% 
  filter(Strain=="rCRISPR", Timepoint!="0")

## Use this dataframe to analyse just the polyclonal treatments
# d.CRISPR <- d %>% 
#   na.exclude %>% 
#   filter(Strain=="rCRISPR", Timepoint!="0", 
#          Treatment %in% c("3-clone", "6-clone", "12-clone", "24-clone", "24-clone_control"))

d.CRISPR$Treatment %<>% relevel(ref="24-clone_control")
d.CRISPR$Treatment %<>% relevel(ref="1-clone_control")
d.CRISPR$Treatment %<>% relevel(ref="24-clone")
d.CRISPR$Treatment %<>% relevel(ref="12-clone")
d.CRISPR$Treatment %<>% relevel(ref="6-clone")
d.CRISPR$Treatment %<>% relevel(ref="3-clone")
d.CRISPR$Treatment %<>% relevel(ref="1-clone")

d.CRISPR$Timepoint %<>% relevel(ref="3")
d.CRISPR$Timepoint %<>% relevel(ref="2")
d.CRISPR$Timepoint %<>% relevel(ref="1")
d.CRISPR$Timepoint %<>% relevel(ref="0")

m1 <- lmer(r~Treatment + (1|Replicate), data=d.CRISPR)
m2 <- lmer(r~Timepoint + (1|Replicate), data=d.CRISPR)
m3 <- lmer(r~Treatment + Timepoint +(1|Replicate), data=d.CRISPR)
m4 <- lmer(r~Treatment * Timepoint + (1|Replicate), data=d.CRISPR)

AIC(m1, m2, m3, m4) %>% 
  compare_AICs()

plot(m1)
plot(m2)
plot(m3)
plot(m4)

anova(m4, type="marginal")

summary(m1)
confint(m1, parm="beta_", level=0.95, method="boot")
ggcoefstats(m1)

#### ---- CRISPR: Extract model coefficients and CIs ---- ####
d.CRISPR$Treatment %<>% relevel(ref="24-clone_control")
d.CRISPR$Treatment %<>% relevel(ref="24-clone")
d.CRISPR$Treatment %<>% relevel(ref="12-clone")
d.CRISPR$Treatment %<>% relevel(ref="6-clone")
d.CRISPR$Treatment %<>% relevel(ref="3-clone")
d.CRISPR$Treatment %<>% relevel(ref="1-clone_control")
d.CRISPR$Treatment %<>% relevel(ref="1-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.CRISPR)

CRISPR.coefs <- data.frame(term = factor(7), 
                           mean = numeric(7), 
                           l.95 = numeric(7), h.95 = numeric(7), 
                           l.89 = numeric(7), h.89 = numeric(7), 
                           l.67 = numeric(7), h.67 = numeric(7))

CRISPR.coefs$term <- c("1-clone", '1-clone\n(ancestral phage)',
                       "3-clone", "6-clone", "12-clone",
                       "24-clone", "24-clone\n(ancestral phage)") %>% 
  as.factor

CRISPR.coefs$mean <- fixef(m1)

for(i in 2:7){
  CRISPR.coefs$mean[i] <- CRISPR.coefs$mean[1]+CRISPR.coefs$mean[i]
}

d.CRISPR$Treatment %<>% relevel(ref="1-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.CRISPR)

CRISPR.coefs$l.95[1] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.95[1] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.89[1] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.89[1] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.67[1] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.67[1] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.CRISPR$Treatment %<>% relevel(ref="1-clone_control")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.CRISPR)

CRISPR.coefs$l.95[2] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.95[2] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.89[2] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.89[2] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.67[2] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.67[2] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.CRISPR$Treatment %<>% relevel(ref="3-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.CRISPR)

CRISPR.coefs$l.95[3] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.95[3] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.89[3] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.89[3] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.67[3] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.67[3] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.CRISPR$Treatment %<>% relevel(ref="6-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.CRISPR)

CRISPR.coefs$l.95[4] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.95[4] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.89[4] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.89[4] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.67[4] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.67[4] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.CRISPR$Treatment %<>% relevel(ref="12-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.CRISPR)

CRISPR.coefs$l.95[5] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.95[5] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.89[5] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.89[5] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.67[5] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.67[5] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.CRISPR$Treatment %<>% relevel(ref="24-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.CRISPR)

CRISPR.coefs$l.95[6] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.95[6] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.89[6] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.89[6] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.67[6] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.67[6] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.CRISPR$Treatment %<>% relevel(ref="24-clone_control")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.CRISPR)

CRISPR.coefs$l.95[7] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.95[7] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.89[7] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.89[7] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
CRISPR.coefs$l.67[7] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
CRISPR.coefs$h.67[7] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
#### ---- lme4 BIM models ---- ####
d.BIM <- d %>% 
  na.exclude %>% 
  filter(Strain=="rBIM", Timepoint!="0")

## Use this dataframe to analyse just the polyclonal treatments
# d.BIM <- d %>% 
#   na.exclude %>% 
#   filter(Strain=="rBIM", Timepoint!="0", 
#          Treatment %in% c("3-clone", "6-clone", "12-clone", "24-clone", "24-clone_control"))

d.BIM$Treatment %<>% relevel(ref="24-clone_control")
d.BIM$Treatment %<>% relevel(ref="1-clone_control")
d.BIM$Treatment %<>% relevel(ref="24-clone")
d.BIM$Treatment %<>% relevel(ref="12-clone")
d.BIM$Treatment %<>% relevel(ref="6-clone")
d.BIM$Treatment %<>% relevel(ref="3-clone")
d.BIM$Treatment %<>% relevel(ref="1-clone")

d.BIM$Timepoint %<>% relevel(ref="3")
d.BIM$Timepoint %<>% relevel(ref="2")
d.BIM$Timepoint %<>% relevel(ref="1")
d.BIM$Timepoint %<>% relevel(ref="0")

m1 <- lmer(r~Treatment + (1|Replicate), data=d.BIM)
m2 <- lmer(r~Timepoint + (1|Replicate), data=d.BIM)
m3 <- lmer(r~Treatment + Timepoint +(1|Replicate), data=d.BIM)
m4 <- lmer(r~Treatment * Timepoint + (1|Replicate), data=d.BIM)

AIC(m1, m2, m3, m4) %>% 
  compare_AICs()

plot(m1)
plot(m2)
plot(m3)
plot(m4)

anova(m3, type="marginal")

summary(m1)
confint(m1, parm="beta_", method="boot")
ggcoefstats(m3)

#### ---- BIM: Extract model coefficients and CIs ---- ####
d.BIM$Treatment %<>% relevel(ref="24-clone_control")
d.BIM$Treatment %<>% relevel(ref="24-clone")
d.BIM$Treatment %<>% relevel(ref="12-clone")
d.BIM$Treatment %<>% relevel(ref="6-clone")
d.BIM$Treatment %<>% relevel(ref="3-clone")
d.BIM$Treatment %<>% relevel(ref="1-clone_control")
d.BIM$Treatment %<>% relevel(ref="1-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.BIM)

BIM.coefs <- data.frame(term = factor(7), 
                           mean = numeric(7), 
                           l.95 = numeric(7), h.95 = numeric(7), 
                           l.89 = numeric(7), h.89 = numeric(7), 
                           l.67 = numeric(7), h.67 = numeric(7))

BIM.coefs$term <- c("1-clone", '1-clone\n(ancestral phage)',
                       "3-clone", "6-clone", "12-clone",
                       "24-clone", "24-clone\n(ancestral phage)") %>% 
  as.factor

BIM.coefs$mean <- fixef(m1)

for(i in 2:7){
  BIM.coefs$mean[i] <- BIM.coefs$mean[1]+BIM.coefs$mean[i]
}

d.BIM$Treatment %<>% relevel(ref="1-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.BIM)

BIM.coefs$l.95[1] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.95[1] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.89[1] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.89[1] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.67[1] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.67[1] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.BIM$Treatment %<>% relevel(ref="1-clone_control")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.BIM)

BIM.coefs$l.95[2] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.95[2] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.89[2] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.89[2] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.67[2] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.67[2] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.BIM$Treatment %<>% relevel(ref="3-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.BIM)

BIM.coefs$l.95[3] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.95[3] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.89[3] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.89[3] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.67[3] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.67[3] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.BIM$Treatment %<>% relevel(ref="6-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.BIM)

BIM.coefs$l.95[4] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.95[4] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.89[4] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.89[4] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.67[4] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.67[4] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.BIM$Treatment %<>% relevel(ref="12-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.BIM)

BIM.coefs$l.95[5] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.95[5] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.89[5] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.89[5] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.67[5] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.67[5] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.BIM$Treatment %<>% relevel(ref="24-clone")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.BIM)

BIM.coefs$l.95[6] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.95[6] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.89[6] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.89[6] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.67[6] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.67[6] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]

d.BIM$Treatment %<>% relevel(ref="24-clone_control")
m1 <- lmer(r~Treatment + (1|Replicate), data=d.BIM)

BIM.coefs$l.95[7] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.95[7] <- confint(m1, parm="beta_", level=0.95) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.89[7] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.89[7] <- confint(m1, parm="beta_", level=0.89) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
BIM.coefs$l.67[7] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,1]
BIM.coefs$h.67[7] <- confint(m1, parm="beta_", level=0.67) %>% 
  as.data.frame() %>% 
  slice(1:7) %>% 
  .[1,2]
#### ---- Join dataframes and save ---- ####
selection_rate_summary <- bind_rows(CRISPR.coefs, BIM.coefs)
selection_rate_summary$strain <- c(rep("CRISPR", 7), rep("BIM", 7)) %>% 
  as.factor

write.csv(x = selection_rate_summary, file = "./summary_data/selection_rate_summary.csv",
          row.names = F)

#### ---- Summary figure (facet version)---- ####
# selection_rate_summary$term %<>% relevel(ref="24-clone\n(ancestral phage)")
# selection_rate_summary$term %<>% relevel(ref="24-clone")
# selection_rate_summary$term %<>% relevel(ref="12-clone")
# selection_rate_summary$term %<>% relevel(ref="6-clone")
# selection_rate_summary$term %<>% relevel(ref="3-clone")
# selection_rate_summary$term %<>% relevel(ref="1-clone\n(ancestral phage)")
# selection_rate_summary$term %<>% relevel(ref="1-clone")
# selection_rate_summary$strain %<>% relevel(ref="CRISPR")
# 
# strain_names <- list(
#   "CRISPR" = "All CRISPR clones",
#   "BIM" = "Susceptible clone"
# )
# 
# strain_labeller <- function(variable, value) {
#   return(strain_names[value])
# }
# 
# p1 <- ggplot(aes(x=term, y=mean), data=selection_rate_summary)+
#   labs(x="CRISPR allele diversity", y="Selection rate")+
#   geom_hline(yintercept = 0, linetype=2, colour="black")+
#   geom_errorbar(aes(ymin=l.67, ymax=h.67), width=0, size=4, alpha=0.5)+
#   geom_errorbar(aes(ymin=l.89, ymax=h.89), width=0, size=2, alpha=0.5)+
#   geom_errorbar(aes(ymin=l.95, ymax=h.95), width=0)+
#   geom_point(aes(position=term), pch=21, fill="white", 
#              colour="black", size=2)+
#   facet_wrap(~strain, nrow=2, labeller = strain_labeller)+
#   cowplot::theme_cowplot()+
#   theme(axis.text = element_text(size=12),
#         axis.title = element_text(face="bold", size=16),
#         strip.text = element_text(face="bold"),
#         strip.background = element_rect(size=11))+
#   NULL
# 
# p1
# 
# # Make the strip rectangles wider
# library(grid)
# g <- ggplotGrob(p1)
# g$heights
# pos =  c(subset(g$layout, grepl("panel", g$layout$name), select = t))
# for(i in pos) g$heights[i-1] = unit(0.5,"cm")
# grobs = which(grepl("strip", g$layout$name))
# for(i in grobs) g$grobs[[i]]$heights <-  unit(1, "npc")      
# grid.newpage()
# grid.draw(g)

# ggsave("new rate plot.png", g, path="./figs", device="png",
#        dpi=600, width=28, height=25, units=c("cm"))


#### ---- Summary figure (plot_grid version) ---- ####
raw_CRISPR <- d %>% 
  mutate(Treatment = fct_recode(Treatment, "1-clone\n(ancestral phage)" = "1-clone_control",
                                "24-clone\n(ancestral phage)" = "24-clone_control")) %>% 
  filter(Timepoint!="0") %>% 
  filter(Strain=="rCRISPR")
raw_BIM <- d %>% 
  mutate(Treatment = fct_recode(Treatment, "1-clone\n(ancestral phage)" = "1-clone_control",
                                "24-clone\n(ancestral phage)" = "24-clone_control")) %>% 
  filter(Timepoint!="0") %>% 
  filter(Strain=="rBIM")

CRISPR.coefs <- read.csv("./summary_data/selection_rate_summary.csv") %>% 
  filter(strain=="CRISPR")
BIM.coefs <- read.csv("./summary_data/selection_rate_summary.csv") %>% 
  filter(strain=="BIM") %>% 
  filter(term%in%c("3-clone", "6-clone", "12-clone", "24-clone", "24-clone\n(ancestral phage)"))

CRISPR.coefs$term %<>% relevel(ref="24-clone\n(ancestral phage)")
CRISPR.coefs$term %<>% relevel(ref="24-clone")
CRISPR.coefs$term %<>% relevel(ref="12-clone")
CRISPR.coefs$term %<>% relevel(ref="6-clone")
CRISPR.coefs$term %<>% relevel(ref="3-clone")
CRISPR.coefs$term %<>% relevel(ref="1-clone\n(ancestral phage)")
CRISPR.coefs$term %<>% relevel(ref="1-clone")

BIM.coefs$term %<>% relevel(ref="24-clone\n(ancestral phage)")
BIM.coefs$term %<>% relevel(ref="24-clone")
BIM.coefs$term %<>% relevel(ref="12-clone")
BIM.coefs$term %<>% relevel(ref="6-clone")
BIM.coefs$term %<>% relevel(ref="3-clone")

CRISPR.plot <- ggplot(aes(x=term, y=mean), data=CRISPR.coefs)+
  labs(x="CRISPR diversity", y="Selection rate of\nall CRISPR clones")+
  geom_hline(yintercept = 0, linetype=2, colour="black")+
  # geom_jitter(aes(y=r, x=Treatment), data=raw_CRISPR,
  #             width=0.2, pch=21, alpha=0.5, size=0.8)+
  geom_errorbar(aes(ymin=l.67, ymax=h.67), width=0, size=3, alpha=0.5)+
  geom_errorbar(aes(ymin=l.89, ymax=h.89), width=0, size=1.5, alpha=0.5)+
  geom_errorbar(aes(ymin=l.95, ymax=h.95), width=0)+
  geom_point(aes(position=term), pch=21, fill="white", 
             colour="black", size=1.5)+
  cowplot::theme_cowplot()+
  theme(axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=10),
        strip.text = element_text(face="bold", size=10),
        strip.background = element_rect(size=11))+
  coord_cartesian(ylim=c(-2.612728, 1.235529))+
  NULL

last_plot()

BIM.plot <- ggplot(aes(x=term, y=mean), data=BIM.coefs)+
  labs(x="CRISPR diversity", y="Selection rate of\nsusceptible CRISPR clones")+
  geom_hline(yintercept = 0, linetype=2, colour="black")+
   # geom_jitter(aes(y=r, x=Treatment), data=raw_BIM,
   #             width=0.2, pch=21, alpha=0.5, size=0.8)+
  geom_errorbar(aes(ymin=l.67, ymax=h.67), width=0, size=3, alpha=0.5)+
  geom_errorbar(aes(ymin=l.89, ymax=h.89), width=0, size=1.5, alpha=0.5)+
  geom_errorbar(aes(ymin=l.95, ymax=h.95), width=0)+
  geom_point(aes(position=term), pch=21, fill="white", 
             colour="black", size=1.5)+
  cowplot::theme_cowplot()+
  theme(axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=10),
        strip.text = element_text(face="bold"),
        strip.background = element_rect(size=11))+
  coord_cartesian(ylim=c(-2.612728, 1.235529))+
  NULL
last_plot()

# Fig_S3 <- plot_grid(CRISPR.plot+xlab(""),
#                    BIM.plot, nrow=2, rel_heights = c(1,1),
#                    labels = c("A", "B"), label_size = 10)
last_plot()

# ggsave("Figure_S3.tif", Fig_S3, path="./figs/", device="tiff", compression="lzw",
#                dpi=300, width=12, height=12, units=c("cm"))



#### --- Models with Time as a continuous variable --- ####
# A reviewer after the initial submission suggested that Time (that is, dpi) be
# treated as a continuous variable in order to provide a clearer indication of
# trends by presenting a directional regression slope, rather than a non-directional
# F-statistic

## Data
d.cont <- read.csv("./original_data/selection rates.csv", header=T)  %>% 
  select(Treatment, Timepoint, Replicate, BIM_mix, Escape_phage, Tracked_BIM,
         rCRISPR, rBIM) %>% 
  gather("rCRISPR", "rBIM", key="Strain", value="r", factor_key = T)

d.cont$Replicate %<>% as.factor
d.cont$Tracked_BIM %<>% as.factor()
d.cont$Escape_phage %<>% as.factor

# Relevel the datasets so things make sense
d$Treatment %<>% relevel(ref="24-clone_control")
d$Treatment %<>% relevel(ref="24-clone")
d$Treatment %<>% relevel(ref="12-clone")
d$Treatment %<>% relevel(ref="6-clone")
d$Treatment %<>% relevel(ref="3-clone")
d$Treatment %<>% relevel(ref="1-clone_control")
d$Treatment %<>% relevel(ref="1-clone")

# CRISPR model 
# d.cont.CRISPR <- d.cont %>% 
#   na.exclude %>% 
#   filter(Strain=="rCRISPR", Timepoint!="0")

# Use this dataframe to analyse just the polyclonal treatments
d.cont.CRISPR <- d.cont %>%
  na.exclude %>%
  filter(Strain=="rCRISPR", Timepoint!="0",
         Treatment %in% c("3-clone", "6-clone", "12-clone", "24-clone", "24-clone_control"))

d.cont.CRISPR$Treatment %<>% relevel(ref="24-clone_control")
d.cont.CRISPR$Treatment %<>% relevel(ref="1-clone_control")
d.cont.CRISPR$Treatment %<>% relevel(ref="24-clone")
d.cont.CRISPR$Treatment %<>% relevel(ref="12-clone")
d.cont.CRISPR$Treatment %<>% relevel(ref="6-clone")
d.cont.CRISPR$Treatment %<>% relevel(ref="3-clone")
d.cont.CRISPR$Treatment %<>% relevel(ref="1-clone")

m1 <- lmer(r~Treatment + (1|Replicate), data=d.cont.CRISPR)
m2 <- lmer(r~Timepoint + (1|Replicate), data=d.cont.CRISPR)
m3 <- lmer(r~Treatment + Timepoint +(1|Replicate), data=d.cont.CRISPR)
m4 <- lmer(r~Treatment * Timepoint + (1|Replicate), data=d.cont.CRISPR)

AIC(m1, m2, m3, m4) %>% 
  compare_AICs()

plot(m1)
plot(m2)
plot(m3)
plot(m4)

anova(m3, type="marginal")
summary(m2)
confint(m2, parm="beta_", method="boot")

# BIM model 
# d.cont.BIM <- d.cont %>% 
#   na.exclude %>% 
#   filter(Strain=="rBIM", Timepoint!="0")

# Use this dataframe to analyse just the polyclonal treatments
d.cont.BIM <- d.cont %>%
  na.exclude %>%
  filter(Strain=="rBIM", Timepoint!="0",
         Treatment %in% c("3-clone", "6-clone", "12-clone", "24-clone", "24-clone_control"))

d.cont.BIM$Treatment %<>% relevel(ref="24-clone_control")
d.cont.BIM$Treatment %<>% relevel(ref="1-clone_control")
d.cont.BIM$Treatment %<>% relevel(ref="24-clone")
d.cont.BIM$Treatment %<>% relevel(ref="12-clone")
d.cont.BIM$Treatment %<>% relevel(ref="6-clone")
d.cont.BIM$Treatment %<>% relevel(ref="3-clone")
d.cont.BIM$Treatment %<>% relevel(ref="1-clone")

m1 <- lmer(r~Treatment + (1|Replicate), data=d.cont.BIM)
m2 <- lmer(r~Timepoint + (1|Replicate), data=d.cont.BIM)
m3 <- lmer(r~Treatment + Timepoint +(1|Replicate), data=d.cont.BIM)
m4 <- lmer(r~Treatment * Timepoint + (1|Replicate), data=d.cont.BIM)

AIC(m1, m2, m3, m4) %>% 
  compare_AICs()

plot(m1)
plot(m2)
plot(m3)
plot(m4)

anova(m3, type="marginal")
summary(m3)
confint(m3, parm="beta_", method="boot")

