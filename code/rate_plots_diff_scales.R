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

#### ---- CRISPR boxplot ---- ####
CRISPR_left_raw <- raw_CRISPR %>% 
  filter(Treatment%in%c("1-clone", "1-clone\n(ancestral phage)"))

CRISPR_right_raw <- raw_CRISPR %>% 
  filter(Treatment%in%c("3-clone", "6-clone", "12-clone", 
                        "24-clone","24-clone\n(ancestral phage)"))

pd <- position_dodge(0.6)
library(wesanderson)
pal <- wes_palette("IsleofDogs1", 3)

CRISPR_left<- ggplot(aes(x=Treatment, y=r), data=CRISPR_left_raw)+
  geom_boxplot(width=0.5, outlier.shape = NA, aes(colour=Timepoint), position = pd)+
  geom_point(aes(colour=Timepoint), alpha=0.5, pch=21, size=0.8, position = pd)+
  geom_hline(yintercept=0, linetype=2, size=0.5)+
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
        legend.key.width = unit(0.6, "cm"),
        legend.key = element_rect(fill="grey95"))+
  NULL

last_plot()

CRISPR_right <- ggplot(aes(x=Treatment, y=r), data=CRISPR_right_raw)+
  geom_boxplot(width=0.5, outlier.shape = NA, aes(colour=Timepoint), position = pd)+
  geom_point(aes(colour=Timepoint), alpha=0.5, pch=21, size=0.8, position = pd)+
  geom_hline(yintercept=0, linetype=2, size=0.5)+
  #stat_summary(fun.y="mean", geom="point", colour="red")+
  #facet_wrap(~Timepoint, scales="free", labeller = timepoint_labeller, ncol=2)+
  labs(y="Selection rate of\nall CRISPR clones", x="CRISPR diversity")+
  scale_x_discrete(labels = treatment_list[3:8])+
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
        legend.key.width = unit(0.6, "cm"),
        legend.key = element_rect(fill="grey95"))+
  NULL

last_plot()

CRISPR_test_top <- plot_grid(CRISPR_left+theme(legend.position = "none"), CRISPR_right+ylab(""),
                         rel_widths = c(0.5,1), align="h",
                         labels = c("A", "B"), label_size = 10)

CRISPR_test <- plot_grid(CRISPR_test_top, CRISPR.plot,
                         nrow=2, rel_widths = c(1.7,1),
                         labels=c("", "C"), label_size = 10)
last_plot()

ggsave("Figure_3.tif", CRISPR_test, path = "./figs/", device="tiff",
       dpi=300, width=21, height=15, units=c("cm"))

#### ---- BIM boxplot ---- ####
BIM_left_raw <- raw_BIM %>% 
  filter(Treatment%in%c("3-clone", "6-clone"))

BIM_right_raw <- raw_BIM %>% 
  filter(Treatment%in%c("12-clone", 
                        "24-clone","24-clone\n(ancestral phage)"))

BIM_left<- ggplot(aes(x=Treatment, y=r), data=BIM_left_raw)+
  geom_boxplot(width=0.5, outlier.shape = NA, aes(colour=Timepoint), position = pd)+
  geom_point(aes(colour=Timepoint), alpha=0.5, pch=21, size=0.8, position = pd)+
  geom_hline(yintercept=0, linetype=2, size=1)+
  #stat_summary(fun.y="mean", geom="point", colour="red")+
  #facet_wrap(~Timepoint, scales="free", labeller = timepoint_labeller, ncol=2)+
  labs(y="Selection rate of\nsusceptible CRISPR clones", x="CRISPR diversity")+
  scale_x_discrete(labels = treatment_list[3:4])+
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
        legend.key.width = unit(0.6, "cm"),
        legend.key = element_rect(fill="grey95"))+
  NULL

last_plot()

BIM_right <- ggplot(aes(x=Treatment, y=r), data=BIM_right_raw)+
  geom_boxplot(width=0.5, outlier.shape = NA, aes(colour=Timepoint), position = pd)+
  geom_point(aes(colour=Timepoint), alpha=0.5, pch=21, size=0.8, position = pd)+
  geom_hline(yintercept=0, linetype=2, size=1)+
  #stat_summary(fun.y="mean", geom="point", colour="red")+
  #facet_wrap(~Timepoint, scales="free", labeller = timepoint_labeller, ncol=2)+
  labs(y="Selection rate of\nsusceptible CRISPR clones", x="CRISPR diversity")+
  scale_x_discrete(labels = treatment_list[5:8])+
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
        legend.key.width = unit(0.6, "cm"),
        legend.key = element_rect(fill="grey95"))+
  NULL

last_plot()

BIM_test_top <- plot_grid(BIM_left+theme(legend.position = "none"), 
                      BIM_right+ylab(""),
                      rel_widths = c(0.6,1), align="h", 
                      label_size=10, labels = c("A", "B"))


BIM_test <- plot_grid(BIM_test_top, BIM.plot,
                         nrow=2, rel_widths = c(1.7,1),
                         labels=c("", "C"), label_size = 10)

ggsave("Figure_4.tif", BIM_test, path = "./figs/", device="tiff",
       dpi=300, width=21, height=15, units=c("cm"))

