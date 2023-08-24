# Header ---------------------------
#
# Script name: PPFanalysis.R
#
# Purpose of script: Statistical analysis and plot creation for PPF paper
#
# Author: Ethan N. Okoshi
#
# Date Created: 2023-03-01
#
# Copyright (c) Ethan Okoshi, 2023
# Email: ethanokoshi@gmail.com
#
# Notes:
#
# set working directory
#
# setwd()
#
# set options
#
options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
Sys.setenv(LANG = "EN")
 extrafont::loadfonts(device = "win")
#
# load functions from another file
#
# source()
#
# load data
#
input_sheet_path <- "input/PPF for USCAP230309_R.xlsx"
library(readxl)
sheet <- read_excel(input_sheet_path, sheet = "Sheet3")
sheet[c("Dlco", "%Dlco")] <- sapply(sheet[c("Dlco", "%Dlco")], as.numeric)
#
# load packages ---------------
require(tidyverse)
library(magrittr)
library(broom)
library(survival)
library(survminer)
library(ggsurvfit)
library(ggExtra)
library(geomtextpath)
library(gridExtra)
library(grid)
library(glue)
library(gt)
library(gtsummary)
library(labelled)


# fig 3 survial plot path UIP ----
png("output/pathUIPsurv_v3.png", res = 300, width = 6, height = 4, units = "in")
survfit2(Surv(obs_months, death_or_lung_transplant) ~ pathUIP, data = sheet) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm"))
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Pathological UIP Pattern -", "Pathological UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate(
    "text",
    x = 90,
    y = 0.94,
    label = "Pathological UIP -",
    color = "#e6ac00",
    family = "Barlow Medium",
    hjust = 0
  ) +
  annotate(
    "text",
    x = 70,
    y = 0.38,
    label = "Pathological UIP +",
    color = "#286d8f",
    family = "Barlow Medium",
    hjust = 0
  ) +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = "log-rank test, p = 0.0005",
    color = "grey30",
    family = "Barlow Medium",
    hjust = 0,
    vjust = 0
  )
dev.off()

# fig 4 suvival plot 2 focaluip individually --------
png("output/focUIPsurv_consensus_v1.png", res = 600, width = 6, height = 4, units = "in")
survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP_consensus, data = sheet) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm"))
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 90, y = 0.97, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 65, y = 0.43, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p < 0.0001", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)
dev.off()
# fig 4 suvival plot 2 focaluip by
png("output/focUIPsurv_patho1_v1.png", res = 600, width = 6, height = 4, units = "in")
survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP1, data = sheet) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm"))
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 60, y = 0.92, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 65, y = 0.41, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p < 0.0001", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)
dev.off()

png("output/focUIPsurv_patho2_v1.png", res = 600, width = 6, height = 4, units = "in")
survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP2, data = sheet) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm"))
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 85, y = 0.97, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 65, y = 0.41, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p < 0.0001", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)
dev.off()

png("output/focUIPsurv_patho3_v1.png", res = 600, width = 6, height = 4, units = "in")
survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP3, data = sheet) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm"))
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 85, y = 0.94, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 65, y = 0.46, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p = 0.0002", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)
dev.off()

# fig 4 gridarrange ------
p1 <- 
  survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP1, data = sheet) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm"))
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 60, y = 0.92, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 65, y = 0.41, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p < 0.0001", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)

p2 <- 
  survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP2, data = sheet) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm"))
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 85, y = 0.97, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 65, y = 0.41, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p < 0.0001", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)

p3 <- 
  survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP3, data = sheet) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm"))
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 85, y = 0.94, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 65, y = 0.46, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p = 0.0002", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)

pc <- survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP_consensus, data = sheet) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm"))
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 90, y = 0.97, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 80, y = 0.43, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p < 0.0001", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)

png("output/fig4_sharedaxistitles.png", res = 300, width = 20, height = 20, units = "cm")
grid.arrange(
  arrangeGrob(p1 + labs(title = "Pathologist 1") + theme(axis.title.x.bottom = element_blank()),
    p2 + labs(title = "Pathologist 2") + theme(axis.title.x.bottom = element_blank()),
    p3 + labs(title = "Pathologist 3"),
    ncol = 1,
    top = textGrob("A", gp = gpar(fontfamily = "Barlow SemiBold", fontsize = 24), x = 0.1)
  ),
  arrangeGrob(pc + labs(title = "Consensus"),
    top = textGrob("B", gp = gpar(fontfamily = "Barlow SemiBold", fontsize = 24), x = 0.1)
  ),
  layout_matrix = cbind(c(1, 1, 1), c(1, 1, 1), c(2, 2, 2), c(2, 2, 2), c(2, 2, 2)),
  padding = unit(0.5, "cm")
)
dev.off()
# fig 5 focalUIP by etiology ------

p1 <- survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP_consensus, data = filter(sheet, disease == "cHP")) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm")),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate", title = "cHP (n=14)") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 75, y = 0.94, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 60, y = 0.43, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p = 0.107", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)

p2 <- survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP_consensus, data = filter(sheet, disease == "CTD-ILD")) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm")),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate", title = "CTD-ILD (n=82)") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 102, y = 0.92, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 55, y = 0.58, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p = 0.0375", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)

p3 <- survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP_consensus, data = filter(sheet, disease == "UC-ILD")) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm")),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate", title = "UC-ILD (n=82)") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 75, y = 0.94, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 80, y = 0.43, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p = 0.0075", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)

p4 <- survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP_consensus, data = filter(sheet, disease == "iNSIP")) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm")),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate", title = "iNSIP (n=19)") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0005"
  ) +
  annotate("text", x = 85, y = 0.92, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 70, y = 0.45, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p = 0.0131", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)

png("output/fig5_focalUIPbyetiology_v1.png", res = 300, width = 12, height = 8, units = "in")
grid.arrange(p1,
  p2,
  p3,
  p4,
  padding = unit(0.5, "cm")
)
dev.off()

# fig 6 focalUIP within pathUIP- cases ------
png("output/focalUIPinpathUIPneg_v1.png", res = 300, width = 6, height = 4, units = "in")
survfit2(Surv(obs_months, death_or_lung_transplant) ~ focalUIP_consensus, data = filter(sheet, pathUIP == 0)) %>%
  ggsurvfit(linewidth = 1, show.legend = FALSE) +
  add_censor_mark(shape = "|", size = 3, alpha = 0.7, show.legend = FALSE) +
  add_confidence_interval(alpha = 0.15, show.legend = FALSE) +
  xlim(0, 125) +
  ylim(0, 1) +
  theme(
    text = element_text(family = "Barlow Medium"),
    panel.grid.major = element_line(colour = "gray85", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray85", linetype = "dotted"),
    axis.title.y = element_text(margin = margin(0, 0.2, 0, 0, "cm")),
    axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, "cm"))
  ) +
  labs(x = "Months after Diagnosis", y = "Survival Rate") +
  scale_fill_manual(values = c("#e6ac00", "#286d8f")) +
  scale_color_manual(
    labels = c("Focal UIP Pattern -", "Focal UIP Pattern +"),
    values = c("#e6ac00", "#286d8f"),
    name = "log-rank test, p = 0.0018"
  ) +
  annotate("text", x = 90, y = 0.97, label = "Focal UIP -", color = "#e6ac00", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 75, y = 0.55, label = "Focal UIP +", color = "#286d8f", family = "Barlow Medium", hjust = 0) +
  annotate("text", x = 0, y = 0, label = "log-rank test, p = 0.0018", color = "grey30", family = "Barlow Medium", hjust = 0, vjust = 0)
dev.off()

# cox univariate testing -------------
variables <- c("focalUIP_consensus", "disease", "Age", "Male", "nonsmoker", "FVC", "`%FVC`", "Dlco", "`%Dlco`")

univariate_cox_results <-
  tibble(variables) %>%
  mutate(formula = map(variables, \(x) formula(glue("Surv(obs_months, death_or_lung_transplant) ~ {x}")))) %>%
  mutate(coxresults = map(formula, \(x) coxph(x, data = sheet))) %>%
  mutate(summary = map(coxresults, tidy, exponentiate = TRUE, conf.int = TRUE)) %>%
  unnest_wider(summary) %>%
  select(term, estimate, p.value, conf.low, conf.high) %>%
  separate_rows(c(term, estimate, p.value, conf.low, conf.high), sep = ", ") %>%
  mutate(
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3)
  ) %>%
  mutate(
    conf.low = as.character(conf.low),
    conf.high = as.character(conf.high)
  ) %>%
  unite(conf.low, conf.high, col = "CI", sep = ", ")

write_csv(univariate_cox_results, "output/univariate_cox_results.csv")

# cox multivariate testing ------------------

# this one uses UC-ILD as the "reference disease"
res.cox <- coxph(Surv(obs_months, death_or_lung_transplant) ~ focalUIP_consensus + disease + Age + Male + smoking_hist + FVC + `%FVC` + Dlco + `%Dlco`, data = sheet)
summary(res.cox)


multivariate_cox_results <- res.cox %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3)
  ) %>%
  mutate(
    conf.low = as.character(conf.low),
    conf.high = as.character(conf.high)
  ) %>%
  unite(conf.low, conf.high, col = "CI", sep = ", ") %>%
  mutate(
    std.error = NULL,
    statistic = NULL
  )

write_csv(multivariate_cox_results, "output/multivariate_cox_results.csv")


test.ph <- cox.zph(res.cox)
test.ph

ggcoxdiagnostics(res.cox,
  type = "dfbeta",
  linear.predictions = FALSE, ggtheme = theme_bw()
)

ggcoxfunctional(Surv(obs_months, death_or_lung_transplant) ~ Age + log(Age) + sqrt(Age), data = sheet)


# cox permutation testing ---------------

n_perm <- 2000

cox_perm <- function() {
  sheet.iter <- filter(sheet, disease != "iPPFE")
  sheet.iter$permdisease <- sample(sheet.iter$disease, replace = FALSE)
  coxph.results <- sheet.iter %>%
    group_by(permdisease) %>%
    tidyr::nest() %>%
    mutate(data = map(data, ~ coxph(Surv(obs_months, death_or_lung_transplant) ~ focalUIP_consensus, data = .x))) %>%
    mutate(data = map(data, \(x) data.frame(summary(x)$coef))) %>%
    unnest(data)
  colnames(coxph.results)[3] <- "HR"
  colnames(coxph.results)[6] <- "p.value"
  return(coxph.results)
}

coxph_results <- data.frame(1:n_perm) %>%
  mutate(perm = map(1:n_perm, \(i) cox_perm(), .progress = TRUE)) %>%
  unnest(perm)

coxph_real <- filter(sheet, disease != "iPPFE") %>%
  group_by(disease) %>%
  tidyr::nest() %>%
  mutate(data = map(data, ~ coxph(Surv(obs_months, death_or_lung_transplant) ~ focalUIP_consensus, data = .x))) %>%
  mutate(data = map(data, \(x) data.frame(summary(x)$coef))) %>%
  unnest(data)
colnames(coxph_real)[3] <- "HR"
colnames(coxph_real)[6] <- "p.value"

# fig 7 cox permutation test -------------

create_marginal_plot_func <- function(coxph_results, permdisease_str, n) {
  if (coxph_real[coxph_real["disease"] == permdisease_str, ][["HR"]] > 1000) {
    vlabel <- paste("HR =", format(coxph_real[coxph_real["disease"] == permdisease_str, ][["HR"]], scientific = TRUE, digits = 4))
  } # this statement checks if the HR needs scientific formatting
  else {
    vlabel <- paste("HR =", format(coxph_real[coxph_real["disease"] == permdisease_str, ][["HR"]], scientific = FALSE, digits = 4))
  } # sets the vertical line label
  coxph_results %>%
    filter(permdisease == permdisease_str) %>%
    mutate(res.ecdf = ecdf(.$coef)(.$coef)) %>%
    filter(.$res.ecdf < 0.99 & .$res.ecdf > 1 - 0.99 & .$p.value < 0.99, ) %>%
    {ggMarginal(
      ggplot(aes(x = HR, y = -log10(p.value)), data = .) +
        geom_point(alpha = 0.2, stroke = NA, size = 1.5) +
        labs(caption = glue("{permdisease_str} (n={n})")) +
        xlab("Hazard Ratio") +
        ylab(expression(paste(-log[10]("p-value")))) +
        scale_x_log10() +
        geom_labelhline(aes(
          yintercept = -log10(coxph_real[coxph_real["disease"] == permdisease_str, ][["p.value"]]),
          label = paste("p =", round(coxph_real[coxph_real["disease"] == permdisease_str, ][["p.value"]], 5)),
        )) +
        geom_labelvline(aes(
          xintercept = (coxph_real[coxph_real["disease"] == permdisease_str, ][["HR"]]),
          label = vlabel,
          hjust = 0.8)) +
        ggpubr::theme_pubr(base_family = "Barlow Medium") +
        theme(
          text = element_text(family = "Barlow Medium"),
          panel.grid.major = element_line(colour = "gray70", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "gray70", linetype = "dotted"),
        ),
      type = "densigram"
    )} #the brackets here are so that the data doesnt get piped into the first argument spot, rather where the . is
}

p1 <- create_marginal_plot_func(coxph_results, "cHP", 14)

p2 <- create_marginal_plot_func(coxph_results, "CTD-ILD", 82)

p3 <- create_marginal_plot_func(coxph_results, "UC-ILD", 84)

p4 <- create_marginal_plot_func(coxph_results, "iNSIP", 19)

png("output/fig7_permutationtestcox.png", res = 500, width = 12, height = 8, units = "in")
grid.arrange(p1, p2, p3, p4)
dev.off()

# table 1 population statistics ---------------

var_label(sheet) <- list(
  focalUIP_consensus = "Focal UIP",
  disease = "Disease",
  obs_months = "Time to Event",
  death_or_lung_transplant = "Event",
  smoking_hist = "Smoking history",
  Male = "Sex"
)

t1 <- tbl_merge(
  tbls = list(sheet %>% #pathUIP
                select(pathUIP, Male, Age, smoking_hist, FVC, `%FVC`, Dlco, `%Dlco`, `KL-6`, disease) %>%
                mutate(pathUIP = if_else(pathUIP == 1, "Positive", "Negative")) %>%
                tbl_summary(by = pathUIP,
                            digits = Age ~ 2,
                            statistic = list(all_continuous() ~ "{mean} ± {sd}"),
                            missing = "ifany",
                            missing_text = "Missing",
                            percent = "row") %>%
                add_p() %>%
                bold_p() %>%
                separate_p_footnotes() %>%
                modify_table_styling(columns = all_of("p.value"), 
                                     text_format = "italic") %>%
                modify_header(p.value = "***p-value***"),
              sheet %>% #focalUIP
                select(focalUIP_consensus, Male, Age, smoking_hist, FVC, `%FVC`, Dlco, `%Dlco`, `KL-6`, disease) %>%
                mutate(focalUIP_consensus = if_else(focalUIP_consensus == 1, "Positive", "Negative")) %>%
                tbl_summary(by = focalUIP_consensus,
                            digits = Age ~ 2,
                            statistic = list(all_continuous() ~ "{mean} ± {sd}"),
                            missing = "ifany",
                            missing_text = "Missing",
                            percent = "row") %>%
                add_p() %>%
                bold_p() %>%
                separate_p_footnotes() %>%
                modify_table_styling(columns = all_of("p.value"), 
                                     text_format = "italic") %>%
                modify_header(p.value = "***p-value***"),
              sheet %>% #Total
                select(Male, Age, smoking_hist, FVC, `%FVC`, Dlco, `%Dlco`, `KL-6`, disease) %>%
                tbl_summary(
                  digits = Age ~ 2,
                  statistic = list(all_continuous() ~ "{mean} ± {sd}"),
                  missing = "ifany",
                  missing_text = "Missing",
                  percent = "row")
  ),
  tab_spanner = c("**Pathological UIP**", "**Focal UIP**", "**Total**")
) %>% as_gt()

# was having a glitch where multiple footnotes were showing up for one p-value
t1[["_footnotes"]] %<>% filter(!(colname == "p.value_1" & rownum == 5 & footnotes == "Fisher's exact test"),
                               !(colname == "p.value_2" & rownum == 5 & footnotes == "Pearson's Chi-squared test"))
gtsave(t1, "output/table1.png", zoom = 10, delay = 0.5, expand = 20)


# table 2 disease stats ---------------------

table2 <- tribble(~disease,~n, ~def, ~prob, ~poss, ~other, ~rate,
                  "iNSIP",19,0,0,0,4,"21%",
                  "cHP",14,5,1,2,1,"64%",
                  "CTD-ILD",82,8,16,3,8,"43%",
                  "iPPFE",2,0,0,0,1,"50%",
                  "UC-ILD",84,13,16,11,8,"57%")

gt(table2, rowname_col = "disease") %>%
  tab_stubhead(label = md("**Etiology**")) %>%
  tab_header(title = md("**PPF**; Progressive Pulmonary Fibrosis (n = 201)")) %>%
  cols_label(
    n = "N",
    def = "Definite UIP",
    prob = "Probable UIP",
    poss = "Possible UIP",
    other = "UIP vs. Other",
    rate = "Percent UIP+"
  ) %>%
  tab_options(
    heading.title.font.size = "medium",
    table.width = px(720),
    table.align = "left"
  ) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body()
  ) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels()
  ) %>%
  gtsave("output/table2.png", zoom = 10, delay = 0.5)

# table 3 cox results in gt --------------

var_label(sheet) <- list(
  focalUIP_consensus = "Focal UIP",
  disease = "Disease",
  obs_months = "Time to Event",
  death_or_lung_transplant = "Event",
  smoking_hist = "Smoking history",
  Male = "Sex"
)

multivariate_tbl <-
  tbl_regression(
    res.cox,
    exponentiate = TRUE,
    show_single_row = "Male"
  ) %>%
  bold_p()

univariate_tbl <- sheet %>%
  select(focalUIP_consensus, obs_months, death_or_lung_transplant, disease, FVC, `%FVC`, Dlco, `%Dlco`, smoking_hist, Male, Age) %>%
  tbl_uvregression(
    method = coxph,
    y = Surv(obs_months, death_or_lung_transplant),
    exponentiate = TRUE,
    show_single_row = "Male"
  ) %>%
  bold_p()

table3 <-
  tbl_merge(
    tbls = list(univariate_tbl, multivariate_tbl),
    tab_spanner = c("**Univariate**", "**Multivariate**")
  ) %>%
  as_gt() %>%
  text_replace(
    locations = cells_body(columns = c("ci_1", "ci_2")),
    pattern = "(^\\d.*)",
    replacement = ("\\[\\1\\]")
  ) %>%
  gtsave("output/table3.png", zoom = 10)
