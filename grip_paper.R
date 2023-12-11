library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(psych)

options(scipen = 999)

path <- "~/Documents/PhD/Grip/"
path_metal <- "/Users/filippoabbondanza/Documents/PhD/Grip/"
filename <- "grip.sample"

### GS measures for ALSPAC ####

# Nomenclature
# grip.best === GSD
# grip.worst === GSND

# This has the zscore transformed values
grip_zscore <- read.table(paste0(path,'/grip.20022022.zscore.sample'), header = T) %>%
  mutate(sex = case_when(
    sex == '1' ~ 'Male',
    sex == '2' ~ 'Female'
  )) %>%
  mutate(sex = as.factor(sex))

# This has the raw values
grip_raw <- read.table(paste0(path,'/grip.20022022.raw.sample'), header = T) %>%
  mutate(sex = case_when(
    sex == '1' ~ 'Male',
    sex == '2' ~ 'Female'
  )) %>%
  mutate(sex = as.factor(sex))


# Check grip strength by sex and hand
grip_raw %>% group_by(sex) %>% summarise(across(.cols = c("grip.best.hand", "grip.worst.hand"), ~mean(., na.rm=T)))

psych::describe(grip_raw$grip.best.hand)
psych::describe(grip_raw$grip.best.hand[grip$sex=="Male"])
psych::describe(grip_raw$grip.best.hand[grip$sex=="Female"])

psych::describe(grip_raw$grip.worst.hand)
psych::describe(grip_raw$grip.worst.hand[grip$sex=="Male"])
psych::describe(grip_raw$grip.worst.hand[grip$sex=="Female"])


# Get distribution for GS traits ####
grip.best <- ggplot(data = grip_raw) + 
  geom_histogram(aes(grip.best.hand, fill = sex, alpha = 0.5), position = 'identity') +
  xlab('GSD (kg)') + 
  scale_fill_manual(values = c("blue", "red")) +
  scale_alpha(guide = 'none')

grip.worst <- ggplot(data = grip_raw) + 
  geom_histogram(aes(grip.worst.hand, fill = sex, alpha = 0.5), position = 'identity') +
  xlab('GSND (kg)') + 
  scale_fill_manual(values = c("blue", "red")) +
  scale_alpha(guide = 'none')

plots <- ggarrange(grip.best, grip.worst,
                   labels = c("A", "B"),
                   ncol = 3, nrow = 1, common.legend = T, legend = "right")

annotate_figure(plots, top = text_grob("Distribution of grip strength by sex", size = 14))

# Plots for age effect ####
GSD_age <- ggplot(grip_raw, aes(x = grip.best.hand, y = age)) +
  geom_point() +
  xlab('GSD') +
  geom_smooth(method = "lm") +
  # ggtitle('Age effect on GSR') +
  stat_poly_eq(formula = ' y ~ x',
               # eq.with.lhs = "italic(y)~`=`~",
               aes(label = paste(..p.value.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)

GSND_age <- ggplot(grip_raw, aes(x = grip.worst, y = age)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab('GSL') +
  # ggtitle('Age effect on GSL') +
  stat_poly_eq(formula = ' y ~ x',
               # eq.with.lhs = "italic(y)~`=`~",
               aes(label = paste(..p.value.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)

plot_age <- ggarrange(GSR_age, GSL_age,
                      labels = c("A", "B"))
plot_age



# Sex effects ####
t.test(grip_raw$grip.best, as.numeric(grip_raw$sex))$p
pvalue <- '2.2e-16' # Will be used for the plots

means <- aggregate(grip.right ~sex, grip, mean)
GSD_boxplot <- ggplot(grip, aes(y = grip.right, x = sex)) +
  geom_boxplot() +
  ylab('GSR') +
  stat_summary(fun=mean, colour="darkred", geom="point", 
               shape=18, size=3, show.legend=FALSE) +
  geom_text(data = means, aes(label = paste0('Mean: ', round(grip.right,2)), y = grip.right + 1))

means <- aggregate(grip.left ~sex, grip, mean)
GSND_boxplot <- ggplot(grip, aes(y = grip.left, x = sex)) +
  geom_boxplot() + 
  ylab('GSL') +
  stat_summary(fun=mean, colour="darkred", geom="point", 
               shape=18, size=3, show.legend=FALSE) +
  geom_text(data = means, aes(label = paste0('Mean: ', round(grip.left,2)), y = grip.left + 1))

boxplots <- ggarrange(GSR_boxplot, GSL_boxplot, 
                      labels = c("A"),
                      ncol = 3, nrow = 1)



# Make plots for genetic correlations ####
rg_to_plot = c("ADHD", "BIP", "CAD", "Overall fracture risk",
               "HBD.left", "HBD.right", "Heart attack", "ASD",
               "SCZ", "AD")


gen.corr <- read.table(paste0(path_metal,'genetic_correlation_best_worst_METAL_reanalysis_CAD.tsv'), sep = '\t', header = T) %>%
  mutate(
    rgplus = rg+se,
    rgminus = rg-se,
    Trait.2 = as.factor(Trait.2),
    Trait.1 = as.factor(Trait.1)
  ) %>%
  filter(Trait.2 %in% rg_to_plot) %>%
  mutate(Significant = case_when(
    p <= 0.0036 ~ 'Significant (p<0.0036)',
    p <= 0.05 ~ 'Nominally significant (p < 0.05)',
    TRUE ~ 'Non significant',
  )) %>%
  mutate(Trait.1.ordered = factor(Trait.1, levels = c("GSD","GSND"))) %>%
  group_by(Trait.1) %>%
  mutate(Trait.2.ordered = fct_reorder(Trait.2, rg, .fun='mean'))

ggplot(gen.corr, aes(x = rg, y = Trait.2.ordered,  xmin = rgminus, xmax = rgplus, label = round(rg,2), color=Significant)) +
  geom_point() +
  theme_minimal() +
  theme_bw() +
  geom_vline(xintercept = 0.0) +
  ylab(label = '') +
  xlab('rg (+- SE)') +
  geom_errorbarh(height=0.2) +
  geom_text(nudge_y = 0.3) +
  facet_wrap(~Trait.1)


# Get overlap of loci for replication study ####
gsb_fuma <- read.table("~/Documents/PhD/Grip/FUMA_job171616_GSB.metal.MAF0.05.common.reanalysis/snps.txt", h=T) %>% rename(MarkerName = "rsID") %>% select(MarkerName)
gsb_METAL <- data.table::fread("~/Documents/PhD/Grip/GS.best.MAF0.05.reanalysis.common.final.rsid.gz")
gsb_join <- inner_join(gsb_fuma, gsb_METAL)
nrow(gsb_join)

gsw_fuma <- read.table("~/Documents/PhD_data/Grip_GWAS/Analysis_GSBest_GSWorst/METAL_BOLT/Replication/FUMA_job162935_GSW/snps.txt", h=T) %>% rename(MarkerName = "rsID") %>% select(MarkerName)
gsw_METAL <- data.table::fread("~/Documents/PhD/Grip/GS.worst.MAF0.05.reanalysis.common.final.rsid.gz")
gsw_join <- inner_join(gsw_fuma, gsw_METAL)
nrow(gsw_join)


gs_join <- inner_join(gsb_METAL, gsw_METAL, by = "SNPID")
