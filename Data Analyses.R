# data analysis for predictive ability of single hills
# kayla altendorf
# last updated: 2/5/25

#### load packages ####
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
library(readxl)
library(stringr)
library(Matrix)
library(lme4)
library(lmerTest)
library(ggplot2)
library(emmeans)
library(multcomp)
library(robustbase)
library(ggResidpanel)
library(cowplot)
library(Hmisc)
library(gtools)

# set dir 
setwd("/users/kayla.altendorf/OneDrive - USDA/Documents/2024/Predictive Ability of Single Hills/")

#### read in data ####
spacing_trial_dat <- read.csv("./Data Analysis Results/spacing_trial_dat.csv", na.strings = c(".", "", "NA"))

#### correlation between observed and predicted dry cone yields #### 
# calculate predicted variable
dat22 <- spacing_trial_dat %>% filter(year == "2022") 

# pearson correlation test
cor.test(dat22$predicted_dry_cone_weight_g, dat22$dry_cone_weight_g) 
summary(lm(sqrt(dat22$predicted_dry_cone_weight_g) ~ sqrt(dat22$dry_cone_weight_g)))

sqrt(0.9571) # calculate r

# correlation figure
ggplot(dat22, aes(x = dry_cone_weight_g, y = predicted_dry_cone_weight_g)) + 
  geom_point() + 
  ylim(0, 225) + 
  xlim(0, 225) +
  geom_smooth(method = "lm", se = FALSE, color = "#CC79A7") + 
  annotate("text", x = 70, y = 225, parse = TRUE, label = "paste(italic(B) [1], \" = 1.03; \", italic(r), \" = 0.978; \", italic(p), \" ≤ 0.001 \")", size = 5) +
  xlab("Observed Dry Cone Yield (g)") + 
  ylab("Predicted Dry Cone Yield (g)") + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16, face = "bold"))

#### prepare for linear models ####
# make variables factors
spacing_trial_dat$rep <- as.factor(spacing_trial_dat$rep)
spacing_trial_dat$year <- as.factor(spacing_trial_dat$year)
spacing_trial_dat$spacing <- as.factor(spacing_trial_dat$spacing)

## calculate average cone weight in hop box
spacing_trial_dat <- spacing_trial_dat %>% mutate(hop_box_weight_g_per_cone = hop_box_weight_g / cone_ct)
hist(spacing_trial_dat$hop_box_weight_g_per_cone)
spacing_trial_dat <- spacing_trial_dat %>% mutate(hop_box_weight_g_per_cone_adj_dm = hop_box_weight_g_per_cone * dry_matter_percent)

#### calculate the means of the 3.5 foot spacings ####
spacing_trial_dat_mean <- spacing_trial_dat %>% 
  unite(spacing_trial_means, genotype, year, spacing, rep, sep = "_") %>% 
  group_by(spacing_trial_means) %>% 
  summarise(mean_percent_up_string = mean(percent_up_string, na.rm = T), 
            mean_arm_length = mean(arm_length, na.rm = T), 
            mean_hop_box_weight_g_per_cone = mean(hop_box_weight_g_per_cone_adj_dm, na.rm = T), 
            mean_harvest_index = mean(harvest_index, na.rm = T), 
            mean_area = mean(area, na.rm = T), 
            mean_perimeter = mean(perimeter, na.rm = T), 
            mean_length = mean(length, na.rm = T), 
            mean_width = mean(width, na.rm = T), 
            mean_openness = mean(openness, na.rm = T), 
            mean_UV_AA = mean(UV_AA, na.rm = T), 
            mean_UV_BA = mean(UV_BA, na.rm = T), 
            mean_UV_HSI = mean(UV_HSI, na.rm = T), 
            mean_predicted_dry_cone_weight_g = mean(predicted_dry_cone_weight_g, na.rm = T), 
            mean_cone_density = mean_hop_box_weight_g_per_cone / mean_area) %>%
  separate(spacing_trial_means, into = c("genotype", "year", "spacing", "rep"), sep = "_")

#### evaluate how the genotypes in the trial performed compared with their commercial historical values for yield #### 
space_dry_cone_weight_g <- lmer(log(mean_predicted_dry_cone_weight_g + 0.1) ~ genotype*year + (1|rep) + (1|year:rep), data = spacing_trial_dat_mean %>% filter(spacing == 3.5))
emm_space_dry_cone_weight_g <- as.data.frame(emmeans(space_dry_cone_weight_g, ~ genotype*year, type = "response")) %>% 
  dplyr::select(genotype, year, response)

# read in historical yield blups
yield_lmer_pred <- read.csv("./Data Analysis Results/historical_blups.csv")
# do some minor editing to the genotype names
yield_lmer_pred <- yield_lmer_pred %>% mutate(genotype = gsub("MountRainier", "Mount Rainier", grp), 
                           genotype = gsub("MountHood", "Mount Hood", genotype)) %>% 
  dplyr::select(genotype, est.mean)

# join two datasets
yield_comparison <- left_join(emm_space_dry_cone_weight_g, yield_lmer_pred, by = "genotype") %>%
  mutate(trial_yield = (response * 889)/454)

# filter out the years
yield22 <- yield_comparison %>% filter(year == 2022)
yield23 <- yield_comparison %>% filter(year == 2023)

# calculate spearman rank correlations
cor.test(x = yield22$response, y = yield22$est.mean, method = "spearman") # rho = -0.0857, p = 0.9194
cor.test(x = yield23$response, y = yield23$est.mean, method = "spearman") # rho = 0.942; p = 0.01667

#### determine the correlation between cone traits to determine if we need to actually analyze them separately ####
spacing_corr <- spacing_trial_dat_mean %>% dplyr::select(year, mean_percent_up_string, mean_arm_length, 
                                                         mean_predicted_dry_cone_weight_g, mean_harvest_index, 
                                                         mean_UV_AA, mean_UV_BA, mean_area, mean_length, mean_width, 
                                                         mean_perimeter, mean_openness, mean_hop_box_weight_g_per_cone, 
                                                         mean_cone_density)

# separate out by years for consistency
spacing_corr_22 <- spacing_corr %>% filter(year == "2022") %>% dplyr::select(-year)
spacing_corr_23 <- spacing_corr %>% filter(year == "2023") %>% dplyr::select(-year)

# first analyze 22, then edit year to 23 below
cor <- rcorr(as.matrix(spacing_corr_22), type = "pearson")
cor.p <- as.data.frame(cor$P)
cor.r <- as.data.frame(cor$r)
cor.r <- cor.r %>% mutate_all(funs(round(., 2))) # round r to two digits
cor.p.stars <- cor.p
for (k in 1:ncol(cor.p)) {
  for (j in 1:nrow(cor.p)) {
    cor.p.stars[k,j] <- stars.pval(cor.p[k,j])
  }
}
cor.r.p <- cor.r

for (k in 1:ncol(cor.p)) {
  for (j in 1:nrow(cor.r)) {
    cor.r.p[k,j] <- paste(cor.r[k,j], cor.p.stars[k,j], sep = "")
  }
}

cor.r.p[upper.tri(cor.r.p)] <- NA


# change year here too
write.csv(cor.r.p, "./Data Analysis Results/correlations_rp_22.csv", row.names = F)

#### evaluate the effect of position within the 3.5' spacing ####
# add a column to designate "hill"
spacing_trial_dat <- spacing_trial_dat %>% filter(spacing == 3.5) %>% arrange(year, row, position) %>% mutate(hill = as.factor(rep(1:7, 56)))

# work through traits (transformation)
# percent_up_string (none), arm_length (none), predicted_dry_cone_weight_g (none), harvest_index (none), UV_AA, UV_BA, area (log + 0.1), openness, hop_box_weight_g_per_cone

# percent_up_string – NS, NS 
# arm_length - NS, cubic 
# predicted_dry_cone_weight_g - NS, 0.002 (quadratic) 
# harvest_index – NS, NS 
# UV_AA – NS, NS 
# UV_BA – NS, NS 
# area – NS, 0.07 (quadratic) 
# openness – NA, NA 
# hop_box_weight_g_per_cone – 0.03 (quadratic), NS 

# 2022
lmer_yield_3.5 <- lmer(predicted_dry_cone_weight_g ~ genotype*hill + (1|rep), data = (spacing_trial_dat %>% filter(year == "2022")))
resid_panel(lmer_yield_3.5)
anova(lmer_yield_3.5) # no effect of hill, no interaction between hill and genotype
plot(lmer_yield_3.5)
yield_em_3.5 <- emmeans(lmer_yield_3.5, ~ hill)
contrast(yield_em_3.5, "poly")

# 2023
lmer_yield_3.5 <- lmer(predicted_dry_cone_weight_g ~ genotype*hill + (1|rep), data = (spacing_trial_dat %>% filter(year == "2023")))
resid_panel(lmer_yield_3.5)
anova(lmer_yield_3.5) # significant effect of hill, no interaction
yield_em_3.5 <- emmeans(lmer_yield_3.5, ~ hill)
contrast(yield_em_3.5, "poly") # looking at contrasts, significant quadratic relationship

# figure for 2022/2023 (replace year)
dat22 <- spacing_trial_dat %>% filter(year == "2023")
dat22$hill <- as.numeric(dat22$hill)
spacing_trial_dat_mean_year <- dat22 %>% group_by(genotype, hill) %>% summarise(mean_yield = mean(predicted_dry_cone_weight_g, na.rm = T))

yield_plot_23 <- ggplot(spacing_trial_dat_mean_year, aes(x=hill, y=mean_yield)) + 
  geom_point() + 
  stat_smooth(se = FALSE, color = "#009E73") + 
  theme_bw() + 
  labs(y = "Average Predicted Dry Cone Yield (g)",
       x = "Hill Position", 
       title = "2023") +
  scale_y_continuous(limits = c(0, 800)) +
  theme(plot.subtitle = element_text(size = 15), 
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold"))


cowplot::plot_grid(yield_plot_22, yield_plot_23)
# export as 1500 x 900  

#### mean_predicted_dry_cone_weight_g ####
# first create a dataframe to convert feet to meters
meter_conv <- data.frame(spacing = c(1.5, 2, 3.5, 5, 7), spacing_m = c(0.46, 0.61, 1.07, 1.52, 2.13))
meter_conv <- meter_conv %>% mutate(spacing = as.factor(spacing), 
                      spacing_m = as.factor(spacing_m))

# 2022
lmer_mean_predicted_dry_cone_weight_g_22 <- lmer(log(mean_predicted_dry_cone_weight_g + 0.01) ~ genotype * spacing + (1|rep) +  (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2022"))
#resid_panel(lmer_mean_predicted_dry_cone_weight_g_22)
anova(lmer_mean_predicted_dry_cone_weight_g_22)
write.csv(anova(lmer_mean_predicted_dry_cone_weight_g_22), "./Data Analysis Results/Anovas/anova_mean_predicted_dry_cone_weight_g_22.csv")

# effect of spacing
em_spacing <- emmeans(lmer_mean_predicted_dry_cone_weight_g_22, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_predicted_dry_cone_weight_g_spacing_22.csv", row.names = F)

# effect of year
lmer_mean_predicted_dry_cone_weight_g <- lmer(log(mean_predicted_dry_cone_weight_g + 0.01) ~ year*spacing + (1|year:rep), data = spacing_trial_dat_mean %>% filter(spacing %in% c(3.5, 5, 7)))
anova(lmer_mean_predicted_dry_cone_weight_g)
em_year <- emmeans(lmer_mean_predicted_dry_cone_weight_g, ~ year | spacing, type = "response")
cld <- data.frame(cld(em_year, Letters = c(letters)))
write.csv(cld, "./Data Analysis Results/mean_predicted_dry_cone_weight_g_year.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_predicted_dry_cone_weight_g_22, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(response) %>% dplyr::select(genotype)
neworder <- c("Sterling", "Mount Rainier", "Vanguard", "Mount Hood", "Magnum", "Bitter Gold", "Comet")
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, response) %>% pivot_wider(names_from = spacing, values_from = response)
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.25; P = NS (0.59)
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.50; P = NS (0.27)
#plot(dat_cor$`3.5`, dat_cor$`2`) # p = 0.2357; NS

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.96; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`5`) # p = 0.002778 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 1; P < 0.001
#plot(dat_cor$`3.5`, dat_cor$`7`) # p = 0.0003968

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue = c("ρ = 0.25; P = NS", "ρ = 0.25; P = NS", " ", "ρ = 0.96; P < 0.01", "ρ = 1; P < 0.001"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

# plot yield bar plot for 2022
yield22 <- ggplot(em_all2, aes(genotype, response, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "Dry Cone Yield (g)",
       x = "Genotype", fill = "Genotype", 
       title = "2022") +
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 165, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 175)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "bottom",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) + 
  scale_fill_manual(values = c("#009E73", "#CC79A7", "#E69F00", "#D55E00", "#F0E442", "#56B4E9", "#0072B2")) #1200 x 900


neworder <- c("Sterling", "Mount Rainier", "Vanguard", "Mount Hood", "Magnum", "Bitter Gold", "Comet")

# 2023
lmer_mean_predicted_dry_cone_weight_g_23 <- lmer(sqrt(mean_predicted_dry_cone_weight_g) ~ genotype * spacing + (1|rep) + (1|spacing:rep), data = spacing_trial_dat_mean %>% filter(year == "2023"))
#resid_panel(lmer_mean_predicted_dry_cone_weight_g_23)
anova(lmer_mean_predicted_dry_cone_weight_g_23)
write.csv(anova(lmer_mean_predicted_dry_cone_weight_g_23), "./Data Analysis Results/Anovas/anova_mean_predicted_dry_cone_weight_g_23.csv")

em_spacing <- emmeans(lmer_mean_predicted_dry_cone_weight_g_23, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_predicted_dry_cone_weight_g_spacing_23.csv", row.names = F)


# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_predicted_dry_cone_weight_g_23, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(response) %>% dplyr::select(genotype)
# set new order, and arrange data so the bars show up in the correct order in the figure
neworder <- c("Vanguard", "Mount Hood", "Mount Rainier", "Magnum", "Comet", "Sterling", "Bitter Gold") # edit new order
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, response) %>% pivot_wider(names_from = spacing, values_from = response) 
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.18; P = NS
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.39; P = NS
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.11; P = NS
#plot(dat_cor$`3.5`, dat_cor$`5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.64; P = NS
#plot(dat_cor$`3.5`, dat_cor$`7`) 

# create a dataframe of the correlation values for the figure
pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue = c("ρ = 0.18; P = NS", "ρ = 0.39; P = NS", " ", "ρ = 0.11; P = NS", "ρ = 0.64; P = NS"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

# plot for 2023
yield23 <- ggplot(em_all2, aes(genotype, response, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "",
       x = "Genotype", fill = "Genotype",
       title = "2023") +
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 1400, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 1500)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "none",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#E69F00", "#D55E00", "#CC79A7", "#F0E442", "#0072B2", "#009E73", "#56B4E9")) #1200 x 900

# reminder of the colors!
# Sterling: #009E73
# Mount Rainier: #CC79A7
# Vanguard: #E69F00
# Mount Hood: #D55E00
# Magnum: #F0E442
# Bitter Gold: #56B4E9
# Comet: #0072B2

legend <- get_legend(yield22)
yield22 <- yield22 + theme(legend.position = "none")
top_row <- plot_grid(yield22, yield23)
plot_grid(top_row, legend, nrow = 2, rel_heights = c(1, 0.1))
# 2000 x 800 TIFF



#### UV_AA #####
# 2022
lmer_mean_UV_AA_22 <- lmer(mean_UV_AA ~ genotype * spacing + (1|rep) +  (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2022")) #%>%
  #filter(! genotype %in% c("Magnum", "Sterling"))) ### ALERT: filtered out mangum and sterling because they were missing data in seedlings!!
resid_panel(lmer_mean_UV_AA_22)
anova(lmer_mean_UV_AA_22)
write.csv(anova(lmer_mean_UV_AA_22), "./Data Analysis Results/Anovas/anova_mean_UV_AA_22.csv")


# effect of spacing
em_spacing <- emmeans(lmer_mean_UV_AA_22, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_UV_AA_spacing_22.csv", row.names = F)

# effect of year
lmer_mean_UV_AA <- lmer(mean_UV_AA ~ year*spacing + (1|year:rep), data = spacing_trial_dat_mean %>% filter(spacing %in% c(3.5, 5, 7)))
anova(lmer_mean_UV_AA)
em_year <- emmeans(lmer_mean_UV_AA, ~ year | spacing, type = "response")
cld <- data.frame(cld(em_year, Letters = c(letters)))
write.csv(cld, "./mean_UV_AA_year.csv", row.names = F)


# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_UV_AA_22, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Mount Rainier", "Mount Hood", "Vanguard", "Sterling", "Bitter Gold", "Magnum", "Comet")
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)), genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean)
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") #  ρ = 0.43; P = NS
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.30; P = NS
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.68; P = NS
#plot(dat_cor$`3.5`, dat_cor$`5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.96; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`7`) 

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue = c("ρ = 0.43; P = NS", "ρ = 0.30; P = NS", " ", "ρ = 0.68; P = NS", "ρ = 0.96; P < 0.01"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")


AA_22 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "Alpha Acid (%)",
       x = "Genotype", fill = "Genotype", 
       title = "2022") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 17, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 18)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "bottom",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) + 
  scale_fill_manual(values = c("#CC79A7", "#D55E00", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#0072B2")) #1200 x 900


neworder <- c("Mount Rainier", "Mount Hood", "Vanguard", "Sterling", "Bitter Gold", "Magnum", "Comet")

# 2023
lmer_mean_UV_AA_23 <- lmer(mean_UV_AA ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2023"))
resid_panel(lmer_mean_UV_AA_23)
anova(lmer_mean_UV_AA_23)
write.csv(anova(lmer_mean_UV_AA_23), "./Data Analysis Results/Anovas/anova_mean_UV_AA_23.csv")

em_spacing <- emmeans(lmer_mean_UV_AA_23, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_UV_AA_spacing_23.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_UV_AA_23, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Mount Hood", "Vanguard", "Sterling", "Mount Rainier", "Comet", "Bitter Gold", "Magnum") # edit new order
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean) 
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.89; P = 0.01
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.86; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`2`) # p = 0.0123

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.93; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`5`) # p = 0.006746 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 1; P < 0.001
#plot(dat_cor$`3.5`, dat_cor$`7`) # p = 0.0003968

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue = c("ρ = 0.89; P = 0.01", "ρ = 0.86; P < 0.05", " " , "ρ = 0.93; P < 0.01", "ρ = 1; P < 0.001"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")


AA_23 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "",
       x = "Genotype", fill = "Genotype",
       title = "2023") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 17, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 18)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "none",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#D55E00", "#E69F00", "#009E73", "#CC79A7", "#0072B2", "#56B4E9", "#F0E442")) #1200 x 900

neworder <- c("Mount Hood", "Vanguard", "Sterling", "Mount Rainier", "Comet", "Bitter Gold", "Magnum") # edit new order

# legend <- get_legend(yield22) # use the same legend as before
AA_22 <- AA_22 + theme(legend.position = "none")
top_row <- plot_grid(AA_22, AA_23)
plot_grid(top_row, legend, nrow = 2, rel_heights = c(1, 0.1))
# 2000 x 800 TIFF


#### percent up string ####
# 2022
lmer_mean_percent_up_string_22 <- lmer(mean_percent_up_string ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2022"))
resid_panel(lmer_mean_percent_up_string_22)
anova(lmer_mean_percent_up_string_22)
write.csv(anova(lmer_mean_percent_up_string_22), "./Data Analysis Results/Anovas/anova_mean_percent_up_string_22.csv")

# effect of spacing
em_spacing <- emmeans(lmer_mean_percent_up_string_22, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_percent_up_string_spacing_22.csv", row.names = F)

# effect of year
lmer_mean_percent_up_string <- lmer(mean_percent_up_string ~ year*spacing + (1|year:rep), data = spacing_trial_dat_mean %>% filter(spacing %in% c(3.5, 5, 7)))
anova(lmer_mean_percent_up_string)
em_year <- emmeans(lmer_mean_percent_up_string, ~ year | spacing, type = "response")
cld <- data.frame(cld(em_year, Letters = c(letters)))
write.csv(cld, "./Data Analysis Results/mean_percent_up_string_year.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_percent_up_string_22, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Mount Rainier", "Sterling", "Mount Hood", "Magnum", "Vanguard", "Comet", "Bitter Gold")
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean)
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.64; P = NS
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = -.04; P = NS
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.86; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 1; P < 0.001
#plot(dat_cor$`3.5`, dat_cor$`7`) 

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue = c("ρ = 0.64; P = NS", "ρ = -.04; P = NS", " ", "ρ = 0.86; P < 0.05", "ρ = 1; P < 0.001"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

percent_up_string_22 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "Vigor (%)",
       x = "Genotype", fill = "Genotype", 
       title = "2022") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 110, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 120)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "bottom",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) + 
  scale_fill_manual(values = c("#CC79A7", "#009E73", "#D55E00", "#F0E442", "#E69F00", "#0072B2", "#56B4E9")) #1200 x 900


neworder <- c("Mount Rainier", "Sterling", "Mount Hood", "Magnum", "Vanguard", "Comet", "Bitter Gold")
# Sterling: #009E73
# Mount Rainier: #CC79A7
# Vanguard: #E69F00
# Mount Hood: #D55E00
# Magnum: #F0E442
# Bitter Gold: #56B4E9
# Comet: #0072B2

# 2023
lmer_mean_percent_up_string_23 <- lmer(mean_percent_up_string ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2023"))
hist(spacing_trial_dat_mean$mean_percent_up_string)
resid_panel(lmer_mean_percent_up_string_23)
anova(lmer_mean_percent_up_string_23)
write.csv(anova(lmer_mean_percent_up_string_23), "./Data Analysis Results/Anovas/anova_mean_percent_up_string_23.csv")

em_spacing <- emmeans(lmer_mean_percent_up_string_23, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_percent_up_string_spacing_23.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_percent_up_string_23, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Mount Rainier", "Magnum", "Sterling", "Mount Hood", "Vanguard", "Bitter Gold", "Comet") # edit new order
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean) 
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") #NA
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman")
#plot(dat_cor$`3.5`, dat_cor$`2`) # NA

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.83; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`5`) # p = 0.04802 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.29; P = NS
#plot(dat_cor$`3.5`, dat_cor$`7`) # p = 0.556

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m)[3:5], pvalue =  c( " " ,"ρ = 0.83; P < 0.05", "ρ = 0.29; P = NS"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

percent_up_string_23 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "",
       x = "Genotype", fill = "Genotype",
       title = "2023") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 110, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 120)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "none",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442", "#009E73", "#D55E00", "#E69F00", "#56B4E9", "#0072B2")) #1200 x 900

neworder <- c("Mount Rainier", "Magnum", "Sterling", "Mount Hood", "Vanguard", "Bitter Gold", "Comet") # edit new order

# legend <- get_legend(yield22) # use the same legend as before
percent_up_string_22 <- percent_up_string_22 + theme(legend.position = "none")
top_row <- plot_grid(percent_up_string_22, percent_up_string_23, rel_widths = c(1, .6))
plot_grid(top_row, legend, nrow = 2, rel_heights = c(1, 0.1))
# 2000 x 800 TIFF


#### arm length ####

# 2022
lmer_mean_arm_length_22 <- lmer(sqrt(mean_arm_length) ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2022"))
resid_panel(lmer_mean_arm_length_22)
anova(lmer_mean_arm_length_22)
write.csv(anova(lmer_mean_arm_length_22), "./Data Analysis Results/Anovas/anova_mean_arm_length_22.csv")

# effect of spacing
em_spacing <- emmeans(lmer_mean_arm_length_22, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_arm_length_spacing_22.csv", row.names = F)

# effect of year
lmer_mean_arm_length <- lmer(sqrt(mean_arm_length) ~ year*spacing + (1|year:rep), data = spacing_trial_dat_mean %>% filter(spacing %in% c(3.5, 5, 7)))
anova(lmer_mean_arm_length)
em_year <- emmeans(lmer_mean_arm_length, ~ year | spacing, type = "response")
cld <- data.frame(cld(em_year, Letters = c(letters)))
write.csv(cld, "./Data Analysis Results/mean_arm_length_spacing.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_arm_length_22, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(response) %>% dplyr::select(genotype)
neworder <- c("Mount Rainier", "Magnum", "Mount Hood", "Sterling", "Vanguard", "Bitter Gold", "Comet")
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, response) %>% pivot_wider(names_from = spacing, values_from = response)

cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman", exact=FALSE) # ρ = 0.77; P < 0.05

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman", exact = FALSE) # ρ = 0.68; P < 0.10

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman", exact = FALSE) # ρ = 0.75; P < 0.10

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman", exact = FALSE) # ρ = 0.71; P < 0.10

font_size <- 20

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue = c("ρ = 0.77; P < 0.05", "ρ = 0.68; P < 0.10", " ", "ρ = 0.75; P < 0.10", "ρ = 0.71; P < 0.10"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

arm_length_22 <- ggplot(em_all2, aes(genotype, response, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) +
  theme_bw() + 
  labs(y = "Lateral Length (Scale)",
       x = "Genotype", fill = "Genotype", 
       title = "2022") +
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9)) +
  geom_text(pval, mapping = aes(x = 4, y = 3.25, label = pvalue), size = 4, inherit.aes = FALSE) +
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 3.5)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "bottom",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#0072B2")) #1200 x 900


neworder <- c("Mount Rainier", "Magnum", "Mount Hood", "Sterling", "Vanguard", "Bitter Gold", "Comet")

# 2023
lmer_mean_arm_length_23 <- lmer(sqrt(mean_arm_length) ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2023"))
resid_panel(lmer_mean_arm_length_23)
anova(lmer_mean_arm_length_23)
write.csv(anova(lmer_mean_arm_length_23), "./Data Analysis Results/Anovas/anova_mean_arm_length_23.csv")

em_spacing <- emmeans(lmer_mean_arm_length_23, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_arm_length_spacing_23.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_arm_length_23, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(response) %>% dplyr::select(genotype)
neworder <- c("Magnum", "Mount Hood", "Sterling", "Vanguard", "Mount Rainier", "Bitter Gold", "Comet") # edit new order
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, response) %>% pivot_wider(names_from = spacing, values_from = response) 
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.72; P < 0.10
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.86; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`2`) # NA

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman", exact = FALSE) # ρ = 0.50; P = NS
#plot(dat_cor$`3.5`, dat_cor$`5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman", exact = FALSE) # ρ = 0.32; P = NS
#plot(dat_cor$`3.5`, dat_cor$`7`) 

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue = c( "ρ = 0.72; P < 0.10", "ρ = 0.86; P < 0.05", "", "ρ = 0.50; P = NS", "ρ = 0.32; P = NS"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

arm_length_23 <- ggplot(em_all2, aes(genotype, response, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "",
       x = "Genotype", fill = "Genotype",
       title = "2023") +
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 3.25, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 3.5)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "none",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#F0E442", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#0072B2")) #1200 x 900

neworder <- c("Magnum", "Mount Hood", "Sterling", "Vanguard", "Mount Rainier", "Bitter Gold", "Comet") # edit new order

legend <- get_legend(arm_length_22) # use the same legend as before
arm_length_22 <- arm_length_22 + theme(legend.position = "none")
top_row <- cowplot::plot_grid(arm_length_22, arm_length_23)
cowplot::plot_grid(top_row, legend, nrow = 2, rel_heights = c(1, 0.1))

# 2000 x 800 TIFF


#### UV_BA ####

# 2022
lmer_mean_mean_UV_BA_22 <- lmer(mean_UV_BA ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2022"))
resid_panel(lmer_mean_mean_UV_BA_22)
anova(lmer_mean_mean_UV_BA_22)
write.csv(anova(lmer_mean_mean_UV_BA_22), "./Data Analysis Results/Anovas/anova_mean_UV_BA_22.csv")

# effect of spacing
em_spacing <- emmeans(lmer_mean_mean_UV_BA_22, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_UV_BA_spacing_22.csv", row.names = F)

# effect of year
lmer_mean_UV_BA <- lmer(mean_UV_BA ~ year*spacing + (1|year:rep), data = spacing_trial_dat_mean %>% filter(spacing %in% c(3.5, 5, 7)))
anova(lmer_mean_UV_BA)
em_year <- emmeans(lmer_mean_UV_BA, ~ year | spacing, type = "response")
cld <- data.frame(cld(em_year, Letters = c(letters)))
write.csv(cld, "./Data Analysis Results/mean_UV_BA_year.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_mean_UV_BA_22, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Sterling", "Mount Rainier", "Bitter Gold", "Vanguard", "Mount Hood", "Magnum", "Comet")
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean)
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 1; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.90; P < 0.10
#plot(dat_cor$`3.5`, dat_cor$`2`) # p = 0.3536; NS

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.93; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`5`) # p = 0.006746

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.85; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`7`) # 0.1389

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue =  c("ρ = 1; P < 0.01", "0.90; P < 0.10", "", "ρ = 0.93; P < 0.01", "ρ = 0.85; P < 0.05"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

UV_BA_22 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "UV_BA",
       x = "Genotype", fill = "Genotype", 
       title = "2022") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 4.75, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 5)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "bottom",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) + 
  scale_fill_manual(values = c("#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#D55E00", "#F0E442", "#0072B2")) #1200 x 900


neworder <- c("Sterling", "Mount Rainier", "Bitter Gold", "Vanguard", "Mount Hood", "Magnum", "Comet")
# Sterling: #009E73
# Mount Rainier: #CC79A7
# Vanguard: #E69F00
# Mount Hood: #D55E00
# Magnum: #F0E442
# Bitter Gold: #56B4E9
# Comet: #0072B2

# 2023
lmer_mean_UV_BA_23 <- lmer(mean_UV_BA ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2023"))
resid_panel(lmer_mean_UV_BA_23)
anova(lmer_mean_UV_BA_23)
write.csv(anova(lmer_mean_UV_BA_23), "./Data Analysis Results/Anovas/anova_mean_UV_BA_23.csv")

em_spacing <- emmeans(lmer_mean_UV_BA_23, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_UV_BA_spacing_23.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_UV_BA_23, ~ genotype * spacing, type = "response"))
ba_23 <- em_all
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Comet", "Sterling", "Mount Hood", "Vanguard", "Magnum", "Mount Rainier", "Bitter Gold") # edit new order
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean) 
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.07; P = NS
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.25; P = NS
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.96; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`5`)

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 1; P < 0.001
#plot(dat_cor$`3.5`, dat_cor$`7`) # p = 0.08967

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue =   c( "ρ = 0.07; P = NS", "ρ = 0.25; P = NS", "", "ρ = 0.96; P < 0.01", "ρ = 1; P < 0.001"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

UV_BA_23 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "",
       x = "Genotype", fill = "Genotype",
       title = "2023") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 7.5, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 8)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "none",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#F0E442", "#CC79A7", "#56B4E9")) #1200 x 900

neworder <- c("Comet", "Sterling", "Mount Hood", "Vanguard", "Magnum", "Mount Rainier", "Bitter Gold") # edit new order

# legend <- get_legend(yield22) # use the same legend as before
UV_BA_22 <- UV_BA_22 + theme(legend.position = "none")
top_row <- cowplot::plot_grid(UV_BA_22, UV_BA_23, rel_widths = c(1, 1))
cowplot::plot_grid(top_row, legend, nrow = 2, rel_heights = c(1, 0.1))
# 2000 x 800 TIFF


#### harvest_index ####

# 2022
lmer_mean_harvest_index_22 <- lmer(mean_harvest_index ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2022"))
resid_panel(lmer_mean_harvest_index_22)
anova(lmer_mean_harvest_index_22)
write.csv(anova(lmer_mean_harvest_index_22), "./Data Analysis Results/Anovas/anova_mean_harvest_index_22.csv")

# effect of spacing
em_spacing <- emmeans(lmer_mean_harvest_index_22, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_harvest_index_spacing_22.csv", row.names = F)

# effect of year
lmer_mean_harvest_index <- lmer(mean_harvest_index ~ year*spacing + (1|year:rep), data = spacing_trial_dat_mean %>% filter(spacing %in% c(3.5, 5, 7)))
anova(lmer_mean_harvest_index)
em_year <- emmeans(lmer_mean_harvest_index, ~ year | spacing, type = "response")
cld <- data.frame(cld(em_year, Letters = c(letters)))
write.csv(cld, "./Data Analysis Results/mean_harvest_index_year.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_harvest_index_22, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Mount Rainier", "Sterling", "Vanguard", "Mount Hood", "Magnum", "Bitter Gold", "Comet")
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean)
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.11; P = NS
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.43; P = NS
#plot(dat_cor$`3.5`, dat_cor$`2`) # 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.92; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`5`) # 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.64; P = 0.10
#plot(dat_cor$`3.5`, dat_cor$`7`) # 

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue =  c("ρ = 0.11; P = NS", "ρ = 0.43; P = NS", " ", "ρ = 0.92; P < 0.01", "ρ = 0.64; P = 0.10"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

harvest_index_22 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "Harvest Index (%)",
       x = "Genotype", fill = "Genotype", 
       title = "2022") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = .45, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, .5)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "bottom",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) + 
  scale_fill_manual(values = c("#CC79A7", "#009E73", "#E69F00", "#D55E00", "#F0E442", "#56B4E9", "#0072B2")) #1200 x 900


neworder <- c("Mount Rainier", "Sterling", "Vanguard", "Mount Hood", "Magnum", "Bitter Gold", "Comet")


# 2023
lmer_mean_harvest_index_23 <- lmer(mean_harvest_index ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2023"))
resid_panel(lmer_mean_harvest_index_23)
anova(lmer_mean_harvest_index_23)
write.csv(anova(lmer_mean_harvest_index_23), "./Data Analysis Results/Anovas/anova_mean_harvest_index_23.csv")

em_spacing <- emmeans(lmer_mean_harvest_index_23, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_harvest_index_spacing_23.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_harvest_index_23, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Vanguard", "Comet", "Mount Hood", "Magnum", "Bitter Gold", "Sterling", "Mount Rainier") # edit new order
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean) 
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.86; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.21; P = NS
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.57; P = NS
#plot(dat_cor$`3.5`, dat_cor$`5`)

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.93; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`7`) # p = 0.08967

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue =  c("ρ = 0.86; P < 0.05", "ρ = 0.21; P = NS", "", "ρ = 0.57; P = NS", "ρ = 0.93; P < 0.01"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

harvest_index_23 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "",
       x = "Genotype", fill = "Genotype",
       title = "2023") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = .45, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, .5)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "none",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2", "#D55E00", "#F0E442", "#56B4E9", "#009E73", "#CC79A7")) #1200 x 900

neworder <- c("Vanguard", "Comet", "Mount Hood", "Magnum", "Bitter Gold", "Sterling", "Mount Rainier") # edit new order

# legend <- get_legend(yield22) # use the same legend as before
harvest_index_22 <- harvest_index_22 + theme(legend.position = "none")
top_row <- cowplot::plot_grid(harvest_index_22, harvest_index_23, rel_widths = c(1, 1))
cowplot::plot_grid(top_row, legend, nrow = 2, rel_heights = c(1, 0.1))
# 2000 x 800 TIFF

#### area ####
# 2022
colnames(spacing_trial_dat_mean)
lmer_mean_area_22 <- lmer(mean_area ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2022"))
resid_panel(lmer_mean_area_22)
anova(lmer_mean_area_22)
write.csv(anova(lmer_mean_area_22), "./Data Analysis Results/Anovas/anova_mean_area_22.csv")

# effect of spacing
em_spacing <- emmeans(lmer_mean_area_22, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_cone_area_spacing_22.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_area_22, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Vanguard", "Sterling", "Mount Hood", "Bitter Gold", "Mount Rainier", "Magnum", "Comet")
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean)
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.94; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.93; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 1; P < 0.001
#plot(dat_cor$`3.5`, dat_cor$`5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.96; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`7`) 

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue =  c("ρ = 0.94; P < 0.05", "ρ = 0.93; P < 0.01", "", "ρ = 1; P < 0.001", "ρ = 0.96; P < 0.01"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

cone_area_22 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "Cone Area (cm2)",
       x = "Genotype", fill = "Genotype", 
       title = "2022") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 11.5, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 12)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "bottom",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) + 
  scale_fill_manual(values = c("#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7", "#F0E442", "#0072B2")) #1200 x 900


neworder <- c("Vanguard", "Sterling", "Mount Hood", "Bitter Gold", "Mount Rainier", "Magnum", "Comet")

# 2023
lmer_mean_cone_area_23 <- lmer(mean_area ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2023"))
resid_panel(lmer_mean_cone_area_23)
anova(lmer_mean_cone_area_23)
write.csv(anova(lmer_mean_cone_area_23), "./Data Analysis Results/Anovas/anova_mean_area_23.csv")

# effect of spacing
em_spacing <- emmeans(lmer_mean_cone_area_23, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_cone_area_spacing_23.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_cone_area_23, ~ genotype * spacing, type = "response"))
cone_area_23 <- em_all
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Sterling", "Vanguard", "Mount Hood", "Mount Rainier", "Bitter Gold", "Magnum", "Comet") # edit new order
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean) 
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.93; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.71; P < 0.10
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.75; P < 0.10
#plot(dat_cor$`3.5`, dat_cor$`5`)

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.71; P < 0.10
#plot(dat_cor$`3.5`, dat_cor$`7`) # p = 0.08967

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue =  c("ρ = 0.93; P < 0.01", "ρ = 0.71; P < 0.10", "", "ρ = 0.75; P < 0.10", "ρ = 0.71; P < 0.10"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

cone_area_23 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "",
       x = "Genotype", fill = "Genotype",
       title = "2023") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 11.5, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 12)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "none",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#009E73", "#E69F00", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#0072B2")) #1200 x 900

neworder <- c("Sterling", "Vanguard", "Mount Hood", "Mount Rainier", "Bitter Gold", "Magnum", "Comet") # edit new order

# legend <- get_legend(yield22) # use the same legend as before
cone_area_22 <- cone_area_22 + theme(legend.position = "none")
top_row <- cowplot::plot_grid(cone_area_22, cone_area_23, rel_widths = c(1, 1))
cowplot::plot_grid(top_row, legend, nrow = 2, rel_heights = c(1, 0.1))
# 2000 x 800 TIFF


#### cone density  ####

# 2022
lmer_mean_cone_density_22 <- lmer(log10(mean_cone_density) ~ genotype * spacing + (1|rep) + (1|rep:spacing),  data = spacing_trial_dat_mean %>% filter(year == "2022"))
resid_panel(lmer_mean_cone_density_22)
anova(lmer_mean_cone_density_22)
write.csv(anova(lmer_mean_cone_density_22), "./Data Analysis Results/Anovas/anova_mean_cone_density_22.csv")

# effect of spacing
em_spacing <- emmeans(lmer_mean_cone_density_22, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_cone_density_spacing_22.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_cone_density_22, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(response) %>% dplyr::select(genotype)
neworder <- c("Mount Rainier", "Sterling", "Mount Hood", "Vanguard", "Comet", "Bitter Gold", "Magnum")
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, response) %>% pivot_wider(names_from = spacing, values_from = response)
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.77; P = 0.10
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.61; P = NS
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.93; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.79; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`7`) 

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue =  c("ρ = 0.77; P = 0.10", "ρ = 0.61; P = NS", "", "ρ = 0.93; P < 0.01", "ρ = 0.79; P < 0.05"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")
em_all2 <- em_all2 %>% mutate(response = response * 1000)

cone_density_22 <- ggplot(em_all2, aes(genotype, response, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = bquote('Adjusted Cone Density ' (mg / cm^2)),
       x = "Genotype", fill = "Genotype", 
       title = "2022") + 
  geom_errorbar(aes(ymin=(response-(SE * 1000)), ymax=(response+(SE * 1000))), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = (55), label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 60)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "bottom",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) + 
  scale_fill_manual(values = c("#CC79A7", "#009E73", "#D55E00", "#E69F00", "#0072B2", "#56B4E9", "#F0E442")) #1200 x 900


neworder <- c("Mount Rainier", "Sterling", "Mount Hood", "Vanguard", "Comet", "Bitter Gold", "Magnum")

# 2023
lmer_mean_cone_density_23 <- lmer(mean_cone_density ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2023"))
resid_panel(lmer_mean_cone_density_23)
anova(lmer_mean_cone_density_23)
write.csv(anova(lmer_mean_cone_density_23), "./Data Analysis Results/Anovas/anova_mean_cone_density_23.csv")

# effect of spacing
em_spacing <- emmeans(lmer_mean_cone_density_23, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/mean_cone_density_spacing_23.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_cone_density_23, ~ genotype * spacing, type = "response"))
cone_density_23 <- em_all
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Mount Rainier", "Vanguard", "Comet", "Sterling", "Mount Hood", "Bitter Gold", "Magnum") # edit new order
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean) 
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.75; P < 0.10
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.57; P = NS
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.92; P < 0.01
#plot(dat_cor$`3.5`, dat_cor$`5`)

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 1; P < 0.001
#plot(dat_cor$`3.5`, dat_cor$`7`) # p = 0.08967

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue =  c("ρ = 0.75; P < 0.10", "ρ = 0.57; P = NS", "", "ρ = 0.92; P < 0.01", "ρ = 1; P < 0.001"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

em_all2 <- em_all2 %>% mutate(emmean = emmean * 1000)

cone_density_23 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "",
       x = "Genotype", fill = "Genotype",
       title = "2023") +
  geom_errorbar(aes(ymin=(emmean-(SE * 1000)), ymax=(emmean+(SE * 1000))), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 55, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 60)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "none",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#0072B2", "#009E73", "#D55E00", "#56B4E9", "#F0E442")) #1200 x 900

neworder <- c("Mount Rainier", "Vanguard", "Comet", "Sterling", "Mount Hood", "Bitter Gold", "Magnum") # edit new order

# legend <- get_legend(yield22) # use the same legend as before
cone_density_22 <- cone_density_22 + theme(legend.position = "none")
top_row <- cowplot::plot_grid(cone_density_22, cone_density_23, rel_widths = c(1, 1))
cowplot::plot_grid(top_row, legend, nrow = 2, rel_heights = c(1, 0.1))
# 2000 x 800 TIFF


#### openness ####

# 2022
colnames(spacing_trial_dat_mean)
lmer_openness_area_22 <- lmer(mean_openness ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2022"))
resid_panel(lmer_openness_area_22)
anova(lmer_openness_area_22)
write.csv(anova(lmer_openness_area_22), "./Data Analysis Results/Anovas/anova_mean_cone_openness_22.csv")

# effect of year
lmer_mean_openness <- lmer(mean_openness ~ year*spacing + (1|year:rep), data = spacing_trial_dat_mean %>% filter(spacing %in% c(3.5, 5, 7)))
anova(lmer_mean_openness)
em_year <- emmeans(lmer_mean_openness, ~ year | spacing, type = "response")
cld <- data.frame(cld(em_year, Letters = c(letters)))
write.csv(cld, "./Data Analysis Results/mean_openness_year.csv", row.names = F)

# effect of spacing
em_spacing <- emmeans(lmer_openness_area_22, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/openness_spacing_22.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_openness_area_22, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Sterling", "Bitter Gold", "Magnum", "Vanguard", "Mount Hood", "Mount Rainier", "Comet")
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean)
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.83; P < 0.10
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.86; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.89; P = 0.01
#plot(dat_cor$`3.5`, dat_cor$`5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.82; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`7`) 

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue =  c("ρ = 0.83; P < 0.10", "ρ = 0.86; P < 0.05", "", "ρ = 0.89; P = 0.01", "ρ = 0.82; P < 0.05"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

cone_openness_22 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "Cone Openness",
       x = "Genotype", fill = "Genotype", 
       title = "2022") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 5, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 5.5)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "bottom",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) + 
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#F0E442", "#E69F00", "#D55E00", "#CC79A7", "#0072B2")) #1200 x 900


neworder <- c("Sterling", "Bitter Gold", "Magnum", "Vanguard", "Mount Hood", "Mount Rainier", "Comet")

# 2023
lmer_mean_cone_openness_23 <- lmer(mean_openness ~ genotype * spacing + (1|rep) + (1|rep:spacing), data = spacing_trial_dat_mean %>% filter(year == "2023"))
resid_panel(lmer_mean_cone_openness_23)
anova(lmer_mean_cone_openness_23)
write.csv(anova(lmer_mean_cone_openness_23), "./Data Analysis Results/Anovas/anova_mean_cone_openness_23.csv")

# effect of spacing
em_spacing <- emmeans(lmer_mean_cone_openness_23, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 
write.csv(cld, "./Data Analysis Results/openness_spacing_23.csv", row.names = F)

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_cone_openness_23, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
neworder <- c("Sterling", "Bitter Gold", "Comet", "Magnum", "Mount Hood", "Vanguard", "Mount Rainier") # edit new order
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

# determine spearman rank correlations with the 3.5 "gold standard"
dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean) 
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # ρ = 0.71; P < 0.10
#plot(dat_cor$`3.5`, dat_cor$`1.5`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # ρ = 0.71; P < 0.10
#plot(dat_cor$`3.5`, dat_cor$`2`) 

cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # ρ = 0.86; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`5`)

cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # ρ = 0.86; P < 0.05
#plot(dat_cor$`3.5`, dat_cor$`7`) # p = 0.08967

pval <- data.frame(spacing_m = unique(meter_conv$spacing_m), pvalue =  c( "ρ = 0.71; P < 0.10", "ρ = 0.71; P < 0.10", "", "ρ = 0.86; P < 0.05", "ρ = 0.86; P < 0.05"))
em_all2 <- em_all2 %>% left_join(meter_conv, by = "spacing")

cone_openness_23 <- ggplot(em_all2, aes(genotype, emmean, fill = genotype)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~ spacing_m) + 
  theme_bw() + 
  labs(y = "",
       x = "Genotype", fill = "Genotype",
       title = "2023") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                position=position_dodge(.9)) + 
  geom_text(pval, mapping = aes(x = 4, y = 5, label = pvalue), size = 4, inherit.aes = FALSE) + 
  scale_y_continuous(expand = c(.01, 0), limits = c(0, 5.5)) +
  theme(legend.text = element_text(size = 20), 
        legend.position = "none",
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size), 
        plot.title = element_text(size = font_size, hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#0072B2", "#F0E442", "#D55E00", "#E69F00", "#CC79A7")) #1200 x 900

neworder <- c("Sterling", "Bitter Gold", "Comet", "Magnum", "Mount Hood", "Vanguard", "Mount Rainier") # edit new order

# legend <- get_legend(yield22) # use the same legend as before
cone_openness_22 <- cone_openness_22 + theme(legend.position = "none")
top_row <- cowplot::plot_grid(cone_openness_22, cone_openness_23, rel_widths = c(1, 1))
cowplot::plot_grid(top_row, legend, nrow = 2, rel_heights = c(1, 0.1))
# 2000 x 800 TIFF


#### compile all anovas #####
anovas <- list.files("./Data Analysis Results/Anovas/", full.names = T)
anova_list <- list() 

for (i in 1:length(anovas)) {
  anova_list[[i]] <- as.data.frame(read.csv(anovas[i]))
}

for (i in 1:length(anova_list)) {
  write.table(str_split(anovas[i], "/")[[1]][11], "./Data Analysis Results/Anovas_Concat.csv", row.names = F, append = TRUE)
  write.table(data.frame(anova_list[[i]]), "./Data Analysis Results/Anovas_Concat.csv", append = TRUE, sep = ",", col.names = TRUE, row.names = FALSE)
}

#### compile all spacing means #####
means <- list.files("./Data Analysis Results/", full.names = T)
means <- means[grepl("spacing_2", means)]

means_list <- list() 

for (i in 1:length(means)) {
  means_list[[i]] <- read.csv(means[i])
}

for (i in 1:length(means_list)) {
  write.table(str_split(means[i], "/")[[1]][10], "./Data Analysis Results/Means_Concat.csv", row.names = F, append = TRUE)
  write.table(data.frame(means_list[[i]]), "./Data Analysis Results/Means_Concat.csv", append = TRUE, sep = ",", col.names = TRUE, row.names = FALSE)
}


#### extrapolating yield - what would happen? ####
yield22 <- read.csv("./Data Analysis Results/mean_predicted_dry_cone_weight_g_spacing_22.csv") %>% mutate(year = 2022)
yield23 <- read.csv("./Data Analysis Results/mean_predicted_dry_cone_weight_g_spacing_23.csv") %>% mutate(year = 2023)
yield <- rbind(yield22, yield23)

# preparing the variables needed to evaluate this
yield <- yield %>% mutate(row_width_ft = 14, # 14 feet between plant rows
                 area_occupied_ft = spacing * row_width_ft, 
                 yield_g_per_sq_ft = response * area_occupied_ft, 
                 sq_ft_in_acre = 43560, 
                 plants_per_acre = sq_ft_in_acre/area_occupied_ft, 
                 yield_per_acre_g = plants_per_acre*response, 
                 yield_per_acre_lbs_adj = yield_per_acre_g / 453, 
                 yield_per_acre_lbs_non_adj = (response*889) / 453, # number of plants per acre
                 yield_per_hectare_kg_adj = yield_per_acre_lbs_adj * 1.1209, 
                 yield_per_hectare_kg_non_adj = yield_per_acre_lbs_non_adj * 1.1209)

yield <- yield %>% dplyr::select(spacing, year, yield_per_hectare_kg_adj, yield_per_hectare_kg_non_adj, .group) 

yield_long <- yield %>% pivot_longer(cols = c(yield_per_hectare_kg_adj, yield_per_hectare_kg_non_adj), names_to = "yield_type", values_to = "yield_kg_ha") %>%
  mutate(group = gsub(" ", "", .group))

non <- yield_long %>% filter(yield_type == "yield_per_hectare_kg_non_adj")
adj <- yield_long %>% filter(yield_type == "yield_per_hectare_kg_adj") %>% mutate(group = NA)

yield_long <- rbind(non, adj)

ggplot(data = yield_long, aes(x = as.factor(spacing), y = yield_kg_ha, colour = as.factor(year))) + 
         geom_point(stat="identity", aes(shape = yield_type), size = 5) +
  geom_text(aes(label = group), nudge_y = 300, size = 5) + 
  scale_color_manual(values = c("#009E73", "#CC79A7")) + 
  scale_shape_discrete(labels = c("Adjusted", "Non-Adjusted")) + 
  theme_bw() + 
  labs(y = "Yield (kg/ha)",
       x = "Spacing", shape = "Yield Type", color = "Year") +
  theme(legend.text = element_text(size = 20), 
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.y = element_text(size = font_size), 
        axis.text.x = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size)) 

#1200 x 900

#### how do our alpha acid estimates compare ####
alpha22 <- em_all %>% mutate(year = 2022)
alpha23 <- em_all %>% mutate(year = 2023) # acquire these from earlier in the script

alpha <- rbind(alpha22, alpha23) %>% filter(spacing == "3.5")

# publication records
pub_record <- data.frame(genotype = c("Comet", "Bitter Gold", "Magnum", "Mount Hood", "Mount Rainier", "Sterling", "Vanguard"), pub_alpha = c(median(c(8,11)), 17, median(c(12, 14)), 8, median(c(7,10)), median(c(6, 12)), median(c(4,6))))

alpha <- left_join(alpha, pub_record, by = "genotype")
alpha22 <- alpha %>% filter(year == 2022)
alpha23 <- alpha %>% filter(year == 2023)

cor.test(alpha22$emmean, alpha22$pub_alpha, method = "spearman")
cor.test(alpha23$emmean, alpha23$pub_alpha, method = "spearman")

plot(alpha22$emmean, alpha22$pub_alpha)
plot(alpha23$emmean, alpha23$pub_alpha)


#### spearman rank correlations figure ####
cor <- read.csv("./Data Analysis Results/Spearman Rank Correlations Table.csv", na = c("NA", ".", ""))
# this table was manually created from running all the spearman tests - see files
cor <- cor %>% mutate(r = as.numeric(r), spacing = as.factor(spacing))

cor$sig <- NA

for (i in 1:nrow(cor)) {
  if (! is.na(cor$p[i])) {
    if (cor$p[i] <= 0.1) {
      cor$sig[i] <- "yes"
      }
    else (cor$sig[i] <- "no")
  }
  else {cor$sig[i] <- NA }
}

cor <- cor %>% mutate(sig = as.factor(sig))

unique(cor$trait)

ggplot(cor) + 
  geom_col(width=0.7, color = "black", aes(x = spacing, y = r, fill = sig)) + 
  facet_grid(year ~ factor(trait, levels = c("Vigor", "Lateral Length", "Dry Cone Yield", "Harvest Index", "Alpha Acid", 
                                             "Beta Acid", "Cone Area", "Cone Density", "Cone Openness"))) + 
  coord_flip() + 
  scale_fill_manual(values = c("#9edcff", "#0072B2")) + 
  scale_y_continuous(limits = c(-0.25, 1.1), breaks = c(-0.25, 0.0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() + geom_hline(yintercept = 0) +
  labs(y = "Spearman Rank Correlation Coefficient",
       x = "Spacing (m)", fill = "Significance (α = 0.1)") +
  theme(legend.text = element_text(size = 20), 
        plot.subtitle = element_text(size = 15), 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = 15), 
        strip.text.y = element_text(size = font_size), 
        axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size),
        panel.spacing.x = unit(0.75, "lines"), 
        aspect.ratio = 2/1, 
        legend.position = "bottom") + 
  geom_text(data = sig_by_trait, mapping = aes(x = 0.54, y = 0.425, label = frequency), size = 4, fontface = "bold") + 
  geom_text(data = sig_by_year, mapping = aes(x = spacing, y = 1.02, label = frequency), size = 4, fontface = "bold")

### size = 2000 x 1000 TIFF

#### determine the frequency of significance for traits and spacings ####
cor_wide <- cor %>% dplyr::select(-sig) %>% pivot_wider(names_from = spacing, values_from = c(r, p))

sig_by_trait <- cor %>% 
  filter(p <= 0.1) %>% 
  group_by(trait, year) %>% 
  tally() %>% 
  mutate(n_env = 4, frequency = round(n/4, 2)) 

sig_by_trait <- rbind(sig_by_trait, data.frame(trait = "Dry Cone Yield", year = 2023, n = 0, n_env = 4, frequency = 0))

# adjust for missing data for vigor
sig_by_trait[17,5] <- 0.5

sig_by_year <- cor %>% 
  filter(p <= 0.1) %>% 
  group_by(spacing, year) %>% 
  tally() %>% 
  mutate(n_traits_years = 9, frequency = round(n/9, 2)) %>%
  mutate(trait = "Cone Openness")

sig_by_year[2, 4] <- 8
sig_by_year[4, 4] <- 8

sig_by_year <- sig_by_year %>% mutate(frequency = round(n / n_traits_years, 2))

# summarizing across single hill spacings
cor %>% filter(! spacing %in% c(0.46, 0.61)) %>% filter(sig == "yes") %>% summarise(mean = mean(r, na.rm = T))

# summarizing across seedling spacing
# which seedling env was most predictive
cor %>% filter(spacing %in% c(0.46, 0.61)) %>% group_by(year, spacing) %>% summarise(mean = mean(r, na.rm = T))
cor %>% filter(spacing %in% c(0.46, 0.61)) %>% group_by(year, spacing, sig) %>% tally()
cor %>% filter(spacing %in% c(0.46, 0.61)) %>% summarise(mean_r = mean(r, na.rm = T))
cor %>% filter(spacing %in% c(0.46, 0.61)) %>% group_by(spacing, year) %>% summarise(mean_r = mean(r, na.rm = T))
cor %>% filter(spacing %in% c(0.46, 0.61)) %>% group_by(year, sig) %>% tally()

# in 2022
8/18
# in 2023
10/16

# summarizing across all cone traits
cor_cones_22 <- cor %>% 
  filter(year == 2022) %>%
  filter(trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>% 
  filter(sig == "yes")

15/20 # 75% 
mean(cor_cones_22$r)

cor_cones_23 <- cor %>% 
  filter(year == 2023) %>%
  filter(trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>% 
  filter(sig == "yes")

17/20 # 85%
mean(cor_cones_23$r)

# summarizing across agronomic traits 
cor_agronomic_22 <- cor %>% 
  filter(year == 2022) %>%
  filter(!trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>% 
  filter(sig == "yes")
9/16
mean(cor_agronomic_22$r)

cor_agronomic_23 <- cor %>% 
  filter(year == 2023) %>%
  filter(!trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>% 
  filter(sig == "yes")
6/16
mean(cor_agronomic_23$r)

# 2022 for seedlings
cor_seedling_cones_22 <- cor %>% 
  filter(year == 2022) %>%
  filter(spacing %in% c(0.46, 0.61)) %>% 
  filter(trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>% 
  filter(sig == "yes")
  
nrow(cor_seedling_cones_22) / 10 # 0.7 frequency of significance
mean(cor_seedling_cones_22$r) # r = 0.91

# 2023
cor_seedling_cones_23 <- cor %>% 
  filter(year == 2023) %>%
  filter(spacing %in% c(0.46, 0.61)) %>% 
  filter(trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>% 
  filter(sig == "yes")

nrow(cor_seedling_cones_23) / 10 # 0.7 frequency of significance
mean(cor_seedling_cones_23$r) # r = 0.80

# agronomic traits
# 2022
cor_seedling_agronomic_22 <- cor %>% 
  filter(year == 2022) %>%
  filter(spacing %in% c(0.46, 0.61)) %>% 
  filter(! trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>%
  filter(sig == "yes")

nrow(cor_seedling_agronomic_22) / 8 # 0.125 frequency of significance 
mean(cor_seedling_agronomic_22$r) # r = 0.87

# 2023
cor_seedling_agronomic_23 <- cor %>% 
  filter(year == 2023) %>%
  filter(spacing %in% c(0.46, 0.61)) %>% 
  filter(! trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>%
  filter(sig == "yes")

nrow(cor_seedling_agronomic_23) / 8 # 0.375 frequency of significance 
mean(cor_seedling_agronomic_23$r) # r = 0.79

#### single-hill spacings
# 2022
cor_single_hill_agronomic_22 <- cor %>% 
  filter(year == 2022) %>%
  filter(! spacing %in% c(0.46, 0.61)) %>% 
  filter(! trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>%
  filter(sig == "yes")
# frequency of significance 
7 / 8 # 88%

# 2023
cor_single_hill_agronomic_23 <- cor %>% 
  filter(year == 2023) %>%
  filter(! spacing %in% c(0.46, 0.61)) %>% 
  filter(! trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>%
  filter(sig == "yes")
# frequency of significance 
3/8 # 37.5%

# 2022
cor_single_hill_cone_22 <- cor %>% 
filter(year == 2022) %>%
  filter(! spacing %in% c(0.46, 0.61)) %>% 
  filter(trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>%
  filter(sig == "yes")
# frequency of significance 
9/10 # 90%

# 2023
cor_single_hill_cone_23 <- cor %>% 
  filter(year == 2023) %>%
  filter(! spacing %in% c(0.46, 0.61)) %>% 
  filter(trait %in% c("Alpha Acid", "Beta Acid", "Cone Area", "Cone Density", "Cone Openness")) %>%
  filter(sig == "yes")
# frequency of significance 
10/10 # 100%

category <- data.frame(trait = unique(cor$trait), category = c("Cone", "Cone", "Cone", "Cone", "Cone", "Agronomic", "Agronomic", "Agronomic", "Agronomic"))
environemnt <- data.frame(spacing = unique(cor$spacing), environment = c("seedling", "seedling", "single hill", "single hill"))

cor <- cor %>% left_join(category, by = "trait")
cor <- cor %>% left_join(environemnt, by = "spacing")

r <- cor %>% group_by(year, category, environment) %>% summarise(mean_r = mean(r, na.rm = T))
n <- cor %>% group_by(year, category, environment) %>% tally() # total
sig <- cor %>% group_by(year, category, environment) %>% filter(sig == "yes") %>% tally()
sig$total <- n$n
sig <- sig %>% mutate(frequency_of_sig = n / total)

sig$rho <- r$mean_r

sig_seedlings <- sig %>% filter(environment == "seedling") %>% summarise(mean_sig = mean(frequency_of_sig), mean_rho = mean(rho))
mean(sig_seedlings$mean_rho)
cor %>% group_by(category, year) %>% summarise(mean_r = mean(r, na.rm = T))
cor_sig_year_type <- cor %>% filter(sig == "yes") %>% group_by(category, year) %>% tally()
cor_sig_year_type2 <- cor %>% group_by(category, year) %>% tally() #total
cor_sig_year_type$freq <-  cor_sig_year_type$n/cor_sig_year_type2$n
write.csv(cor, "./Data Analysis Results/correlation_summary.csv", row.names = F)

#### exploring why beta acid predictive ability goes away ####
# 2022
spacing_trial_dat_mean <- spacing_trial_dat_mean %>% mutate(UV_BA_adj = (mean_UV_BA/100) * mean_hop_box_weight_g_per_cone)
lmer_mean_UV_BA_adj_23 <- lmer(UV_BA_adj ~ genotype * spacing + (1|rep), data = spacing_trial_dat_mean %>% filter(year == "2023"))
resid_panel(lmer_mean_UV_BA_adj_23)
anova(lmer_mean_UV_BA_adj_23)

# effect of spacing
em_spacing <- emmeans(lmer_mean_UV_BA_adj_23, ~ spacing, type = "response")
cld <- data.frame(cld(em_spacing,  adjust = "none", Letters = c(letters))) 

# determine the gold standard order
em_all <- data.frame(emmeans(lmer_mean_UV_BA_adj_23, ~ genotype * spacing, type = "response"))
ord <- em_all %>% filter(spacing == 3.5) %>% arrange(emmean) %>% dplyr::select(genotype)
em_all2 <- arrange(transform(em_all,
                             genotype=factor(genotype,levels=neworder)),genotype)

dat_cor <- em_all %>% dplyr::select(genotype, spacing, emmean) %>% pivot_wider(names_from = spacing, values_from = emmean) 
cor.test(x = dat_cor$`3.5`, y = dat_cor$`1.5`,  method = "spearman") # r = 0.43; p = 0.3536
cor.test(x = dat_cor$`3.5`, y = dat_cor$`2`,  method = "spearman") # r = 0.77; p = 0.1028
cor.test(x = dat_cor$`3.5`, y = dat_cor$`5`,  method = "spearman") # r = 0.6; p = 0.1667
cor.test(x = dat_cor$`3.5`, y = dat_cor$`7`,  method = "spearman") # r = 0714; p = 0.088
cor %>% filter(spacing %in% c(0.46, 0.61)) %>% filter(trait == "Beta Acid") %>% filter(year == 2023)


               