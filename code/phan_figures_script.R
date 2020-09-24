#### Figures for perfusion manuscript

# Figure 1 Schematic for the perfusion wound meat model with stacked barchart of microbiome
# Figure 2 Perfusion and static model growth comparisons for the MB control community across antimicrobial treatments
# Figure 3 Growth, in CFU/mL, of MB alone, P. aeruginosa strains [PA14, PaFLR01] and their corresponding MB communities [MB(PA14), MB(PaFLR01)] in the perfused media meat model
# Figure 4 Total ion intensity of all volatile compounds. The PA14+MB community represents the microbial community on a meat sample with both PA14 and MB(PA14) populations. 
# Figure 5 Nonmetric multidimensional scaling (NMDS) plot of volatile signatures from the perfusion model from A) hydrogen peroxide conditions and the control and B) the antibiotic conditions and the control. 
# Figure 6 Heatmap showing clustering of relative abundances of volatile metabolites. PA14+MB represents the microbial community on a meat sample with both PA14 and MB(PA14) populations
# Figure 7 Relative abundances of metabolites significantly different between PA14+MB and MB control communities.  

# Supplemental Figure 2 Growth summary of additional treatments for the perfusion model with MB control, PA14 and MB(PA14) populations
# Supplemental Figure 3 Metabolites significantly different in the PA14+MB and MB control communities across treatments
# Supplemental Figure 4 MIC Pseudomonas PA14 and PaFLR01

########################################################################################
########################################################################################
# Figure 1 : stacked barcharts
# Schematic for the perfusion wound meat model with stacked barchart of microbiome

library(tidyverse)
library(reshape)
library(stringr)

otu <- read.csv("level-7.csv", row.names = 1)
otu <- otu[,colSums(otu) > 0]

## stacked barplot Figure 1D

otu$sample <- rownames(otu)
otu.melt <- melt(otu)
otu.melt$sample <- factor(otu.melt$sample, levels = c('MB1staticDNA', 'MB2StaticDNA', 'MB3staticDNA', 'MB1plateDNA', 'MB2PlateDNA', 'MB3PlateDNA'))

cbPalette <- c("#56B4E9","#CC79A7","#999999", "#F0E442", "#0072B2", "#D55E00", '#000000', "#009E73", "#E69F00",
               "#009E73","#E69F00","#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#version sideways
ggplot(subset(otu.melt, sample %in% c('MB1staticDNA', 'MB2StaticDNA', 'MB3staticDNA'))) +
  geom_bar(aes(x = sample, y = value, fill = variable), stat = 'identity', position = 'fill') + 
  theme_bw(base_size = 20) +
  theme(legend.position = 'bottom', 
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(fill=guide_legend(nrow=10,byrow=TRUE)) +
  scale_fill_manual(values = cbPalette) +
  labs(y = 'Relative Abundance') +
  coord_flip()

########################################################################################
########################################################################################
# Figure 2 : barchart of static vs perfusion
# Perfusion and static model growth comparisons for the MB control community across antimicrobial treatments

library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

avg <- read.csv("Perfusion_CFUs_summary_avg2.csv")
cfus.2 <- subset(avg, Condition %in% c("Control","Carbenicillin","Gentamicin", "H2O2_327mM") & 
                   Strain %in% c("PA14", "MB_PA14", "MB control", "PaFLR01", "MB_PaFLR01"))
cfus.2.melt <- melt(cfus.2, id.vars = c("Date", "Strain", "Condition", "Model"))
cfus.2.melt$Strain <- ordered(cfus.2.melt$Strain, levels = c("PA14", "MB_PA14","MB control", "PaFLR01", "MB_PaFLR01"))
cfus.2.melt$Condition <- ordered(cfus.2.melt$Condition, levels = c('Control', 'Carbenicillin', 'Gentamicin','H2O2_327mM'))

#comparing the perfusion and static models
ttest_table <-compare_means(value~Model, data = cfus.2.melt, group.by = c("Condition", 'Strain'), method = "t.test", p.adjust.method = "BH")

static_perfusion <- c("#000000", "#999999")
static_perfusion2 <- c("#999999", "#0072B2")

#t.tests comparing MB control for perfusion and static models
ttests_fig2 <-compare_means(value~Model, data = subset(cfus.2.melt, Strain == 'MB control'), group.by = c("Condition"), method = "t.test", p.adjust.method = "BH")

## anova for comparing models
# normal distribution of residuals
residuals <- resid(aov(value~Model, data = subset(cfus.2.melt, Strain == 'MB control' & Condition == 'H2O2_327mM')))
shapiro.test(residuals)
#control W = 0.90703, p-value = 0.3336
#Carbenicillin W = 0.96186, p-value = 0.8176
#Gentamicin W = 0.8925, p-value = 0.3315
#W = 0.9573, p-value = 0.7376

# homogeneity of variances
library(car)
leveneTest(value~Model, data = subset(cfus.2.melt, Strain == 'MB control' & Condition == 'H2O2_327mM'))
#Control 0.3699
#Carbenicillin 0.5189
#Gentamicin 0.646
#H2O2_327mM 0.4933

# anova for each treatment condition
anova <- aov(value~Model, data = subset(cfus.2.melt, Strain == 'MB control' & Condition == 'H2O2_327mM'))
summary(anova)
# control 0.0281 *
# Carbenicillin 2.66e-05 ***
# Gentamicin 0.00185 **
# H2O2_327mM 0.0226 *

# Figure 2 comparing MB in perfusion and static models
colblindpal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data = subset(cfus.2.melt, Strain == 'MB control'), 
       aes(x = Condition, y = value, color = Model, shape = Model)) + 
  geom_boxplot(position = "identity", fill = NA) + 
  geom_point(size = 3, position = position_jitter(width = 0.2)) +
  #geom_jitter(size = 3) +
  scale_shape_manual(values = c(16,6)) +
  scale_color_manual(values = colblindpal) +
  #facet_wrap(Condition~.) + 
  labs(x = "", y = "log10 CFU/ml") +
  theme_bw(base_size = 25) +
  theme(aspect.ratio = 0.5, 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #theme(axis.text.x = element_text(angle = 90, h = 1, v = .5)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  ylim(7, 11)

########################################################################################
########################################################################################
# Figure 3 : barchart of perfusion and treatments 
# Growth, in CFU/mL, of MB alone, P. aeruginosa strains [PA14, PaFLR01] and their corresponding MB communities [MB(PA14), MB(PaFLR01)] in the perfused media meat model

avg <- read.csv("Perfusion_CFUs_summary_avg2.csv")
avg1 <- avg %>%
  subset(., Model == "Perfusion" & 
           Strain %in% c("PA14", "MB_PA14", "MB control", 'PaFLR01', 'MB_PaFLR01') & 
           Condition %in% c("Control", "Carbenicillin", "Gentamicin", 
                            "H2O2_327mM")) %>%
  mutate(Condition = ordered(.$Condition, 
                             levels = c("Control", "Carbenicillin", "Gentamicin", "Gentamicin_H2O2_327mM", "H2O2_50uM",
                                        "H2O2_100uM","H2O2_200uM", "H2O2_500uM", "H2O2_1mM", "H2O2_10mM", 
                                        "H2O2_327mM", "H2O2_980mM")),
         Strain = ordered(.$Strain, levels = c("MB control","PA14", "MB_PA14", "PaFLR01", "MB_PaFLR01","phz", "MB_phz")))

#########################################
# anova tests with post-hoc tukey test
# http://www.sthda.com/english/wiki/one-way-anova-test-in-r
# assumptions for anova
## 1. normal distribution of residuals
## 2. homogeneity of variance 

##### treatment: Control, Carbenicillin, Gentamicin, H2O2_327mM 

# do for each treatment: check for normal distribution of residuals
residuals <- resid(aov(Avg ~ Strain, data = subset(avg1, Condition == 'H2O2_327mM')))
shapiro.test(residuals) 
# -> all residuals are normally distributed 

# do for each treatment: check for homogeneity of variances
library(car)
leveneTest(Avg ~ Strain, data = subset(avg1, Condition == 'Gentamicin'))
# H0 accepted - variances are homogeneous 

# anova (both assupmtions are met for Control, Carbenicillin, Gentamicin
aov <- aov(Avg ~ Strain, data = subset(avg1, Condition == 'H2O2_327mM'))
summary.aov(aov)
TukeyHSD(aov)

# Figure 3 
colblindpal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data = avg1, aes(x = Strain, y = Avg, col = Strain)) + 
  geom_boxplot(aes(color = Strain)) + 
  geom_point(position = "jitter", size = 3, aes(color = Strain)) + 
  #scale_color_manual(values = PA14_3) +
  scale_color_manual(values = colblindpal) +
  facet_wrap(~Condition, nrow = 2) + 
  labs(x = "", y = "log10 CFU/ml") +
  theme_bw(base_size = 20) +
  theme(aspect.ratio = .5, legend.position = 'na') +
  stat_compare_means(comparisons = list1, method = "t.test", label = "p.signif") + 
  stat_compare_means(comparisons = list2, method = "t.test", label = "p.signif") +
  ylim(6, 11)

########################################################################################
########################################################################################
# Figure 4 : dotplot of TIC
# Total ion intensity of all volatile compounds. The PA14+MB community represents the microbial community on a meat sample with both PA14 and MB(PA14) populations. 

tic <- read.csv("perfusion_headspace_tic.csv", row.names = "Compound", header = T)
tic_map <- read.csv("perfusion_headspace_tic_map.csv") %>%
  subset(., Strain %in% c('PA14+MB', 'MB control'))

column_sums <- colSums(tic)

column_sums.df <- as.data.frame(column_sums) %>%
  merge(., tic_map, by.x = 'row.names', by.y = 'Sample') %>%
  mutate(date = .$Row.names)
column_sums.df$date <- apply(column_sums.df, 1, function(x){
  return(unlist(strsplit(as.character(x["date"]), "_"))[1])
})

column_sums.df$date <- gsub("X", "", column_sums.df$date)
column_sums.df$Concentration <- ordered(column_sums.df$Concentration, levels = c('Control', '20ug/ml','50ug/ml', '20ug/ml + 327mM', '50uM', '100uM', '200uM', '327mM', '980mM'))
column_sums.df$Treatment <- ordered(column_sums.df$Treatment, levels = c('Control', 'Carbenicillin', 'Gentamicin', 'Gent + H2O2', 'H2O2'))


variable_labels <- c('Control', 'Carbenicillin', 'Gentamicin', "Gent + H2O[2]", 'H[2]O[2]')
column_sums.df$Strain <- ordered(column_sums.df$Strain, levels = c('PA14+MB', 'MB control', 'phz'))

#ttests
list1 <- list(c("PA14+MB", "MB control"))
compare_means(column_sums~Strain, data = column_sums.df, group.by = "Condition", method = "t.test")

PA14_MBcontrol <- c("#0072B2", "#E69F00")
#rename strain column to Model
colnames(column_sums.df)[3] <- 'Model'

ggplot(data = subset(column_sums.df, Model %in% c("PA14+MB", "MB control")), 
             aes(x = Concentration, y = column_sums, group = Concentration)) + 
  facet_grid(~Treatment, scales = "free") +
  scale_color_manual(values = PA14_MBcontrol) +
  scale_shape_manual(values = c(19,17)) +
  geom_point(stat = "identity", aes(shape = Model, color = Model), size = 4) + 
  ylab("Total ion intensity") +
  stat_compare_means(aes(group = Concentration), comparisons = list1, method = 't.test', label = 'p.signif') +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, h = 1, size = 14, vjust = 0.5), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16), 
        strip.text.x = element_text(size = 13))

### correlation between hydrogen peroxide concentration and TIC for PA14+MB
library(ggpubr)

h2o2_cor <- subset(column_sums.df, H2O2_conc_uM <327000 & Model == 'PA14+MB')
ggscatter(h2o2_cor, x = 'H2O2_conc_uM', y = 'column_sums', 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "H2O2", ylab = "TIC")

cor(h2o2_cor$column_sums, h2o2_cor$H2O2_conc_uM, method = 'pearson')


########################################################################################
########################################################################################
# Figure 5 : NMDS of antibiotcs and hydrogen peroxide 
# Nonmetric multidimensional scaling (NMDS) plot of volatile signatures from the perfusion model from A) hydrogen peroxide conditions and the control and B) the antibiotic conditions and the control. 

#https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/

library(vegan)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)
library(ggrepel)

#this analysis excludes the phz mutant and the data is normalized by metabolite intensity
tic_map <- read.csv("perfusion_headspace_tic_map.csv") %>%
  subset(., Strain %in% c('PA14+MB', 'MB control'))

#df below has rownames as metabolites so that the names stay as they should
#data is normalized by total sample intensity
tic_samplenorm <- read.csv("perfusion_headspace_tic_normalized2.csv", row.names = 1)
#df below changes orientation appropriate for nmds with rows as samples and columns as variables
tic_samplenorm <- as.data.frame(t(tic_samplenorm))

tic.nmds <- metaMDS(tic_samplenorm, distance = "bray", k = 2, autotransform = F, noshare = F)
tic.nmds
#Call:
#metaMDS(comm = tic_samplenorm, distance = "bray", k = 2, autotransform = F, noshare = F) 
#global Multidimensional Scaling using monoMDS
#Data:     tic_samplenorm 
#Distance: bray 
#Dimensions: 2 
#Stress:     0.1774659  
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘tic_samplenorm’

#A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions, < 0.1 is great, < 0.2 is good/ok, and stress < 0.3 provides a poor representation.**

sites <- as.data.frame(scores(tic.nmds, display = "sites"))
species <- as.data.frame(scores(tic.nmds, display = "species"))

nmds.sites <- sites %>%
  merge(x = sites, y = tic_map, by.x = 'row.names', by.y = "Sample") %>%
  mutate(Strain = ordered(.$Strain, levels = c("PA14+MB", "MB control")))

colnames(nmds.sites) <- c("Sample", "NMDS1", "NMDS2", "Model", "Condition", "Concentration", 'Treatment')
nmds.sites$Condition <- ordered(nmds.sites$Condition, levels = c('Control', 'Carbenicillin', 'Gentamicin', 'Gent_H2O2_327mM', 'H2O2_50uM', 'H2O2_100uM', 'H2O2_200uM', 'H2O2_327mM', 'H2O2_980mM'))

cbPalette <- c("#56B4E9", "#D55E00", "#0072B2", "#F0E442", "#009E73", "#CC79A7", "#E69F00", "#999999", 'black')
nmds_all <- ggplot(data = nmds.sites, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(color = Condition, shape = Model), size = 4) +
  scale_color_manual(values = cbPalette)

###############################################
######separate NMDS plots for hydrogen peroxide conditions + control and antibiotic conditions + control
###### hydrogen peroxide

tic_map.h2o2 <- subset(tic_map, Treatment %in% c("Control", "H2O2"))

tic.h2o2 <- tic_samplenorm %>% 
  add_rownames(., "Sample") %>%
  subset(., Sample %in% tic_map.h2o2$Sample) %>%
  column_to_rownames(., var = "Sample")

tic.nmds.h2o2 <- metaMDS(tic.h2o2, distance = "bray", k = 2, autotransform = F, noshare = F)
tic.nmds.h2o2
#Call:
#metaMDS(comm = tic.h2o2, distance = "bray", k = 2, autotransform = F,      noshare = F) 
#global Multidimensional Scaling using monoMDS
#Data:     tic.h2o2 
#Distance: bray 
#Dimensions: 2 
#Stress:     0.1546693 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘tic.h2o2’ 

sites.h2o2 <- as.data.frame(scores(tic.nmds.h2o2, display = "sites"))
species.h2o2 <- as.data.frame(scores(tic.nmds.h2o2, display = "species"))

nmds.sites.h2o2 <- sites.h2o2 %>%
  merge(x = sites.h2o2, y = tic_map.h2o2, by.x = 'row.names', by.y = "Sample") %>%
  mutate(Condition = ordered(.$Condition, levels = c('Control', 'H2O2_50uM', 'H2O2_100uM', 
                                                     'H2O2_200uM','H2O2_327mM', 'H2O2_980mM')),
         Strain = ordered(.$Strain, levels = c("PA14+MB", "MB control")))

colnames(nmds.sites.h2o2) <- c("Sample", "NMDS1", "NMDS2", "Strain", "Condition", "strain_by_condition")
cbPalette <- c("#56B4E9", "#D55E00", "#0072B2", "#F0E442", "#009E73", "#CC79A7", "#E69F00", "#999999")

nmds.plot1 <- ggplot() + 
  geom_point(data = nmds.sites.h2o2, aes(x = NMDS1, y = NMDS2, color = Strain, shape = Condition), size = 7, stroke = 2) + 
  #geom_text(data = species.h2o2, aes(NMDS1, NMDS2, label = rownames(species.h2o2)), 
  #         check_overlap = TRUE, size = 6, color = 'grey50') +
  annotate("text", x = -2, y = 0.8, label = "stress = 0.15", size = 6) +
  scale_shape_manual(values = c(19,0,17,4,1,15)) +
  scale_color_manual(values = c("#0072B2", "#E69F00")) +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  xlim(-3, 2)

## IQR interquartile range of NMDS2 for the hydrogen peroxide NMDS
# control
nmds.control.h2o2 <- subset(nmds.sites.h2o2, Treatment == 'Control' & Strain == 'PA14+MB')
IQR(nmds.control.h2o2$NMDS2) # 0.4195963
nmds.control.h2o2 <- nmds.control.h2o2 %>%
  mutate(out_IQR = if_else(abs(NMDS2) > 1.5*0.4195963, 'true', 'false'))

######separate NMDS plots for hydrogen peroxide conditions + control and antibiotic conditions + control
#### abx

tic_map.abx <- subset(tic_map, Treatment %in% c("Control", "Gentamicin", "Gentamicin + H2O2", "Carbenicillin"))

tic.abx <- tic_samplenorm %>% 
  add_rownames(., "Sample") %>%
  subset(., Sample %in% tic_map.abx$Sample) %>%
  column_to_rownames(., var = "Sample")

tic.nmds.abx <- metaMDS(tic.abx, distance = "bray", k = 2, autotransform = F, noshare = F)
tic.nmds.abx
#Call:
#metaMDS(comm = tic.abx, distance = "bray", k = 2, autotransform = F,      noshare = F) 
#global Multidimensional Scaling using monoMDS
#Data:     tic.abx 
#Distance: bray 
#Dimensions: 2 
#Stress:     0.09300182 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘tic.abx’

sites.abx <- as.data.frame(scores(tic.nmds.abx, display = "sites"))
species.abx <- as.data.frame(scores(tic.nmds.abx, display = "species"))

nmds.sites.abx <- sites.abx %>%
  merge(x = sites.abx, y = tic_map.abx, by.x = 'row.names', by.y = "Sample") %>%
  mutate(Treatment = ordered(.$Treatment, levels = c('Control', 'Carbenicillin', 'Gentamicin', 
                                                     'Gentamicin + H2O2')),
         Strain = ordered(.$Strain, levels = c("PA14+MB", "MB control")))

colnames(nmds.sites.abx)[5] <- "Treatment"
colnames(nmds.sites.abx)[7] <- "Condition"

nmds.plot2 <- ggplot() + 
  geom_point(data = nmds.sites.abx, aes(x = NMDS1, y = NMDS2, shape = Condition, color = Strain), size = 7, stroke = 2) + 
  #geom_text(data = species.abx, aes(NMDS1, NMDS2, label = rownames(species.abx)), 
  #         check_overlap = TRUE, size = 6, color = 'grey50') +
  annotate("text", x = -1, y = 1.7, label = "stress = 0.15", size = 6) +
  scale_shape_manual(values = c(19,0,15,2)) +
  scale_color_manual(values = c("#0072B2", "#E69F00")) +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))

### two panel figures
## 1555 x 550
ggarrange(nmds.plot1, nmds.plot2, labels = c('A', 'B'), widths = c(1,1), heights = c(1,1), align = 'h')
ggsave("test.pdf",plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 40, height = 15, units = c("cm"),
       dpi = 300, limitsize = TRUE)

grid.newpage()
print(proportion_ordered_final, vp = viewport(x = 0.31, y = 0.5, width = 0.6, height = .95))
print(dendro.plot, vp = viewport(x = 0.81, y = 0.513, width = 0.4, height = .88))

#cowplot
ggdraw() +
  draw_plot(nmds.plot1, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(nmds.plot2, x = 0.5, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1))

########################################################################################
########################################################################################
# Figure 6 : heatmap of volatiles
#Heatmap showing clustering of relative abundances of volatile metabolites. PA14+MB represents the microbial community on a meat sample with both PA14 and MB(PA14) populations

library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggdendro)
library(grid)
library(vegan)

###dendrogram help: https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html

#mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
#scale_fill_gradientn(colours = mycol)

tic_map <- read.csv("perfusion_headspace_tic_map.csv")
tic_samplenorm <- read.csv('perfusion_headspace_tic_normalized2.csv', row.names = 1)

###heatmap
tic_samplenorm.hm <- tic_samplenorm %>%
  rownames_to_column(var = 'Sample') %>%
  melt(id.vars = 'Sample') %>%
  rename(., 'Compound' = 'Sample', 'Sample' = 'variable') %>%
  merge(., tic_map, by = "Sample") %>%
  subset(., Strain %in% c("PA14+MB", "MB control"))
colnames(tic_samplenorm.hm) <- c('Sample', 'Compound', 'value', 'Strain', 'Condition', 'Concentration', 'Treatment')

tic_samplenorm.hm$Condition <- ordered(tic_samplenorm.hm$Condition, levels = c('Control','Carbenicillin' ,'Gentamicin', 'Gent_H2O2_327mM', 'H2O2_50uM', 'H2O2_100uM', 'H2O2_200uM', 'H2O2_327mM', 'H2O2_980mM'))

heatmap <- ggplot() + geom_tile(data = tic_samplenorm.hm, aes(x = Strain, y = Compound, fill = log2(value))) + 
  theme(axis.text.x = element_text(size = 12, angle = 90), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12), axis.title.y = element_blank(),
        strip.text.x = element_text(size = 12, angle = 90)) + 
  #scale_fill_gradientn(colors = mycol) +
  scale_fill_continuous(low = "#132B43", high = "#56B1F7", na.value = "#132B43") +
  facet_grid(~Condition, scales = "free")

###dendrogram
tic_samplenorm.matrix <- vegdist(tic_samplenorm, method = "bray")
tic_samplenorm.matrix <- as.matrix(t(tic_samplenorm.matrix))

####rownames(tic_analyte.matrix) <- unique(tic_analyte.hm$Compound) ###this looks wrong
tic_samplenorm.dendro <- as.dendrogram(hclust(d = dist(x = tic_samplenorm.matrix)))
dendro.plot <- ggdendrogram(data = tic_samplenorm.dendro, rotate = TRUE, labels = TRUE, leaf_labels = FALSE) +
  theme(axis.text.y = element_text(hjust = 0), axis.text.x = element_text(angle = 90))

#####putting it all together
tic_samplenorm.order <- order.dendrogram(tic_samplenorm.dendro)
tic_samplenorm.hm$Compound <- factor(tic_samplenorm.hm$Compound, 
                                     levels = tic_samplenorm.hm$Compound[tic_samplenorm.order],
                                     ordered = TRUE)
##reordering what is in dendrogram to reflect in heatmap
## check the order of the dendrogram and the heatmap to make sure all of the metabolties are in the correct order. 
heatmap.reordered <- ggplot(tic_samplenorm.hm, aes(Strain, Compound, fill = log2(value))) +
  geom_tile(color = "gray", size = 0.20) +
  facet_grid(~Condition, scales = "free") +
  scale_fill_viridis_c(na.value = '#440154FF') +
  #scale_fill_continuous(low = "#132B43", high = "#56B1F7", na.value = "#132B43") +
  theme(legend.position = "top", 
        axis.text.y = element_blank(), # can hash at the beginning of there to check order of mx
        #axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 12, h = 0, vjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(), strip.text = element_text(size = 14))

#the values in viewport need some adjustmends based on figure sizes
grid.newpage()
print(heatmap.reordered, vp = viewport(x = 0.35, y = 0.5, width = 0.7, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.815, y = 0.475, width = 0.25, height = 0.865))

########################################################################################
########################################################################################
# Figure 7 : boxplot of metabolites different between MB and control and PA14+MB
# Relative abundances of metabolites significantly different between PA14+MB and MB control communities.  

library(reshape2)
library(tidyverse)
library(ggpubr)

tic_samplenorm <- read.csv("perfusion_headspace_tic_normalized2.csv", row.names = 1)
tic_map <- read.csv("perfusion_headspace_tic_map.csv")

### comparing all metabolites between MB and PA14+MB / across conditions 
tic_all <- tic_samplenorm %>%
  rownames_to_column() %>%
  melt(.) %>%
  merge(., tic_map, by.x = 'variable', by.y = 'Sample')

# comparing PA14+MB and MB control
tic_all_results <- compare_means(formula = value ~ Strain, 
                                 data = tic_all, 
                                 group.by = 'rowname',
                                 method = 't.test', p.adjust.method = 'BH')
write.csv(tic_all_results, 'mx_comparisons_PA14+MB_MBcontrol.csv')

# results: metabolites significantly different between PA14+MB and MB control
significant_mx <- subset(tic_all_results, p.adj < 0.05)

#rowname                   .y.   group1  group2                p     p.adj p.format p.signif method
#<chr>                     <chr> <chr>   <chr>             <dbl>     <dbl> <chr>    <chr>    <chr> 
#1 2-butanone                value PA14+MB MB control 0.00426      0.029     0.00426  **       T-test
#2 2-Nonanone                value PA14+MB MB control 0.000678     0.006     0.00068  ***      T-test
#3 Acetophenone              value PA14+MB MB control 0.0000296    0.00046   3.0e-05  ****     T-test
#4 Butanal, 2-methyl         value PA14+MB MB control 0.000345     0.0041    0.00035  ***      T-test
#5 Butanal, 3-methyl         value PA14+MB MB control 0.00819      0.043     0.00819  **       T-test
#6 Ethanone, 1-(2-aminoph... value PA14+MB MB control 0.0000000410 0.0000019 4.1e-08  ****     T-test
#7 Heptanal                  value PA14+MB MB control 0.000765     0.006     0.00077  ***      T-test
#8 Phenylethyl alcohol       value PA14+MB MB control 0.000000121  0.0000028 1.2e-07  ****     T-test
#9 Propanal, 2-methyl        value PA14+MB MB control 0.00799      0.043     0.00799  **       T-test

### data: subsetting mx data for significantly different metabolites between PA14+MB and MB control
mx <- subset(tic_samplenorm, rownames(tic_samplenorm) %in% significant_mx$rowname)
mx <- mx %>% 
  rownames_to_column() %>%
  melt(.) %>%
  merge(., tic_map, by.x = 'variable', by.y = 'Sample')

### figure 7
## plot showing relative abundances of these metabolites for each model (either with PA14 or MB)
figure_7 <- ggplot(data = mx, aes(x = Strain, y = value, col = Strain)) +
  geom_point(position = position_jitter(width = 0.2), size = 2) +
  scale_color_manual(values = c("#E69F00","#0072B2")) +
  geom_boxplot(fill = NA) +
  facet_wrap(~rowname, scales = 'free') +
  theme_bw(base_size = 20) +
  labs(y = 'Relative Abundance') +
  theme(axis.title.x = element_blank(), legend.position = 'na',
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = 't.test', label = 'p.signif')

## t test for each metabolite between MB and PA14+MB
ttests_figure7 <- compare_means(formula = value ~ Strain, 
                                data = mx, 
                                group.by = 'rowname',
                                method = 't.test', p.adjust.method = 'bonferroni')


########################################################################################
########################################################################################
# Supplemental Figure 2 : boxplot of additional perfusion treatments 
# Growth summary of additional treatments for the perfusion model with MB control, PA14 and MB(PA14) populations

avg <- read.csv("Perfusion_CFUs_summary_avg2.csv")
avg1 <- avg %>%
  subset(., Model == "Perfusion" & 
           Strain %in% c("PA14", "MB_PA14", "MB control")) %>%
  mutate(Condition = ordered(.$Condition, 
                             levels = c("Control", "Carbenicillin", "Gentamicin", "Gentamicin_H2O2_327mM", "H2O2_50uM",
                                        "H2O2_100uM","H2O2_200uM", "H2O2_500uM", "H2O2_1mM", "H2O2_10mM", 
                                        "H2O2_327mM", "H2O2_980mM")),
         Strain = ordered(.$Strain, levels = c("MB control","MB_PA14", "PA14", "phz", "MB_phz", "PaFLR01", 
                                               "_MBPaFLR01")))

list3 <- list(c("PA14", "MB_PA14"), c("MB control", "MB_PA14"))

#########################################
# anova tests with post-hoc tukey test
# http://www.sthda.com/english/wiki/one-way-anova-test-in-r
# assumptions for anova
## 1. normal distribution of residuals
## 2. homogeneity of variance 

##### treatment: Control, Carbenicillin, Gentamicin, Gentamicin_H2O2_327mM, H2O2_50uM, H2O2_100uM, H2O2_200uM, H2O2_500uM, H2O2_1mM, H2O2_10mM, H2O2_327mM, H2O2_980mM
# do for each treatment: check for normal distribution of residuals
residuals <- resid(aov(Avg ~ Strain, data = subset(avg1, Condition == 'H2O2_980mM')))
shapiro.test(residuals) 

# do for each treatment: check for homogeneity of variances
library(car)
leveneTest(Avg ~ Strain, data = subset(avg1, Condition == 'H2O2_980mM'))

# anova 
# Control, Carbenicillin, Gentamicin, Gentamicin_H2O2_327mM, H2O2_50uM, H2O2_100uM, H2O2_200uM, H2O2_500uM, H2O2_1mM, H2O2_10mM, H2O2_327mM, H2O2_980mM
aov <- aov(Avg ~ Strain, data = subset(avg1, Condition == 'Carbenicillin'))
summary.aov(aov)
TukeyHSD(aov)

# figure supplemental 1
ggplot(data = avg1, aes(x = Strain, y = Avg)) + 
  geom_boxplot(aes(color = Strain)) + 
  geom_point(position = "jitter", size = 3, aes(color = Strain)) + 
  scale_color_manual(values = PA14_3) +
  #scale_color_manual(values = c("#000033FF", "#1400FFFF", "#C729D6FF", "#FF9C63FF", "#FFFF60FF")) +
  facet_wrap(~Condition) + 
  labs(x = "", y = "log10 CFU/ml") +
  theme_bw(base_size = 20) +
  theme(legend.position = 'none') +
  stat_compare_means(comparisons = list3, method = "t.test", label = "p.signif") +
  ylim(6.5, 11.5)

########################################################################################
########################################################################################
# Supplemental Figure 3 : boxplot of diff abundant metabolites in MB control or PA14+MB across treatments 
# Metabolites significantly different in the PA14+MB and MB control communities across treatments

## continued from Figure 7 script 

mx$Treatment <- ordered(mx$Treatment, levels = c('Control', 'Carbenicillin','Gentamicin', 'Gent + H2O2', 'H2O2'))
list1 <- list(c('Control', 'Carbenicillin'))
list2 <- list(c('Control', 'Gentamicin'))
list3 <- list(c('Control', 'Gent + H2O2'))
list4 <- list(c('Control', 'H2O2'))

# selecting all conditions and only the h2o2 condition at concentration 327 mM
mx2 <- subset(mx, !(Treatment == "H2O2" & Concentration != "327mM"))

## t test for comparing metabolite across the treatment conditions
# cannot have levels
#PA14+MB ........'Ethanone, 1-(2-aminoph...', 'Acetophenone', '2-Nonanone', '2-butanone'
compare_means(formula = value ~ Treatment, 
              data = subset(mx2, Strain == 'PA14+MB' & rowname ==  '2-butanone'), 
              group.by = 'rowname',
              method = 't.test', p.adjust.method = 'BH', ref.group = 'Control')

#########################################
# anova tests with post-hoc tukey test
# http://www.sthda.com/english/wiki/one-way-anova-test-in-r
# assumptions for anova
## 1. normal distribution of residuals
## 2. homogeneity of variance 

#PA14+MB ........'Ethanone, 1-(2-aminoph...', 'Acetophenone', '2-Nonanone', '2-butanone'
# do for each treatment: check for normal distribution of residuals
residuals <- resid(aov(value ~ Treatment, data = subset(mx2, Strain == 'PA14+MB' & rowname ==  '2-Nonanone')))
shapiro.test(residuals) 


# do for each treatment: check for homogeneity of variances
library(car)
leveneTest(value ~ Treatment, data = subset(mx2, Strain == 'PA14+MB' & rowname ==  '2-butanone'))

### since not all did not meet assumptions: alternative is to perform krusall-wallice test and pairwise.wilcox.test
subset <- subset(mx2, Strain == 'PA14+MB' & rowname == 'Ethanone, 1-(2-aminoph...')
#kruskal.test(value ~ Treatment, data = subset)
#pairwise.wilcox.test(subset$value, subset$Treatment, p.adjust.method = 'BH')
compare_means(value ~ Treatment, data = subset, method = 't.test', p.adjust.method = 'BH', ref.group = 'Control')

#MB control ......'Phenylethyl alcohol', 'Butanal, 2-methyl', 'Butanal, 3-methyl', 'Heptanal', 'Propanal, 2-methyl'
compare_means(formula = value ~ Treatment, 
              data = subset(mx2, Strain == 'MB control' & rowname ==  'Propanal, 2-methyl'), 
              group.by = 'rowname',
              method = 't.test', p.adjust.method = 'BH', ref.group = 'Control')

suppfig2_pa14.mb <- ggplot(data = subset(mx2, rowname %in% c('Ethanone, 1-(2-aminoph...', 'Acetophenone', '2-Nonanone', '2-butanone') & Strain == 'PA14+MB'), aes(x = Treatment, y = value, col = Strain)) +
  scale_color_manual(values = c("#0072B2")) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(fill = NA) +
  facet_wrap(~rowname, scales = 'free', nrow = 1) +
  theme_bw(base_size = 16) +
  labs(y = 'Relative Abundance') +
  stat_compare_means(comparisons = list1, method = "t.test", label = "p.format") +
  stat_compare_means(comparisons = list2, method = "t.test", label = "p.format") +
  stat_compare_means(comparisons = list3, method = "t.test", label = "p.format") +
  stat_compare_means(comparisons = list4, method = "t.test", label = "p.format") +
  theme(legend.position = 'na', axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))

suppfig2_mb <- ggplot(data = subset(mx2, rowname %in% c('Phenylethyl alcohol', 'Butanal, 2-methyl', 'Butanal, 3-methyl', 'Heptanal', 'Propanal, 2-methyl') & Strain == 'MB control'), aes(x = Treatment, y = value, col = Strain)) +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_color_manual(values = c("#E69F00")) +
  geom_boxplot(fill = NA) +
  facet_wrap(~rowname, scales = 'free_y', nrow = 2, ncol = 4) +
  theme_bw(base_size = 16) +
  labs(y = 'Relative Abundance') +
  stat_compare_means(comparisons = list1, method = "t.test", label = "p.format") +
  stat_compare_means(comparisons = list2, method = "t.test", label = "p.format") +
  stat_compare_means(comparisons = list3, method = "t.test", label = "p.format") +
  stat_compare_means(comparisons = list4, method = "t.test", label = "p.format") +
  theme(legend.position = 'na', axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))

ggarrange(suppfig2_pa14.mb, suppfig2_mb, ncol = 1, nrow = 2, heights = c(1.4, 2.2), labels = c('A', 'B'))

########################################################################################
########################################################################################
## supplemental figure 4

library(reshape2)
library(ggplot2)
library(plyr)

df1 <- read.csv("MIC_summary.csv")
df1 <- melt(df1)

df1$Concentration <- ordered(df1$Concentration, levels = c('Control', '20 ug/ml', '50 ug/ml', '75 ug/ml', '100 ug/ml', '150 ug/ml', '200 ug/ml', '250 ug/ml', '300 ug/ml', '400 ug/ml'))

df1.mean <- df1 %>% 
  ddply(., .(Strain, Concentration, Treatment, Timepoint), transform, avg = mean(value, na.rm = T)) %>%
  select(., c(Strain, Concentration, Treatment, Timepoint, avg)) %>%
  unique(.)

colblindpal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df1$Treatment <- factor(df1$Treatment, levels = c('Control', 'Carbenicillin', 'Gentamicin'))
(p1 <- ggplot(data = df1, aes(x = Concentration, y = value, col = Treatment, shape = Timepoint)) + 
    geom_boxplot(aes(shape = Timepoint, col = Treatment), fill = NA) +
    geom_point(size = 3, position = position_jitterdodge(jitter.width = 0.1)) +
    facet_grid(Strain~.) +
    theme_bw(base_size = 16) +
    scale_color_manual(values = colblindpal) +
    theme(axis.text = element_text(size = 11), axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 14), strip.text = element_text(size = 14),
          legend.text = element_text(size = 12), legend.title = element_text(size = 14),
          panel.grid = element_blank()) +
    labs(y = "OD 500 nm"))
