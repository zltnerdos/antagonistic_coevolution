# ---------------------------
# Purpose of script: Process and analyse colony count data from timeshift experiment for the journal article "Community complexity does not weaken pairwise coevolution in a soil bacterial community" (Ecology Letters)
# This script is in addition to timeshift_ze.R. It runs the same process without excluding the outlier point 'VD4PC4'. The reason why this is a standalone script is that the model simplification is not automated, therefore model selection would be different if timeshift_ze were to be run on the dataset including 'VD4PC4'.
# What this script does:
# 1. Read in raw data from timeshift.csv
# 2. Convert variables
# 3. Visualise and remove outliers
# 4. Statistical analysis of contemporaneous interactions, coevolution, and total density
# 5. Visualise contemporaneous interactions, coevolution, and total density
#
# Author: Dr. Zoltan Erdos
#
# Date last modified:  2025-10-21
# -------------------------

# Load required packages
# -----------------------------------------------------------------------------
librarian::shelf(
  tidyverse, janitor, dplyr, patchwork, GGally, ggpubr, RColorBrewer, cowplot,
  lme4, emmeans, car, rstatix, DHARMa, boot, LaplacesDemon, devtools,
  BlakeRMills/MoMAColors, outliers, broom.mixed
)
# Set theme for plots
theme_ze <-  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        strip.background = element_blank(),
        strip.text = element_text(angle = 0, hjust = 0.0, size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))


#  import timeshift data
# -----------------------------------------------------------------------------

timeshift <- read.csv('data/timeshift.csv', header = T)

# Data type conversions and factor level ordering
# =============================================================================

# Convert categorical variables to factors
# -----------------------------------------------------------------------------
timeshift$treatment <- as.factor(timeshift$treatment)
timeshift$community <- as.factor(timeshift$community)
timeshift$rep <- as.factor(timeshift$rep)

# Convert numeric variables 
# -----------------------------------------------------------------------------
timeshift$dilution <- as.numeric(timeshift$dilution)  # Dilution factor
timeshift$P <- as.numeric(timeshift$P)                # P counts (raw)
timeshift$V <- as.numeric(timeshift$V)                # V counts (raw)
timeshift$productivity <- as.numeric(timeshift$productivity)
# include random effect for P and V rep separately
timeshift$repP <- interaction(timeshift$P_time, timeshift$rep)

timeshift$repV <- interaction(timeshift$V_time, timeshift$rep)

timeshift <- timeshift %>%
  mutate(
    repP = case_when(
      repP == "ancestor.1" ~ "ancestor.9",
      repP == "ancestor.2" ~ "ancestor.10", 
      repP == "ancestor.3" ~ "ancestor.11",
      repP == "ancestor.4" ~ "ancestor.12",
      repP == "ancestor.5" ~ "ancestor.13",
      repP == "ancestor.6" ~ "ancestor.14",
      repP == "ancestor.7" ~ "ancestor.15",
      repP == "ancestor.8" ~ "ancestor.16",
      TRUE ~ repP  # Keep all other values unchanged
    ),
    repV = case_when(
      repV == "ancestor.1" ~ "ancestor.9",
      repV == "ancestor.2" ~ "ancestor.10",
      repV == "ancestor.3" ~ "ancestor.11", 
      repV == "ancestor.4" ~ "ancestor.12",
      repV == "ancestor.5" ~ "ancestor.13",
      repV == "ancestor.6" ~ "ancestor.14",
      repV == "ancestor.7" ~ "ancestor.15",
      repV == "ancestor.8" ~ "ancestor.16",
      TRUE ~ repV  # Keep all other values unchanged
    )
  )

#convert to factor
timeshift$repP <- as.factor(timeshift$repP)
timeshift$repV <- as.factor(timeshift$repV)
# Set factor levels in chronological/logical order
# -----------------------------------------------------------------------------
# Time points for V measurements (chronological order)
timeshift$V_time <- factor(
  timeshift$V_time, 
  levels = c("ancestor", "6 weeks", "10 weeks")
) 

# Time points for P measurements (chronological order)
timeshift$P_time <- factor(
  timeshift$P_time, 
  levels = c("ancestor", "6 weeks", "10 weeks")
)

# Treatment levels (from simple to complex culture conditions)
timeshift$treatment <- factor(
  timeshift$treatment, 
  levels = c('ancestor', 'Mono', 'Coculture', 'Community'),
  labels = c('ancestor', 'monoculture', 'coculture', 'community')
)

# Remove outlier
#histogram of V counts to identify outliers
# -----------------------------------------------------------------------------
ggplot(timeshift, aes(x = V)) +
  geom_histogram(bins = 30, color = "black") +
  labs(title = "Histogram of V Counts", x = "V Counts", y = "Frequency") +
  theme_minimal()

#histogram of P counts to identify outliers
ggplot(timeshift, aes(x = P)) +
  geom_histogram(bins = 30, color = "black") +
  labs(title = "Histogram of P Counts", x = "P Counts", y = "Frequency") +
  theme_minimal()


#outlier seems to be VD4PC4, treatment = coculture, V_time = 10 weeks, P_time = 6 weeks
#create a subset for these conditions
outlier_condition <- timeshift %>%
  filter(treatment == 'coculture', V_time == '10 weeks', P_time == '6 weeks')
#make a plot showing the outlier
ggplot(timeshift, aes(x = community, y = V)) +
  geom_point() +
  labs(x = "all Variovorax replicates", y = "Variovorax count") +
  theme_minimal()

ggplot(outlier_condition, aes(x = community, y = V)) +
  geom_point() +
  labs(x = "replicate", y = "Variovorax count") +
  theme_minimal()

# Plot 1: All Variovorax replicates
Sup1 <- ggplot(timeshift, aes(x = community, y = V)) +
  geom_point(color = "Black", size = 2) +
  labs(
    x = "All Variovorax replicates",
    y = "Variovorax count"
  ) +
  theme_minimal(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Outlier condition
Sup2 <- ggplot(outlier_condition, aes(x = community, y = V)) +
  geom_point(color = "Black", size = 2) +
  labs(
    x = "Treatment level replicates",
    y = "Variovorax count"
  ) +
  theme_minimal(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined <- Sup1 + Sup2 +
  plot_annotation(
    tag_levels = "A",   # adds A, B labels automatically
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

combined

#grubbs test for outliers in V counts
grubbs.test(timeshift$V)
# OUTLIERS ARE NOT REMOVED
#timeshift <- timeshift %>%
  #filter(!(community == 'VD4PC4'))

# Calculate log-transformed cell densities per ml
# =============================================================================
# Formula: log10((raw_count * 10^dilution_factor) / sample_volume_ul)
# Note: Sample volume 35 μL

timeshift$Plogperml <- log10((timeshift$P * 10^timeshift$dilution) / 35)
timeshift$Vlogperml <- log10((timeshift$V * 10^timeshift$dilution) / 35)

#convert to long format for plots
time_long <- pivot_longer(timeshift, cols = c("P", "V"), names_to = "strain", values_to = "count")


# =============================================================================
# STATISTICAL ANALYSIS - CONTEMPORANEOUS INTERACTIONS
# =============================================================================

# Filter data for contemporaneous time points (P and V measured at same time)
# -----------------------------------------------------------------------------

contemporaneous_interactions <- timeshift %>%
  filter((P_time == V_time))

# P density 
# =============================================================================

# Fit initial mixed-effects model with interaction term
Pdesc <- lmer(Plogperml ~ treatment * P_time + (1|repP) + (1|repV),
              data = contemporaneous_interactions)

aa <- allFit(Pdesc)
print(aa)
summary(aa)
tidy(aa, conf.int = TRUE) |> 
  arrange(effect, term, estimate) |> 
  select(-c(std.error, statistic))

glance(aa) |> select(optimizer, AIC, NLL_rel) |> arrange(NLL_rel)

# Test significance of interaction term using likelihood ratio test
drop1(Pdesc, test = "Chisq") # no interaction

Pdesc <- update(Pdesc, ~ . - treatment:P_time)

drop1(update(Pdesc, REML = FALSE), test = "Chisq") # no P_time effect

Pdesc <- update(Pdesc, ~ . - P_time) 

drop1(update(Pdesc, REML = FALSE), test = "Chisq") # treatment effect


plot(simulateResiduals(fittedModel = Pdesc, plot = FALSE))
summary(Pdesc)

# Post-hoc comparisons: treatment effects across all time points
emPres <- emmeans(Pdesc, specs = pairwise ~ treatment)
print(emPres) #not sure if we can use this comparison...

#test equality of variances between co and community culture
# Levene's test for homogeneity of variances for only coculture and community

levene_test(Plogperml ~ treatment, data = subset(contemporaneous_interactions, treatment == c('coculture', 'community')))
var.test(Plogperml ~ treatment, data = subset(contemporaneous_interactions, treatment == c('coculture', 'community')))
# V density 
# =============================================================================

# Fit initial mixed-effects model with interaction term
Vdesc <- lmer(Vlogperml ~ treatment * V_time + (1|repP) + (1|repV), data = contemporaneous_interactions)

# Test significance of interaction term using likelihood ratio test
drop1(Vdesc, test = "Chisq")

Vdesc <- update(Vdesc, ~ . - treatment:V_time)

drop1(update(Vdesc, REML = FALSE), test = "Chisq")

Vdesc <- update(Vdesc, ~ . - V_time) # no V_time effect

drop1(update(Vdesc, REML = FALSE), test = "Chisq")

plot(simulateResiduals(fittedModel = Vdesc, plot = FALSE))

summary(Vdesc)

# Post-hoc comparisons: treatment effects within each time point
emVres <- emmeans(Vdesc, specs = pairwise ~ treatment)
print(emVres)

# =============================================================================
# ANALYSIS of P proportion | coevolution in each treatment
# =============================================================================

# proportions in monoculture
prop_pmonoculture <- glmer(cbind(P,V) ~ P_time*V_time+(1|repP) + (1|repV),
                            data = subset(timeshift, treatment %in% c('monoculture')),
                            family = binomial)

drop1(prop_pmonoculture, test = 'Chisq')

prop_pmonoculture <- update(prop_pmonoculture, ~ . - P_time:V_time) #no interaction
drop1(update(prop_pmonoculture), test = "Chisq")
prop_pmonoculture <- update(prop_pmonoculture, ~ . - V_time) #no V_time effect
drop1(update(prop_pmonoculture), test = "Chisq")

# Model diagnostics using DHARMa residual analysis
plot(simulateResiduals(fittedModel = prop_pmonoculture, plot = FALSE)) 

summary(prop_pmonoculture) # 6 week P higher proportion
# Post-hoc comparisons for P proportion in monoculture
emmeans(prop_pmonoculture, specs = pairwise ~ P_time)

# proportions in coculture
prop_pcoculture <- glmer(cbind(P,V) ~ P_time*V_time+(1|repP) + (1|repV),
                          data = subset(timeshift, treatment %in% c('coculture')),
                          family = binomial)

drop1(prop_pcoculture, test = 'Chisq') # significant interaction

plot(simulateResiduals(fittedModel = prop_pcoculture, plot = FALSE))

summary(prop_pcoculture)
emmeans(prop_pcoculture, specs = pairwise ~ P_time)

# proportions in Community
prop_pcommunity <- glmer(cbind(P,V) ~ P_time*V_time+(1|repP) + (1|repV),
                          data = subset(timeshift, treatment %in% c('community')),
                          family = binomial)
drop1(prop_pcommunity, test = 'Chisq')
plot(simulateResiduals(fittedModel = prop_pcommunity, plot = FALSE))

summary(prop_pcommunity)
emmeans(prop_pcommunity, specs = pairwise ~ P_time | V_time)

# =============================================================================
# ANALYSIS of P proportion | does coevolution differ between coculture and community?
# =============================================================================

prop_p.1noa <-  glmer(cbind(P,V) ~ P_time*V_time*treatment+ (1|repP) + (1|repV),
                      data = subset(timeshift, !(treatment %in% c('monoculture','ancestor'))),
                      family = binomial)
drop1(prop_p.1noa, test = 'Chisq') # 3way interaction

# Model diagnostics using DHARMa residual analysis
plot(simulateResiduals(fittedModel = prop_p.1noa, plot = FALSE))
summary(prop_p.1noa)
#post-hoc comparisons for P proportion in coculture and community
coandcom_compare <- emmeans::emmeans(prop_p.1noa, pairwise ~ treatment)

# check optimizers
noanc <- allFit(prop_p.1noa)
print(noanc)
summary(noanc)
tidy(noanc, conf.int = TRUE) |> 
  arrange(effect, term, estimate) |> 
  select(-c(std.error, statistic))

glance(aa) |> select(optimizer, AIC, NLL_rel) |> arrange(NLL_rel)

#comparison of P density between coculture and community
P2.1cc <- lmer(Plogperml~V_time*P_time*treatment+ (1|repP) + (1|repV),
               data = subset(timeshift, !(treatment %in% c('monoculture','ancestor'))))
               
drop1(P2.1cc, test = 'Chisq')
P2.1cc <- update(P2.1cc, ~ . - P_time:V_time:treatment)
drop1(update(P2.1cc, REML = 'F'), test = "Chisq")
P2.1cc <- update(P2.1cc, ~ . - P_time:treatment)
drop1(update(P2.1cc, REML = 'F'), test = "Chisq")
P2.1cc <- update(P2.1cc, ~ . - V_time:treatment)
drop1(update(P2.1cc, REML = 'F'), test = "Chisq")
P2.1cc <- update(P2.1cc, ~ . - V_time:P_time)
drop1(update(P2.1cc, REML = 'F'), test = "Chisq")
P2.1cc <- update(P2.1cc, ~ . - P_time) # no P_time effect
drop1(update(P2.1cc, REML = 'F'), test = "Chisq")
P2.1cc <- update(P2.1cc, ~ . - treatment) # no treatment effect
summary(P2.1cc)
plot(simulateResiduals(fittedModel =P2.1cc, plot = F))
pcount_full<-emmeans(P2.1cc, specs = pairwise ~ V_time)

#comparison of V density between coculture and community
V2.1cc <- lmer(Vlogperml~V_time*P_time*treatment+ (1|repP) + (1|repV),
               data = subset(timeshift, !(treatment %in% c('monoculture','ancestor'))))

drop1(V2.1cc, test = 'Chisq')
V2.1cc <- update(V2.1cc, ~ . - P_time:V_time:treatment)
drop1(update(V2.1cc, REML = 'F'), test = "Chisq")
V2.1cc <- update(V2.1cc, ~ . - V_time:treatment)
drop1(update(V2.1cc, REML = 'F'), test = "Chisq")
V2.1cc <- update(V2.1cc, ~ . - P_time:treatment)
drop1(update(V2.1cc, REML = 'F'), test = "Chisq")
summary(V2.1cc)
plot(simulateResiduals(fittedModel =V2.1cc, plot = F))
vcount_full<-emmeans(V2.1cc, specs = pairwise ~ treatment)


# =============================================================================
# ANALYSIS of total density
# =============================================================================

density <-  lmer(log10(productivity) ~ P_time*V_time*treatment+ (1|repP) + (1|repV), data = (timeshift))
drop1(density, test = 'Chisq')
density <- update(density, ~ . - P_time:V_time:treatment) # no 3way interaction
drop1(density, test = 'Chisq')                 
density <- update(density, ~ . - V_time:treatment) # no interaction
drop1(density, test = 'Chisq')
density <- update(density, ~ . - P_time:treatment) # no interaction
drop1(density, test = 'Chisq')
density <- update(density, ~ . - treatment) # no treatment effect
drop1(density, test = 'Chisq')
density <- update(density, ~ . - P_time:V_time) # no P_time:V_time interaction
drop1(density, test = 'Chisq')
density <- update(density, ~ . - P_time) # no P_time effect
drop1(density, test = 'Chisq')
# Model diagnostics using DHARMa residual analysis
plot(simulateResiduals(fittedModel = density, plot = FALSE))
summary(density)

density_comp <- emmeans(density, list(pairwise ~ V_time), infer = TRUE)
summary(density_comp)
emmeans(density, specs = pairwise ~ V_time)
# =============================================================================
# Figures
# =============================================================================

#dummy dataset to get ancestors for plot
time_long_anc <- time_long %>%
  filter(P_time  %in% c('ancestor') & V_time %in% c('ancestor'))

# Create a summary dataframe with mean values
summary_data_P <- timeshift %>%
  group_by(treatment, P_time, V_time) %>%
  summarize(mean = mean(P_prop,na.rm = TRUE),
            mean_countP = mean(P,na.rm = TRUE),
            maxProp = max(P_prop, na.rm = TRUE),
            minProp = min(P_prop, na.rm = TRUE),
            sd = sd(P_prop, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n))

time_long$treatment <- factor(time_long$treatment, levels = c('ancestor','monoculture', 'coculture', 'community'))

# create a dataset with the ancestor values duplicated for each treatment (to have as reference point)
# Subset the ancestor rows
ancestor_rows <- subset(timeshift, treatment == "ancestor")
# Define the target treatments to duplicate to
new_treatments <- c("monoculture", "coculture", "community")
# Duplicate ancestor rows and assign new treatment labels
duplicated_ancestors <- do.call(rbind, lapply(new_treatments, function(treat) {
  temp <- ancestor_rows
  temp$treatment <- treat
  return(temp)
}))
# Combine original data with duplicated ancestor data
timeshift_with_duplicated_ancestors <- rbind(timeshift, duplicated_ancestors)
timeshift_with_duplicated_ancestors <- subset(timeshift_with_duplicated_ancestors, treatment != ('ancestor'))
time_long_pseudo_anc <- pivot_longer(timeshift_with_duplicated_ancestors, cols = c("P", "V"), names_to = "strain", values_to = "count")

# Create a summary dataframe with mean values
summary_data_pseudo_P <- timeshift_with_duplicated_ancestors %>%
  group_by(treatment, P_time, V_time) %>%
  summarize(mean = mean(P_prop,na.rm = TRUE),
            mean_logcountP = mean(log10(P),na.rm = TRUE),
            maxProp = max(P_prop, na.rm = TRUE),
            minProp = min(P_prop, na.rm = TRUE),
            sd = sd(P_prop, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            countsd = sd(log10(P), na.rm = TRUE),
            countse = countsd / sqrt(n))

# Create a summary dataframe with mean values
summary_data_pseudo_V <- timeshift_with_duplicated_ancestors %>%
  group_by(treatment, V_time, P_time) %>%
  summarize(mean = mean(V_prop,na.rm = TRUE),
            mean_logcountV = mean(log10(V),na.rm = TRUE),
            maxProp = max(V_prop, na.rm = TRUE),
            minProp = min(V_prop, na.rm = TRUE),
            sd = sd(V_prop, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            countsd = sd(log10(V), na.rm = TRUE),
            countse = countsd / sqrt(n))

# Plots for density comparison
# change to data to contemporaneous_interactions_start_end for only anc and 10 week
contemporaryplot_P <- ggplot(contemporaneous_interactions, 
                             aes(x = P_time, y = log10(P*10^dilution/35))) +
  geom_boxplot(aes(fill = treatment), 
               outlier.shape = NA) +
  geom_point(aes(fill = treatment), size = 2, alpha = 0.2, 
             position = position_jitterdodge(jitter.width = 0.1)) +
  scale_fill_moma_d(palette = "VanGogh")+ 
  ylab(expression('density CFU/ml '(log[10]))) +
  scale_x_discrete(name = NULL, 
                   labels = function(x) sub(".","\n",x,fixed=TRUE)) +
  ggtitle("A) Pseudomonas") +
  theme_ze +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

contemporaryplot_V <- ggplot(contemporaneous_interactions, 
                             aes(x = V_time, y = log10(V*10^dilution/35))) +
  geom_boxplot(aes(fill = treatment), 
               outlier.shape = NA) +
  geom_point(aes(fill = treatment), size = 2, alpha = 0.2, 
             position = position_jitterdodge(jitter.width = 0.1)) +
  scale_fill_moma_d(palette = "VanGogh")+ 
  ylab(expression('density CFU/ml '(log[10]))) +
  scale_x_discrete(name = NULL, 
                   labels = function(x) sub(".","\n",x,fixed=TRUE)) +
  ggtitle("B) Variovorax") +
  theme_ze +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank()
  )

contemporary_total <- ggplot(contemporaneous_interactions, 
                             aes(x = V_time, y = log10(productivity*10^dilution/35))) +
  geom_boxplot(aes(fill = treatment), 
               outlier.shape = NA) +
  geom_point(aes(fill = treatment), size = 2, alpha = 0.2, 
             position = position_jitterdodge(jitter.width = 0.1)) +
  scale_fill_moma_d(palette = "VanGogh")+ 
  ylab(expression('density CFU/ml '(log[10]))) +
  scale_x_discrete(name = NULL, 
                   labels = function(x) sub(".","\n",x,fixed=TRUE)) +
  ggtitle("C) Total density") +
  theme_ze +
    theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank()
  )

Fig2 <- contemporaryplot_P + contemporaryplot_V + contemporary_total +
  plot_layout(ncol = 3, nrow = 1, guides = 'collect') 
Fig2
ggsave("Fig2.pdf", height = 8, width = 10)
ggsave("Fig2.tiff",
       width = 6,          # in inches
       height = 6,         # in inches
       dpi = 600,          # recommended for publication
       units = "in",
       compression = "lzw")  # lossless compression



###colour by Ptime
timeshift_pseudo_anc_ordered <- timeshift_with_duplicated_ancestors[order(timeshift_with_duplicated_ancestors$V_time),]

Fig3 <- ggplot(subset(timeshift_pseudo_anc_ordered, treatment  %in% c('monoculture', 'coculture', 'community')), aes(x = V_time, y = V_prop, group = interaction(P_time, rep)))+
  geom_point(aes(x = V_time, y = V_prop, group = interaction(P_time, rep), colour = P_time), size = 2, alpha = 0.2, position = position_dodge(width = 0.06))+
  geom_point(data = subset(time_long_anc, strain == 'V' & treatment != 'ancestor'),
             aes(x = V_time, y = V_prop, group = rep), colour = "#2a2e38", size = 2, alpha = 0.2, position = position_dodge(width = 0.06))+
  geom_path(aes(x = V_time, y = V_prop, group = interaction(P_time, rep), colour = P_time),
            position = position_dodge(width = 0.06),
            alpha = 0.2,
            linewidth = 0.3)+
  scale_colour_moma_d(palette = "VanGogh")+
  ylab('Variovorax proportion')+
  scale_x_discrete(name = "Variovorax timepoint",labels=function(x) sub(".","\n",x,fixed=TRUE))+
  ylim(c(0,1))+
  theme_ze+
  theme(strip.text.y.right = element_blank())+
  facet_grid(~treatment, drop = TRUE,
             labeller = labeller(treatment = c("monoculture" = "A) Monoculture", 
                                               "coculture" = "B) Coculture", 
                                               "community" = "C) Community"),
                                 multi_line = TRUE))+
  geom_errorbar(inherit.aes = FALSE, data = summary_data_pseudo_V,
                aes(x = V_time, ymin = mean-se, ymax = mean+se, group = P_time),
                width = 0.15,
                position = position_dodge(width = 0.08))+
  geom_line(data = summary_data_pseudo_V, aes(x = V_time, y = mean, colour = P_time, group = P_time),
            size = 2,
            position = position_dodge(width = 0.06))+
  geom_point(inherit.aes = FALSE, data = summary_data_pseudo_V,
             aes(x = V_time, y = mean, colour = P_time),
             size = 4,
             position = position_dodge(width = 0.08))+
  guides(colour=(guide_legend(title = 'Pseudomonas')))
Fig3
ggsave("Fig3.pdf", height = 8, width = 10)
ggsave("Fig3.tiff",
       width = 10,          # in inches
       height = 8,         # in inches
       dpi = 600,          # recommended for publication
       units = "in",
       compression = "lzw")  # lossless compression

###colour by vtime
timeshift_pseudo_anc_ordered <- timeshift_with_duplicated_ancestors[order(timeshift_with_duplicated_ancestors$P_time),]

Fig4 <- ggplot(subset(timeshift_pseudo_anc_ordered, treatment  %in% c('monoculture', 'coculture', 'community')), aes(x = P_time, y = P_prop, group = interaction(V_time, rep)))+
  geom_point(aes(x = P_time, y = P_prop, group = interaction(V_time, rep), colour = V_time), size = 2, alpha = 0.2, position = position_dodge(width = 0.06))+
  geom_point(data = subset(time_long_anc, strain == 'P' & treatment != 'ancestor'),
             aes(x = P_time, y = P_prop, group = rep), colour = "#2a2e38", size = 2, alpha = 0.2, position = position_dodge(width = 0.06))+
  geom_path(aes(x = P_time, y = P_prop, group = interaction(V_time, rep), colour = V_time),
            position = position_dodge(width = 0.06),
            alpha = 0.2,
            linewidth = 0.3)+
  scale_colour_moma_d(palette = "VanGogh")+
  ylab('Pseudomonas proportion')+
  scale_x_discrete(name = "Pseudomonas timepoint",labels=function(x) sub(".","\n",x,fixed=TRUE))+
  ylim(c(0,1))+
  theme_ze+
  theme(strip.text.y.right = element_blank())+
  facet_grid(~treatment, drop = TRUE,
             labeller = labeller(treatment = c("monoculture" = "A) Monoculture", 
                                               "coculture" = "B) Coculture", 
                                               "community" = "C) Community"),
                                 multi_line = TRUE))+
  geom_errorbar(inherit.aes = FALSE, data = summary_data_pseudo_P,
                aes(x = P_time, ymin = mean-se, ymax = mean+se, group = V_time),
                width = 0.15,
                position = position_dodge(width = 0.08))+
  geom_line(data = summary_data_pseudo_P, aes(x = P_time, y = mean, colour = V_time, group = V_time),
            size = 2,
            position = position_dodge(width = 0.06))+
  geom_point(inherit.aes = FALSE, data = summary_data_pseudo_P,
             aes(x = P_time, y = mean, colour = V_time),
             size = 4,
             position = position_dodge(width = 0.08))+
  
  guides(colour=(guide_legend(title = 'Variovorax')))
Fig4
ggsave("Fig4.pdf", height = 8, width = 10)
ggsave("Fig4.tiff",
       width = 10,          # in inches
       height = 8,         # in inches
       dpi = 600,          # recommended for publication
       units = "in",
       compression = "lzw")  # lossless compression

#post-hoc comparisons for total density
density_emm <- emmeans(density, specs = pairwise ~ V_time)
print(density_emm)

#plot for total density
total_density <- time_long %>%
  group_by(V_time, treatment, rep) %>%  # include replicate if needed to keep boxplot spread
  summarise(total_count = sum(count, na.rm = TRUE),
            dilution = unique(dilution)) %>%  # assume same dilution per group
  mutate(log_density = log10(total_count * 10^dilution / 35),
         strain = "Pseudomonas+Variovorax")  # synthetic strain label for plotting

#Fig5 <- ggplot(subset(time_long, treatment %in% c('monoculture', 'coculture', 'community')), aes(x = V_time, y = log10(productivity*10^dilution/35)))+
#  geom_boxplot(aes(x=V_time, y=log10(count*10^dilution/35), fill=strain), outliers = F, alpha = 0.2)+
#  geom_boxplot(aes(fill = "Pseudomonas+Variovorax"), outliers = F, alpha = 0.8)+
#  geom_point(aes(x=V_time, y=log10(count*10^dilution/35), fill=strain), size = 2, alpha = 0.2, position = position_jitterdodge(jitter.width = 0.1))+
##  scale_fill_moma_d(palette = "VanGogh", 
#                    labels = c("Pseudomonas", "Pseudomonas+Variovorax", "Variovorax"))+
#  ylab(expression('density CFU/ml '(log[10])))+
#  scale_x_discrete(name = 'evolutionary time',labels=function(x) sub(".","\n",x,fixed=TRUE))+
# theme_ze+
#  theme(legend.title=element_blank())+
#  facet_grid(~treatment, drop = TRUE,
#             labeller = labeller(treatment = c("monoculture" = "A) Monoculture", 
#                                              "coculture" = "B) Coculture", 
#                                               "community" = "C) Community"),
#                                 multi_line = TRUE))
#Fig5
#ggsave("Fig5.pdf", height = 8, width = 10)
#ggsave("Fig5.tiff",
#      width = 10,          # in inches
#       height = 8,         # in inches
#       dpi = 600,          # recommended for publication
#       units = "in",
#       compression = "lzw")  # lossless compression

