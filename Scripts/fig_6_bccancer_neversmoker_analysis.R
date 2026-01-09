library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(patchwork)
library(broom)

# LUCAS: For the version we upload: please remove columns that are not used in analysis 
shade_data = fread("data/Anthracosis_Project_master_20251218.csv") %>%
  janitor::clean_names() %>%
  mutate(birth_region = ifelse(birth_continent%in%c("Asia"),birth_continent,"Other")) %>%
  mutate(birth_region = ifelse(birth_continent%in%c(""),NA,birth_region)) %>%
  mutate(birth_region = factor(birth_region,levels=c("Other","Asia"))) %>%
  mutate(birth_continent = ifelse(birth_continent%in%c(""),NA,birth_continent)) %>%
  mutate(egf_rmut=ifelse(egf_rmut=="",NA,egf_rmut)) %>%
  mutate(ethnicity=ifelse(ethnicity=="","Unknown",ethnicity))

# Removing any outliers in the data
shade_data_no_outlier = shade_data %>%
filter(anth_percent<0.4)

shade_for_analysis = shade_data_no_outlier

# Stating variables to be tested
numeric_vars = c("age","cum_3yr_pm","cum_30yr_pm")
cat_vars = c("sex","upper_vs_lower","birth_region","t2_ormore","egf_rmut")

# Comparing anthracosis levels with numerical variables (Spearman correlation)
num_res_df = lapply(numeric_vars,function(x){
  tt = shade_for_analysis %>%
    reshape::melt(id.vars=c("pt_id","anth_percent"),measure.vars=x)
  res = cor.test(tt$anth_percent,tt$value,method="spearman")
  c(x,as.numeric(res$estimate),res$p.value)
})
num_res_df = data.frame(do.call(rbind,num_res_df))
colnames(num_res_df) = c("variable","estimate","p.value")
num_res_df$estimate = as.numeric(num_res_df$estimate)
num_res_df$p.value = as.numeric(num_res_df$p.value)

# Comparing anthracosis levels with categorical variables (t test)
var_res_df = lapply(cat_vars,function(x){
  tt = shade_for_analysis %>%
    reshape::melt(id.vars=c("pt_id","anth_percent"),measure.vars=x)
  tt %>%
    t.test(anth_percent ~ value, data = .) %>%
    tidy() %>%
    mutate(variable=x) %>%
    dplyr::select(variable,estimate,p.value)
})
var_res_df = bind_rows(var_res_df)

# Combining all analyses and multiple test correction
all_res_df = bind_rows(num_res_df,var_res_df)
all_res_df$p.adj = p.adjust(all_res_df$p.value,method="BH")


# Plotting
birth_region_plot = shade_for_analysis %>%
  filter(!is.na(birth_region)) %>%
  ggplot(aes(x=birth_region,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom(aes(color=birth_continent)) +
  stat_compare_means(method="t.test") +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,
  ) +
  ylab("% Anthracosis") +
  xlab("Birthplace (Continent)")

eth_order = shade_for_analysis %>%
  group_by(ethnicity) %>%
  summarise(median_anth_percent=median(anth_percent)) %>%
  arrange(-median_anth_percent) %>%
  pull(ethnicity)

ethnicity_plot = shade_for_analysis %>%
  filter(!is.na(ethnicity)) %>%
  mutate(ethnicity=factor(ethnicity,levels=eth_order)) %>%
  mutate(birth_continent=ifelse(is.na(birth_continent),"Unknown",as.character(birth_continent))) %>%
  ggplot(aes(x=ethnicity,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom(aes(color=birth_continent)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,
  ) +
  ylab("% Anthracosis") +
  xlab("Ethnicity")

sex_plot = shade_for_analysis %>%
  ggplot(aes(x=sex,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom() +
  stat_compare_means(method="t.test") +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,        # Adjust the width of the crossbar (0 to 1)
  )+
  ylab("% Anthracosis") +
  xlab("Sex")

lobe_plot = shade_for_analysis %>%
  ggplot(aes(x=upper_vs_lower,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom() +
  stat_compare_means(method="t.test") +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,        
  )+
  ylab("% Anthracosis") +
  xlab("Lobe Sampled")

age_plot = shade_for_analysis %>%
  ggplot(aes(x=age,y=anth_percent)) +
  geom_point() +
  geom_smooth(method="lm",se = FALSE) +
  stat_cor(method="spearman") +
  theme_cowplot()+
  ylab("% Anthracosis") +
  xlab("Age")


pm25_30yr_plot = shade_for_analysis %>%
  mutate(birth_continent=ifelse(is.na(birth_continent),"Unknown",as.character(birth_continent))) %>%
  ggplot(aes(x=cum_30yr_pm,y=anth_percent)) +
  geom_point(aes(color=birth_continent)) +
  stat_cor(method="spearman") +
  theme_cowplot()+
  ylab("% Anthracosis") +
  xlab("30yr PM2.5 Exposure")

pm25_3yr_plot = shade_for_analysis %>%
  mutate(birth_continent=ifelse(is.na(birth_continent),"Unknown",as.character(birth_continent))) %>%
  ggplot(aes(x=cum_3yr_pm,y=anth_percent)) +
  geom_point(aes(color=birth_continent)) +
  stat_cor(method="spearman") +
  theme_cowplot()+
  ylab("% Anthracosis") +
  xlab("3yr PM2.5 Exposure")


subtype_order = shade_for_analysis %>%
  group_by(subtype) %>%
  summarise(median_anth_percent=median(anth_percent)) %>%
  arrange(median_anth_percent) %>%
  pull(subtype)

subtype_plot = shade_for_analysis %>%
  mutate(subtype=factor(subtype,levels=rev(subtype_order))) %>%
  ggplot(aes(x=subtype,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,        
  )+
  ylab("% Anthracosis") +
  xlab("Tumour Subtype")

egfr_plot = shade_for_analysis %>%
  filter(!is.na(egf_rmut)) %>%
  ggplot(aes(x=egf_rmut,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom() +
  stat_compare_means(method="t.test") +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,   
  )+
  ylab("% Anthracosis") +
  xlab("Tumour EGFR Mutation")

stage_plot = shade_for_analysis %>%
  ggplot(aes(x=staging,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom() +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,   
  )+
  ylab("% Anthracosis") +
  xlab("Tumour Stage")

#  Regression to identify relative associations of anthracosis levels across variables
all_res_df %>%
  filter(p.adj<=0.1)

summary(ols <- glm(anth_percent ~ age + sex + upper_vs_lower + birth_region, data = shade_for_analysis))

tidy_df <- tidy(ols, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(significance=ifelse(p.value<0.05,"p<0.05","None"))%>%
  mutate(significance=ifelse(p.value<0.01,"p<0.01",significance))


regression_plot = ggplot(tidy_df, aes(x = estimate, y = reorder(term, estimate))) +
  geom_point(size = 3,aes(color=significance)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "Effect estimate (95% CI)",
    y = NULL
  ) +
  theme_cowplot() +
  scale_color_manual(values=c("black","red3","salmon"))


#  MAIN FIGURE
(((age_plot | age_plot ) + plot_layout(widths=c(2,1))) / (sex_plot | lobe_plot | birth_region_plot ) / regression_plot ) +
  plot_layout(guides="collect",heights=c(2,2,1)) +
  plot_annotation(tag_levels = "A")

# SUPPLEMENTARY FIGURE
(
  ((pm25_3yr_plot | pm25_30yr_plot | ethnicity_plot) + 
    plot_layout(guides="collect", widths=c(1,1,2))) / 
    (( egfr_plot | stage_plot | subtype_plot) +
    plot_layout(guides="collect", widths=c(1,1,2)))
  ) +
  plot_annotation(tag_levels = "A") 
