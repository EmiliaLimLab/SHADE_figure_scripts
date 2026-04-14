library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(patchwork)
library(broom)

# Never-Smoking Cohort - FIGURE 4 Analysis
shade_data = fread("../Cohort_data/Never_smoking_Cohort.csv") %>%
  janitor::clean_names() %>%
  mutate(birth_region = ifelse(birth_continent%in%c("Asia"),birth_continent,"Other")) %>%
  mutate(birth_region = ifelse(birth_continent%in%c(""),NA,birth_region)) %>%
  mutate(birth_region = factor(birth_region,levels=c("Other","Asia"))) %>%
  mutate(birth_continent = ifelse(birth_continent%in%c(""),NA,birth_continent)) %>%
  mutate(egf_rmut=ifelse(egf_rmut=="",NA,egf_rmut)) %>%
  mutate(ethnicity=ifelse(ethnicity=="","Unknown",ethnicity)) %>%
  mutate(subtypeGroup=ifelse(subtype%in%c("Micropapillary","Solid","Acinar","Papillary","Cribriform","Squamous"),"HighIntermediateGrade",subtype)) %>%
  mutate(subtypeGroup=ifelse(subtype%in%c("Lepidic","Mucinous"),"LowGrade",subtypeGroup))

shade_for_analysis = shade_data

# Stating variables to be tested
numeric_vars = c("age","cum_3yr_pm","cum_30yr_pm")
cat_vars = c("sex","upper_vs_lower","birth_region", "egf_rmut","subtypeGroup")

# Comparing anthracosis levels with numerical variables (Spearman correlation)
num_res_df = lapply(numeric_vars,function(x){
  tt = shade_for_analysis %>%
    reshape::melt(id.vars=c("study_id","anth_percent"),measure.vars=x)
  res = cor.test(tt$anth_percent,tt$value,method="spearman")
  c(x,as.numeric(res$estimate),res$p.value)
})
num_res_df = data.frame(do.call(rbind,num_res_df))
colnames(num_res_df) = c("variable","estimate","p.value")
num_res_df$estimate = as.numeric(num_res_df$estimate)
num_res_df$p.value = as.numeric(num_res_df$p.value)

# Comparing anthracosis levels with categorical variables (wilcox test)
var_res_df = lapply(cat_vars,function(x){
  tt = shade_for_analysis %>%
    reshape::melt(id.vars=c("study_id","anth_percent"),measure.vars=x)
  tt %>%
    wilcox.test(anth_percent ~ value, data = .) %>%
    tidy() %>%
    mutate(variable=x) %>%
    dplyr::select(variable,statistic,p.value)
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
  stat_compare_means(method="wilcox.test") +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,
  ) +
  ylab("% Anthracosis") +
  xlab("Birthplace (Continent)") +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(0.98, 0.02),   # bottom-right
        legend.justification = c(1, 0))

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
  stat_compare_means(method="wilcox.test") +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,        # Adjust the width of the crossbar (0 to 1)
  )+
  ylab("% Anthracosis") +
  xlab("Sex") +
  coord_cartesian(clip = "off")

lobe_plot = shade_for_analysis %>%
  ggplot(aes(x=upper_vs_lower,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom() +
  stat_compare_means(method="wilcox.test", size = 4.5) +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,        
  )+
  ylab("% Anthracosis") +
  xlab("Lobe Sampled") +
  coord_cartesian(clip = "off")

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
  geom_point() +
  geom_smooth(method="lm",se = FALSE) +
  stat_cor(method="spearman", label.x.npc = "left", label.y.npc = 1,
           vjust = 0.40) +
  theme_cowplot()+
  ylab("% Anthracosis") +
  xlab(expression("30yr " * PM[2.5] * " Exposure")) +
  coord_cartesian(clip = "off")

pm25_3yr_plot = shade_for_analysis %>%
  mutate(birth_continent=ifelse(is.na(birth_continent),"Unknown",as.character(birth_continent))) %>%
  ggplot(aes(x=cum_3yr_pm,y=anth_percent)) +
  geom_point() +
  geom_smooth(method="lm",se = FALSE) +
  stat_cor(method="spearman", label.y.npc = 1, vjust = 0.40, size = 4.5) +
  theme_cowplot()+
  ylab("% Anthracosis") +
  xlab(expression("3yr " * PM[2.5] * " Exposure")) +
  coord_cartesian(clip = "off")


subtype_order = shade_for_analysis %>%
  group_by(subtype) %>%
  summarise(median_anth_percent=median(anth_percent)) %>%
  arrange(median_anth_percent) %>%
  pull(subtype)

subtype_plot = shade_for_analysis %>%
  mutate(subtype=factor(subtype,levels=rev(subtype_order))) %>%
  ggplot(aes(x=subtype,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom(aes(color=subtypeGroup)) +
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

subtype_group_plot = shade_for_analysis %>%
  filter(subtypeGroup%in%c("HighIntermediateGrade","LowGrade")) %>%
  ggplot(aes(x=subtypeGroup,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom() +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,        
  )+
  stat_compare_means(method="wilcox.test") +
  ylab("% Anthracosis") +
  xlab("Tumour Subtype") +
  theme(axis.text.x = element_text(size = 10))+
  coord_cartesian(clip = "off")


egfr_plot = shade_for_analysis %>%
  filter(!is.na(egf_rmut)) %>%
  ggplot(aes(x=egf_rmut,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom() +
  stat_compare_means(method="wilcox.test", size = 4.5) +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,   
  )+
  ylab("% Anthracosis") +
  xlab("Tumour EGFR Mutation") +
  coord_cartesian(clip = "off")

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

pleura_plot = shade_for_analysis %>%
  filter(!is.na(pleura_present)) %>%
  ggplot(aes(x=pleura_present,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom() +
  stat_compare_means(method="wilcox.test", size = 4.5) +
  theme_cowplot() +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar", color="blue",
    width = 0.7,   
  )+
  ylab("% Anthracosis") +
  xlab("Pleura Present")

#  Regression to identify relative associations of anthracosis levels across variables
all_res_df %>%
  filter(p.adj<=0.1)

summary(ols <- glm(anth_percent ~ age + sex + birth_region + subtypeGroup + cum_30yr_pm, data = shade_for_analysis))

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
  scale_color_manual(values=c("black","red3","salmon")) +
  theme(legend.position = c(0.98, 0.02),   # bottom-right
        legend.justification = c(1, 0))


#  MAIN FIGURE
main_figure <- (((age_plot | age_plot ) + plot_layout(widths=c(2,1))) / 
    ((sex_plot | subtype_group_plot | pm25_30yr_plot ) + plot_layout(widths=c(1,1,1))) /
    (( birth_region_plot + regression_plot) + plot_layout(widths=c(1,2))))   +
  plot_layout(guides="collect",heights=c(2,2,2)) +
  plot_annotation(tag_levels = "A")

# SUPPLEMENTARY FIGURE
suppl_figure <- (((pm25_3yr_plot | egfr_plot | ethnicity_plot) + 
    plot_layout(guides="collect", widths=c(1,1,2))) / 
    ((stage_plot | lobe_plot | pleura_plot) +
       plot_layout(guides="collect", widths=c(1,1,2))) +
      (subtype_plot +
         plot_layout(guides="collect", widths=c(1,1,2)))
) +
  plot_annotation(tag_levels = "A") 

#  Smoking Cohort - FIGURE 5 Analysis
ever_smoker_data_df = fread("data/Smoking_Cohort.csv") %>%
  janitor::clean_names() %>%
  mutate(smoking_status = smoking) %>%
  select(slide_id,age,sex,pack_year,smoking_status,upper_vs_lower,anth_percent)

never_smoker_df = shade_data_no_outlier %>% 
  mutate(pack_year=NA) %>%
  select(slide_id,age,sex,pack_year,smoking_status,upper_vs_lower,anth_percent)

comb_df = bind_rows(ever_smoker_data_df,never_smoker_df)

ever_never_smoker_plot = comb_df %>%
  mutate(smoking_status=ifelse(smoking_status%in%c("current","former"),"ever","never")) %>%
  ggplot(aes(x=smoking_status,y=anth_percent)) +
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
  xlab("Smoking Status") +
  stat_compare_means()


pack_year_plot = ever_smoker_data_df %>%
  ggplot(aes(x=pack_year,y=anth_percent)) +
  geom_point() +
  stat_cor(method="spearman") +
  geom_smooth(method="lm") +
  theme_cowplot()+
  ylab("% Anthracosis") +
  xlab("Smoking Pack Years")

# SMOKING HISTORY PLOT
ever_never_smoker_plot | pack_year_plot
