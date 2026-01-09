library(data.table)
library(tidyr)
library(ggplot2)
library(gplots)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggpubr)

# LUCAS: For these 2 input files, can we just include data for the 15 cases included.
bal_data = fread("data/nanostringdata_matchedtoscans.csv")
bal_cohort_info = readxl::read_xlsx("data/clinicalforBALmatchedslide_deID.xlsx") %>%
  janitor::clean_names() %>%
  mutate(smoking=ifelse(smoking%in%c("minimal","former"),"EverSmoked",smoking))%>%
  mutate(smoking=ifelse(smoking%in%c("never"),"NeverSmoked",smoking)) %>%
  mutate(lobe_for_both_slide_and_bal=ifelse(lobe_for_both_slide_and_bal%in%c("RUL","LUL"),"UpperLobe",lobe_for_both_slide_and_bal)) %>%
  mutate(lobe_for_both_slide_and_bal=ifelse(lobe_for_both_slide_and_bal%in%c("RLL","LLL"),"LowerLobe",lobe_for_both_slide_and_bal))

bal_anthracosis_df = fread("data/BAL_anthracosis_measurements.csv")%>%
  janitor::clean_names() %>%
  filter(!image%in%c("BAL_21_N.ndpi","BAL_25_N.ndpi"))

bal_anthracosis_calls_df = bal_anthracosis_df %>%
  filter(classification=="Anthracosis") %>%
  mutate(id_for_scan=gsub("_N.ndpi","",image)) %>%
  mutate(percent_anthracosis=percent) %>%
  select(id_for_scan,percent_anthracosis)

bal_df = as.matrix(bal_data[,20:37])
rownames(bal_df) = bal_data$ProbeName

bal_cohort_info = bal_cohort_info %>%
  left_join(.,bal_anthracosis_calls_df,by="id_for_scan") %>%
  filter(id_for_scan%in%colnames(bal_df))

median_anthracosis_percent = bal_cohort_info %>%
  filter(!is.na(percent_anthracosis)) %>%
  pull(percent_anthracosis) %>%
  median()

bal_cohort_anthra = bal_cohort_info %>%
  filter(!is.na(percent_anthracosis)) %>%
  mutate(anthracosis_group=ifelse(percent_anthracosis>=median_anthracosis_percent,"High","0Low"))
bal_df = bal_df[,bal_cohort_anthra$id_for_scan]



pheno = bal_cohort_anthra %>%
  mutate(condition=anthracosis_group) %>%
  select(id_for_scan,condition) %>%
  data.frame()
rownames(pheno) = bal_cohort_anthra$id_for_scan


# Accounting for control and housekeeping gene expression
counts = bal_df
neg_ctrls <- grep("^NEG", rownames(counts), value = TRUE)
pos_ctrls <- grep("^POS", rownames(counts), value = TRUE)
hk_genes <- setdiff(
  intersect(
    rownames(counts),
    c("ACTB","GAPDH","RPLP0","RPL19","RPL10",
      "HPRT1","GUSB","B2M","TBP","EEF1G","PGK1","SDHA")
  ),
  c(neg_ctrls, pos_ctrls)
)

# Background Correction
neg_means <- colMeans(counts[neg_ctrls, , drop = FALSE])
neg_sds   <- apply(counts[neg_ctrls, , drop = FALSE], 2, sd)

bg_threshold <- neg_means + 2 * neg_sds

# Subtract background, floor at 1
bg_corrected <- sweep(counts, 2, bg_threshold, FUN = "-")
bg_corrected[bg_corrected < 1] <- 1


# Housekeeping normalization
hk_geo_means <- apply(bg_corrected[hk_genes, , drop = FALSE], 2,
                      function(x) exp(mean(log(x))))

hk_geo_means <- hk_geo_means / exp(mean(log(hk_geo_means)))

norm_counts <- sweep(bg_corrected, 2, hk_geo_means, FUN = "/")


# Log Normalization
log_counts <- log2(norm_counts)


# Design Matrix
pheno$condition <- factor(pheno$condition)

design <- model.matrix(~ condition, data = pheno)

#  Limma Diff Exp
library(limma)

fit <- lmFit(log_counts, design)
fit <- eBayes(fit)

res <- topTable(
  fit,
  coef = "conditionHigh",  # change to your contrast
  number = Inf,
  adjust.method = "BH"
)

limma_res_df = data.frame(res)
colnames(limma_res_df) = paste("limma",colnames(limma_res_df),sep="_")
limma_res_df$gene_name = rownames(res)


#  Correlation Test
anthracosis <- bal_cohort_anthra$percent_anthracosis
names(anthracosis) = bal_cohort_anthra$id_for_scan
anthracosis = anthracosis[ colnames(norm_counts)]

cor_results <- apply(norm_counts, 1, function(gene_exp) {
  ct <- cor.test(gene_exp, anthracosis, method = "spearman")
  c(cor = ct$estimate, p = ct$p.value)
})

cor_df <- as.data.frame(t(cor_results))
cor_df$FDR <- p.adjust(cor_df$p, method = "BH")
cor_res_df = data.frame(cor_df)
colnames(cor_res_df) = paste("spearman",colnames(cor_res_df),sep="_")
cor_res_df$gene_name = rownames(cor_df)

#  Expression Plotting
sig_genes <- rownames(res[
  signif(res$adj.P.Val,1) <= 0.1 & abs(res$logFC) >= 2,
])

all_results_df = full_join(limma_res_df,cor_res_df,by="gene_name")

all_results_sig = all_results_df %>%
  # filter(signif(limma_adj.P.Val,1)<=0.1&abs(limma_logFC)>=2.5&spearman_p<0.05) %>%
  filter(spearman_FDR<0.05)

norm_counts_df = norm_counts %>%
  reshape::melt() %>%
  rename(gene_name = X1) %>%
  rename(sample_name = X2)


norm_counts_df %>%
  left_join(.,bal_cohort_anthra,by=c("sample_name"="id_for_scan")) %>%
  filter(gene_name%in%all_results_sig$gene_name) %>%
  ggplot(aes(y=percent_anthracosis,x=value)) +
  geom_point() +
  facet_wrap(~gene_name,scales="free") +
  stat_cor(method = "spearman", 
           label.x.npc = "left", 
           label.y.npc = "top") +
  theme_cowplot() +
  geom_smooth(method="lm") +
  ylab("% Anthracosis") +
  xlab("Normalized gene expression counts")

#  GLM to assess the relative contribution of other variables to the gene expression
gene_names = unique(norm_counts_df$gene_name)
glm_res_df = lapply(gene_names,function(x){
  tt_df = norm_counts_df %>%
    filter(gene_name==x) %>%
    left_join(.,bal_cohort_anthra,by=c("sample_name"="id_for_scan"))
  
  g_res_df = glm(value~percent_anthracosis+age+sex+smoking+lobe_for_both_slide_and_bal,data=tt_df) %>%
    broom::tidy()
  g_res_df$variable = x
  g_res_df
})
glm_res_df = bind_rows(glm_res_df)

glm_res_anthracosis_df = glm_res_df %>%
  filter(term=="percent_anthracosis")

