
#------------------------------------------------------------------------------
# library
#------------------------------------------------------------------------------
library(data.table)
library(tidyr)
library(ggplot2)
library(gplots)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggpubr)

#------------------------------------------------------------------------------
# Figure 6: BAL Cohort Analysis
#------------------------------------------------------------------------------

bal_data = readxl::read_xlsx("../Cohort_data/BAL_nanostringdata.xlsx")
bal_cohort_info = fread("../Cohort_data/BAL_Cohort.csv") %>%
  janitor::clean_names() %>%
  mutate(smoking=ifelse(smoking%in%c("minimal","former"),"EverSmoked",smoking))%>%
  mutate(smoking=ifelse(smoking%in%c("never"),"NeverSmoked",smoking)) %>%
  mutate(lobe_for_both_slide_and_bal=ifelse(lobe_for_both_slide_and_bal%in%c("RUL","LUL"),"UpperLobe",lobe_for_both_slide_and_bal)) %>%
  mutate(lobe_for_both_slide_and_bal=ifelse(lobe_for_both_slide_and_bal%in%c("RLL","LLL"),"LowerLobe",lobe_for_both_slide_and_bal))

bal_df = as.matrix(bal_data[,11:25])
rownames(bal_df) = bal_data$ProbeName

median_anthracosis_percent = bal_cohort_info %>%
  filter(!is.na(anth_percent)) %>%
  pull(anth_percent) %>%
  median()

bal_cohort_anthra = bal_cohort_info %>%
  filter(!is.na(anth_percent)) %>%
  mutate(anthracosis_group=ifelse(anth_percent>=median_anthracosis_percent,"High","Low"))
fbal_df = bal_df[,bal_cohort_anthra$study_id]

bal_cohort_anthra %>%
  ggplot(aes(x=smoking,y=anth_percent)) +
  ggbeeswarm::geom_quasirandom() +
  stat_compare_means() +
  theme_cowplot()

pheno = bal_cohort_anthra %>%
  mutate(condition=anthracosis_group) %>%
  select(condition) %>%
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
  coef = "conditionLow",  # change to your contrast
  number = Inf,
  adjust.method = "BH"
)

limma_res_df = data.frame(res)
colnames(limma_res_df) = paste("limma",colnames(limma_res_df),sep="_")
limma_res_df$gene_name = rownames(res)


#  Correlation Test
anthracosis <- bal_cohort_anthra$anth_percent
names(anthracosis) = bal_cohort_anthra$study_id
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

spearman_dataset_result = all_results_df %>%
  mutate(Significance=ifelse(spearman_p<0.05&spearman_FDR<0.1,"p-value<0.05; FDR<0.1","None")) %>%
  mutate(Significance=ifelse(spearman_FDR<0.05&spearman_FDR<0.05,"p-value<0.05; FDR<0.05",Significance)) %>%
  mutate(gene_name_lab=ifelse(spearman_FDR<0.05,gene_name,"")) %>%
  arrange(spearman_cor.rho) %>%
  mutate(gene_name = factor(gene_name, levels = gene_name)) %>%
  ggplot(aes(x = spearman_cor.rho, y = gene_name, color = Significance)) +
  geom_segment(aes(x = 0, xend = spearman_cor.rho, yend = gene_name)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label=gene_name_lab),    force = 5,                 # default is 1
                           point.padding = 0.5,       # distance from point
                           box.padding = 0.7,         # distance between labels
                           min.segment.length = 0,    # always draw segments
                           max.overlaps = Inf) +
  theme_cowplot() +
  labs(x = "Spearman ρ", y = "Gene Name") +
  scale_color_manual(values=c("grey","red3","salmon")) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  coord_flip()


all_results_sig = all_results_df %>%
  filter(spearman_FDR<0.05)

norm_counts_df = norm_counts %>%
  reshape::melt() %>%
  rename(gene_name = X1) %>%
  rename(sample_name = X2)


correlation_plots = norm_counts_df %>%
  left_join(.,bal_cohort_anthra,by=c("sample_name"="study_id")) %>%
  filter(gene_name%in%all_results_sig$gene_name) %>%
  ggplot(aes(y=anth_percent,x=value)) +
  geom_point() +
  facet_wrap(~gene_name,scales="free",nrow=1) +
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
    left_join(.,bal_cohort_anthra,by=c("sample_name"="study_id"))
  
  g_res_df = glm(value~anth_percent+age+sex+smoking+lobe_for_both_slide_and_bal,data=tt_df) %>%
    broom::tidy()
  g_res_df$variable = x
  g_res_df
})
glm_res_df = bind_rows(glm_res_df)

glm_res_anthracosis_df = glm_res_df %>%
  filter(term=="anth_percent")





#  FGSEA
# https://nanostring.com/ncounterfiles/ for annotations

nanostring_annotations= fread("../Cohort_data/LBL-10043-08_nCounter_PanCancer_Immune_Profiling_Panel_Gene_List_Annotations.csv")
colnames(nanostring_annotations) <- as.character(nanostring_annotations[2,])
nanostring_annotations=nanostring_annotations[-c(1,2),] %>%
janitor::clean_names()

pathways <- nanostring_annotations %>%
  filter(!is.na(annotation)) %>%
  separate_rows(.,annotation,sep=", ") %>%
  filter(annotation!="",annotation!="Cell Type specific") %>%
  group_by(annotation) %>%
  summarise(genes = list(unique(gene_name))) %>%
  tibble::deframe()

library(fgsea)
ranks = all_results_df$spearman_cor.rho
names(ranks) = all_results_df$gene_name
ranks <- sort(ranks, decreasing = TRUE)

fg <- fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 6,
  maxSize = 200
)

fg %>%
  ggplot(aes(x=NES,y=-log(padj))) +
  geom_point(aes(size=log2err)) +
  ggrepel::geom_label_repel(aes(label=pathway))

pathway_order = fg %>%
  arrange(NES) %>%
  pull(pathway)

top_pathways = c(head(pathway_order,n=5),tail(pathway_order,n=5))

pathway_enrichment_plot = fg %>%
  mutate(Significance=ifelse(pval<0.05&padj<0.1,"p-value<0.05; FDR<0.1","None")) %>%
  mutate(Significance=ifelse(pval<0.05&padj<0.05,"p-value<0.05; FDR<0.05",Significance)) %>%
  filter(pathway%in%top_pathways) %>%
  mutate(direction=factor(ifelse(NES>0,"Activated","Suppressed"),levels=c("Suppressed","Activated"))) %>%
  mutate(pathway=factor(pathway,levels=pathway_order)) %>%
  ggplot(aes(x = NES, y = pathway, fill = Significance, size = size)) +
  geom_point(pch=21, color='black') +
  scale_fill_manual(values=c("grey","salmon")) +
  facet_grid(~direction, scales="free_x", space="free_y") +
  labs(x = "fgsea Normalized enrichment score", y = "Gene Set") 
  
figure_6 <- ( ((plot_spacer() | spearman_dataset_result) + plot_layout(widths=c(2,2))) / 
    correlation_plots / 
    (plot_spacer() | pathway_enrichment_plot) ) +
  plot_layout(heights=c(2,2,2))
#ggsave("figure_6_bal_cohort_analysis.png",figure_6,width=20,height=12)

