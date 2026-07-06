library(tidyverse)
library(pheatmap)
library(factoextra)
library(FactoMineR)
library(edgeR)
library(limma)
library(glue)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)
library(emmeans)


options(scipen=0)
set.seed(1)

in_dir <- "~/EpiBORL/git/scripts/rmd_files"
out_dir <- "~/EpiBORL/paper/scripts"

# Description file 

desc <- as_tibble(read.table(glue('{in_dir}/epiborl_sp_102024.txt'),h=T,sep="\t",check.names=F)) 

############################
## RNAseq quality control ##
############################

# FastQC stats

lfiler <- list.files(in_dir,pattern="_multiqc_fastqc_fastqc_trimmed.txt",full.names=T)

fastqc_stats <- map_dfr(lfiler,function(file){
  tab <- read.table(file,h=T,sep="\t",check.names=F) %>% 
    dplyr::select(Sample,`%GC`)
}) %>% 
  dplyr::rename("Barcode_RNA"="Sample")

# Aligmnent stats

star_stats <- read.table(glue('{in_dir}/all_multiqc_start.txt'),h=T) %>% 
  dplyr::select(Sample,total_reads,uniquely_mapped_percent,multimapped_percent) %>% 
  dplyr::rename("Barcode_RNA"="Sample") %>% 
  left_join(statsr)

# Coding gene stats 

counts_stats <- read.table(glue('{in_dir}/all_counts_genetype.txt'),h=T,sep=",",row.names=1) %>% 
  rownames_to_column("Barcode_RNA") %>% 
  dplyr::select(Barcode_RNA,protein_coding)

stats <- fastqc_stats %>% 
  left_join(star_stats) %>% 
  left_join(counts_stats) %>% 
  filter(!is.na(Barcode_RNA)) %>% 
  dplyr::select(Barcode_RNA,total_reads:protein_coding,`%GC`)

# thresholds
min_gc <- 30
max_gc <- 45
min_unique <- 65
min_genes <- 10000

stats <- stats %>% 
  mutate("Category" = case_when(`%GC` < min_gc ~ "Low %GC",
                                `%GC`> max_gc ~ "High %GC", 
                                uniquely_mapped_percent< min_unique ~ "Low %Unique", 
                                protein_coding < min_genes ~ "Low Nb Genes", 
                                TRUE ~ "OK")) 


desc_filter <- desc %>% filter(Barcode_RNA %in% (stats %>% filter(Category=="OK") %>% pull(Barcode_RNA)))


desc_filter <- desc_filter %>% mutate(disease_type=gsub("NK",NA,disease_type),
                                      `Best response`=gsub("NK",NA,`Best response`),
                                      `HPV status`=gsub("NK",NA,`HPV status`)) %>% 
                   mutate("Ages"=if_else(Age<70,"<70",">=70"),
                          "Alcohols"=if_else(Alcohol=="NON CONSUMER" ,"non consumer","consumer"),
                          "Tobaccos"=if_else(Tobacco=="NON CONSUMER","non consumer","consumer"),
                          "CPS"=case_when(`CPS %` < 1 ~ "[0-0.9]", 
                                          `CPS %` >= 1 & `CPS %` <=19 ~ "[1-19]", 
                                          TRUE ~ ">=20")) %>%
                   mutate(`Disease type` = case_when(grepl("Loc",disease_type) ~ "Local",
                                                     grepl("Meta",disease_type) ~ "Metastatic", 
                                                     disease_type == "Both" ~ "Both", 
                                                     TRUE ~ NA )) %>% 
                  mutate(`Best response`=factor(`Best response`,levels=c("CR","PR","SD","PD")), 
                                      Response=case_when(`Best response` %in% c("CR","PR") ~ "R",
                                                         `Best response` %in% c("SD","PD") ~ "NR",
                                                         TRUE ~ "NK"), 
                                      Response=factor(Response,levels=c("R","NR","NK"))) 

# 176 RNAseq samples kept after QC 
  
################
## HPV RNAseq ##
################

hpvkite_dir <- "~/EpiBORL/analyses/hpv-kite"
tversky_thres <- 0.035

lfile <- list.files(hpvkite_dir,"*txt",recursive=T,full.names=T)

lfile_rna <- lfile[grep("/D",lfile)]

res_hpv <- map_dfr(lfile_rna,function(file){
  fastq <- sapply(strsplit(basename(file),".",fixed=T),"[",2)
  name <- sapply(strsplit(basename(file),".",fixed=T),"[",1)
  read.table(file) %>% mutate(ID=name,fastq=fastq)
  }) %>% 
  as_tibble() 

colnames(res_hpv) <- c("HPV1","score1","HPV2","score2","HPV3","score3","read","Barcode","fastq")

res_hpv <- res_hpv %>% dplyr::select(Barcode,fastq,HPV1,score1) %>% arrange(HPV1,-score1) %>%
  mutate("HPV_RNAseq_3.5%"=if_else(score1>=tversky_thres & HPV1=="HPV16REF" ,"P16pos","P16neg"))

write.table(res_hpv,glue('{out_dir}/epiborl_hpvkite.txt'),quote=F,col.names=T,row.names=F,sep="\t")

desc_filter <- desc_filter %>% left_join(res_hpv[,c("Barcode","HPV_RNAseq_3.5%")],by=c("Barcode_RNA"="Barcode")) %>%
  mutate(HPV16=factor(`HPV_RNAseq_3.5%`,levels=c("P16neg","P16pos"))) %>% 
  arrange(`Best response`,HPV16)  %>% 
  mutate(IDAT=glue('{Array}_{Slide}'), IDAT = if_else(IDAT=="NA_NA",NA,IDAT))

desc_filter <- desc_filter %>% dplyr::select(ID,Barcode_RNA, IDAT, Sex, Ages, Alcohols, Tobaccos, `Disease type`,`Primary site of cancer`, HPV16, `Best response`,Response) %>% 
  mutate(`Primary site of cancer`=factor(`Primary site of cancer`,levels=c("Hypopharynx","Larynx","Oral cavity","Oropharynx","Other")))

write.table(desc_filter,glue('{out_dir}/files/epiborl_desc_filter.txt'),quote=F,col.names=T,row.names=F,sep="\t")

#########################
## Expression analysis ##
#########################


# Color code 

disease_color_code <- c("Both"="#E89989","Local"="#E34122","Metastatic"="#9C311C","NA"="#c4c4bf")
br_color_code <- c("CR"="#CA7B3A","PR"="#BE8C64","SD"="#157DC6","PD"="#427231")
hpv_color_code <- c("P16neg"="#315a72","P16pos"="#f3b00c")
alc_color_code <- c("consumer"="#b6ad95","non consumer"="#9e7f32")
tob_color_code <- c("consumer"="#b6ad95","non consumer"="#9e7f32")
prim_color_code <- c("Oral cavity"="#88AD60","Oropharynx"="#54950F","Hypopharynx"="#4D7226","Larynx"="#3B6B09","Other"="#2E5307")
deg_color_code <- c("NS"="gray85","Downregulated"="steelblue","Upregulated"="tomato3")
resp_color_code <- c("R"="#D99E0F", "NR"="#728A69","NK"="#c4c4bf")
pathway_color <- c("Altered"="steelblue")

annocol <- list("Disease_type"=disease_color_code,
                "Best_response"=br_color_code,
                "HPV16"=hpv_color_code,
                "Alcohol"=alc_color_code,
                "Tobacco"=tob_color_code,
                "Primary_tumor_site"=prim_color_code,
                "Response"=resp_color_code,
                "DEG"=deg_color_code)


coldata <- desc_filter %>% remove_rownames() %>% 
  column_to_rownames("Barcode_RNA") %>% 
  rename("Primary_tumor_site"="Primary site of cancer","Best_response"="Best response","Disease_type"="Disease type", "Tobacco"="Tobaccos","Alcohol"="Alcohols") %>% 
  dplyr::select("Response","HPV16","Disease_type","Primary_tumor_site","Best_response", "Tobacco","Alcohol")


# Normalized counts

counts <- read.csv(glue('{in_dir}/all_tablecounts_raw.csv'),row.names=1,check.names=F) %>% 
  rownames_to_column("ENSG") %>% 
  column_to_rownames("ENSG") %>%
  dplyr::select(desc2$Barcode_RNA)

genes_anno <- read.csv(glue('{in_dir}/hg38_tableannot.csv'),row.names=1,stringsAsFactors = FALSE) %>% 
  add_row(gene_id="ENSG00000286219","gene_name"="NOTCH2NLC")

#focus on protein coding genes
allcor_id <- read.table(glue('{in_dir}/gencode_v34.bed'),sep="\t",stringsAsFactors=FALSE)

counts <- counts[intersect(rownames(counts),allcor_id[,1]),]

dge <- DGEList(counts=counts,genes=rownames(counts),group=coldata$Run)
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

v1 <- voom(dge, plot = FALSE)
v1t <- v1$E

rownames(v1t) <- plyr::mapvalues(rownames(v1t),genes_anno$gene_id,genes_anno$gene_name,warn_missing = FALSE)

write.table(v1t,glue('{out_dir}/files/epiborl_rna_normcounts_hugo.txt'),quote=F,col.names=T,row.names=T,sep="\t") #14158 genes

###########################
## Differential analysis ##
########################### 

padj_th <- 0.05
fc_th <- 1

ngb <- 100

nc <- 15
h <- 1000
kegg_species <- "hsa"


# CR/PR vs SD/PD 

coldata_br <- coldata %>% filter(!is.na(Best_response))
counts_br <- counts %>% dplyr::select(coldata_br %>% rownames_to_column("Barcode") %>% pull(Barcode))

dds <- DGEList(counts=counts_br,genes=rownames(counts_br), group=coldata_br$`Best_response`)
design <- model.matrix(~ -1 + `Best_response`, data = coldata_br)
colnames(design) <- make.names(levels(as.factor(coldata_br$`Best_response`)))

keep <- filterByExpr(dds, design)
dds <- dds[keep,,keep.lib.sizes=FALSE]
dds <- calcNormFactors(dds)

vbr <- voom(dds, design = design, plot = FALSE)
vbrt <- vbr$E

contrast.matrix <- makeContrasts((CR+PR)/2-(SD+PD)/2, levels=design)
res_fit_br <- lmFit(vbr,design=design)
fit_br <- contrasts.fit(res_fit_br, contrast.matrix)
fit_br <- eBayes(fit_br)

tt_br <- topTable(fit_br,number=Inf) %>% 
  left_join(genes_anno,by=c("genes"="gene_id")) %>% 
  mutate(log10padj=-log10(`adj.P.Val`)) %>% 
  mutate(DEG="NS")

tt_br <- tt_br  %>%
  mutate(DEG = case_when((`adj.P.Val`<=padj_th & logFC >= fc_th) ~ 'Upregulated',
  (`adj.P.Val`<=padj_th & logFC <= -fc_th ) ~ 'Downregulated',
  TRUE ~ DEG))
 
tt_brs <- tt_br  %>% filter(`adj.P.Val`<=padj_th & abs(logFC)>=fc_th )  %>% arrange(logFC) 

volcano_br <- ggplot(tt_br,aes(x=logFC,y=log10padj,color=DEG)) +
  geom_point(size=1,aes(text=paste0("Gene:",gene_name))) +
  theme_bw()+labs(x="log2(FoldChange)", y ="-log10(pvalue)") + 
  geom_hline(yintercept = -log10(padj_th),lty=2,col="steelblue") + 
  geom_vline(xintercept = c(-fc_th,fc_th),lty=2,col="red") + 
  scale_color_manual(values=annocol$DEG) + 
  ggtitle("Best response: CR/PR vs SD/PD") + 
  geom_label_repel(data=tt_brs,label=tt_brs$gene_name,show.legend = FALSE)

ggsave(glue('{out_dir}/files/volcano_crpr_sdpd.pdf'),volcano_brb,width=8,height=8)

if(nrow(tt_brs)>ngb){val <- c(1:ngb,(nrow(tt_brs)-ngb):nrow(tt_brs))}else{val <- c(1:nrow(tt_brs))}

pbr <- vbrt[which(rownames(vbrt)%in%tt_brs[val,"genes"]),]
rownames(pbr) <- plyr::mapvalues(rownames(pbr),genes_anno$gene_id,genes_anno$gene_name,warn_missing = FALSE)

rbr <- data.frame("DEG"=tt_brs[val,"DEG"])
rownames(rbr) <- tt_brs[val,"gene_name"]

ph1 <- pheatmap(pbr[rownames(rbr),], 
                scale="row",
                show_rownames=TRUE,
                show_colname=FALSE, 
                annotation_col=coldata, 
                clustering_method="ward.D2", 
                annotation_colors=annocol,
                annotation_row=rbr, 
                clustering_distance_rows="correlation",
                clustering_distance_cols="correlation",
                cluster_rows=TRUE,
                cluster_cols = TRUE,
                main=glue::glue('Top DEG genes for Best response: CR/PR vs SD/PD, abs(logFC)>={fc_th}, adj pval <= {padj_th}'),
                fontsize = 7, 
                fontsize_col=4)

ggsave(glue('{out_dir}/files/heatmap_deg_crpr_sdpd.pdf'),ph1,width=12,height=8)


# CR/PR/SD vs PD 


contrast.matrix <- makeContrasts(((CR+PR+SD)/3)-PD, levels=design)
res_fit_br2 <- lmFit(vbr,design=design)
fit_br2 <- contrasts.fit(res_fit_br2, contrast.matrix)
fit_br2 <- eBayes(fit_br2)

tt_br2 <- topTable(fit_br2,number=Inf) %>% 
  left_join(genes_anno,by=c("genes"="gene_id")) %>% 
  mutate(log10padj=-log10(`adj.P.Val`))%>% 
  mutate(DEG="NS")

tt_br2 <- tt_br2  %>%
  mutate(DEG = case_when((`adj.P.Val`<=padj_th & logFC >= fc_th) ~ 'Upregulated',
  (`adj.P.Val`<=padj_th & logFC <= -fc_th ) ~ 'Downregulated', TRUE ~ DEG))

tt_br2s <- tt_br2  %>% filter(`adj.P.Val`<=padj_th & abs(logFC)>=fc_th )  %>% arrange(logFC) 


volcano_br2 <- ggplot(tt_br2,aes(x=logFC,y=log10padj,color=DEG))+
  geom_point(size=1,aes(text=paste0("Gene:",gene_name)))+
  theme_bw()+labs(x="log2(FoldChange)", y ="-log10(pvalue)") + 
  geom_hline(yintercept = -log10(padj_th),lty=2,col="steelblue") + 
  geom_vline(xintercept = c(-fc_th,fc_th),lty=2,col="red") + 
  scale_color_manual(values=annocol$DEG) + 
  ggtitle("Best response: CR/PR/SD vs PD") + 
  geom_label_repel(data=tt_br2s,label=tt_br2s$gene_name, show.legend = FALSE)

ggsave(glue('{out_dir}/files/volcano_crprsd_pd.pdf'),volcano_br2,width=8,height=8)


if(nrow(tt_br2s)>ngb){val <- c(1:ngb,(nrow(tt_br2s)-ngb):nrow(tt_br2s))}else{val <- c(1:nrow(tt_br2s))}

pbr2 <- vbrt[which(rownames(vbrt)%in%tt_br2s[val,"genes"]),]
rownames(pbr2) <- plyr::mapvalues(rownames(pbr2),genes_anno$gene_id,genes_anno$gene_name,warn_missing = FALSE)

rbr2 <- data.frame("DEG"=tt_br2s[val,"DEG"])
rownames(rbr2) <- tt_br2s[val,"gene_name"]


ph2 <- pheatmap(pbr2[rownames(rbr2),],
                scale="row", 
                show_rownames=TRUE,
                show_colname=FALSE, 
                annotation_col=coldata, 
                clustering_method="ward.D2", 
                annotation_colors=annocol,
                annotation_row=rbr2, 
                clustering_distance_rows="correlation",
                clustering_distance_cols="correlation",
                cluster_rows=TRUE,
                cluster_cols = TRUE,
                main=glue::glue('Top DEG genes for Best response: CR/PR/SD vs PD, abs(logFC)>={fc_th}, adj pval <= {padj_th}'),
                fontsize = 7, fontsize_col=4, fontsize_row=8)

ggsave(glue('{out_dir}/files/heatmap_deg_crprsd_pd.pdf'),ph2,width=12,height=5)


######################
## Pathway analysis ##
######################

deg_func <- function(df,outDir,padj_th,fc_th){
  allg <- bitr(df$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  df <- df %>% left_join(allg,by=c("gene_name"="SYMBOL"))
  
  All_keggk <- clusterProfiler::enrichKEGG(gene  = df$ENTREZID, organism="hsa",keyType = "kegg")
  
  if(nrow(as.data.frame(All_keggk))!=0){
    
    p <- plot(dotplot(All_keggk,showCategory=nc,orderBy="GeneRatio",color="p.adjust",decreasing=T,font.size=15, title = glue::glue('KEGG'))) +
      scale_fill_gradient(low = "red", high = "blue")
    ggsave(glue::glue('{outDir}_pathway_KEGG_all_dotplot.png'),p,width=12,height=10,dpi=300)
    
  }
}

# CR/PR vs SD/PD

deg_func(tt_brs,glue('{out_dir}/files/crpr_sdpd'),padj_th,fc_th)

# CR/PR/SD vs PD

deg_func(tt_br2s,glue('{out_dir}/files/crprsd_pd'),padj_th,fc_th)


##########################
## Immune deconvolution ##
##########################


r <- read.table(glue('{in_dir}/immunedeconv/D1490-D1578-D1597_quantiseq_results.txt'),h=T,sep="\t",check.names=FALSE)

r2 <- r %>% tidyr::gather(sample, fraction, -cell_type) %>% 
                      dplyr::group_by(sample) %>% 
                      dplyr::mutate(prop = fraction / (1-fraction[cell_type == "uncharacterized cell"])) %>% 
                      dplyr::ungroup() 

r2$sample <- factor(r2$sample, levels=r %>% tidyr::gather(sample, fraction, -cell_type) %>% 
                       dplyr::filter(cell_type=="uncharacterized cell") %>% arrange(fraction) %>% pull(sample))

r2 <- r2 %>% left_join(coldata %>% rownames_to_column("sample") %>% dplyr::rename("Best response"="Best_response") %>% dplyr::select(sample,`Best response`))%>% arrange(prop)

color_code <- colorRampPalette(brewer.pal(8, "Paired"))(length(unique(r2$cell_type))-1)
names(color_code) <- grep("uncharacterized",unique(r2$cell_type),value=T,invert=T)
color_code <- c(color_code,"uncharacterized cell"="gray90")

# CR/PR vs SD/PD

rq <- as.data.frame(t(r))
colnames(rq) <- rq[1,]
rq <- rq[-1,]

rq <- rq %>% rownames_to_column("Barcode_RNA") %>% 
  left_join(coldata %>% rownames_to_column("Barcode_RNA") %>% dplyr::rename("Best response"="Best_response") %>% dplyr::select(Barcode_RNA,`Best response`)) %>% 
  rename("Barcode"="Barcode_RNA") %>%  
  mutate(Response=if_else(`Best response` %in%c("CR","PR"),"R","NR")) %>% 
  filter(!is.na(`Best response`))

rq$Response <- factor(rq$Response,levels=c("R","NR"))
rq$`Best response` <- factor(rq$`Best response`,levels=c("CR","PR","SD","PD"))
rq[,2:(ncol(rq)-2)] <- apply(rq[,2:(ncol(rq)-2)],2,function(x){as.numeric(x)})


aa <- rq %>% as_tibble() %>% 
  dplyr::select(-contains("counts")) %>%
  pivot_longer(cols = -c(Barcode, Response,`Best response`), names_to = "cell_type", values_to = "prop") %>%
  mutate(prop = 100 * as.numeric(prop), prop_other = 100 - prop) %>%
  group_by(cell_type) %>%
  mutate(prop = cbind(prop, prop_other)) %>%
  nest() %>% 
  filter(cell_type!="uncharacterized cell") %>%
  rowwise() %>%
  summarise(mod = list(glm(prop ~ Response ,family = quasibinomial(), data = data)))

aa_res <- aa %>%
  rowwise() %>%
  mutate(coeffs = list(broom::tidy(mod)),
                       em = list(emmeans(mod,  ~Response ,type="response")),
                       p = list(plot(em) + theme_light() + labs(title = cell_type))
                                )

get_pval <- function(mod){
  broom::tidy(mod) %>% 
    filter(term == "ResponseNR") %>% 
    pull(p.value)
}

get_est <- function(mod){
  broom::tidy(mod) %>% 
    filter(term == "ResponseNR") %>% 
    pull(estimate)
}

aa_res <- aa_res %>% 
  rowwise() %>% 
  mutate(pval = get_pval(mod), 
         est = get_est(mod),  
         p = list(p + theme(plot.title = element_text(face = "bold")) + labs(title = str_glue("{cell_type} {if_else(pval<=0.05, '*', '')}"), subtitle = if_else(est >= 0, "NR > R", "NR < R" ))+ xlab("Estimated marginal mean"))) 


p_all <- patchwork::wrap_plots(arrange(aa_res, est)$p) 

p_all

ggsave(glue('{out_dir}/files/quantiseq_crpr_sdpd_plot_woUncharacterized_cell.pdf'),p_all,width=16,height=12)

aa_dt <- aa_res %>%
  rowwise() %>%
  mutate(em = list(broom::tidy(em))) %>%
  unnest(em) %>% 
  mutate(prob=prob*100,cell_type=factor(cell_type,levels=arrange(aa_res, est)$cell_type)) %>% 
  group_by(cell_type) %>% 
  mutate(ratio=case_when(prob[Response=="R"]>prob[Response=="NR"] ~ round(prob[Response=="R"]/prob[Response=="NR"],2),
                         TRUE ~ round(prob[Response=="NR"]/prob[Response=="R"],2)),
         Group = case_when(prob[Response=="R"]>prob[Response=="NR"] ~ "R/NR", TRUE ~ "NR/R"),
         Proba = case_when(prob[Response=="R"]>prob[Response=="NR"] ~ glue::glue('{round(prob[Response=="R"],2)}/{round(prob[Response=="NR"],2)}'), TRUE ~ glue::glue('{round(prob[Response=="NR"],2)}/{round(prob[Response=="R"],2)}')),
         Signif= if_else(pval<=0.05, '*', '')) %>% 
  filter(Response=="R") %>% 
  dplyr::select(cell_type,Proba, Group, ratio,pval,Signif) %>% 
  arrange(cell_type)

write.table(aa_dt,glue('{out_dir}/files/quantiseq_crpr_sdpd_table_prop.txt'),quote=F,col.names=T,row.names=F,sep="\t")    

rq_all2 <- rq %>% 
  dplyr::select(-contains("counts")) %>% 
  pivot_longer(cols = -c(Barcode, Response, `Best response`), names_to = "cell_type", values_to = "prop") %>% 
  mutate(prop = 100 * as.numeric(prop), cell_type=factor(cell_type,levels=arrange(aa_res, est)$cell_type)) %>% filter(cell_type!="uncharacterized cell") %>% ggplot(aes(x = prop, y =  `Best response`, fill = Response)) + 
  geom_boxplot() + 
  facet_wrap(~cell_type, scales = "free_x") + 
  theme_light() + 
  scale_fill_brewer(palette= "Set2") + 
  theme(strip.text = element_text(color = "gray10",face="bold"), strip.background = element_rect(fill = "grey95")) + 
  labs(title="Proportion of immune cells") + xlab("Proportion") 

ggsave(glue('{out_dir}/files/quantiseq_crpr_sdpd_table_prop_proportion_woUncharacterized_max50.pdf'),rq_all2 + facet_wrap(~cell_type,scales="free_y")+coord_cartesian(xlim = c(NA, 50)),width=16,height=12)


# CR/PR/SD vs PD

rqb <- as.data.frame(t(r))
colnames(rqb) <- rqb[1,]
rqb <- rqb[-1,]

rqb <- rqb %>% rownames_to_column("Barcode_RNA") %>% 
  left_join(coldata %>% rownames_to_column("Barcode_RNA") %>% dplyr::rename("Best response"="Best_response") %>% dplyr::select(Barcode_RNA,`Best response`)) %>% 
  rename("Barcode"="Barcode_RNA") %>%  
  mutate(Response=if_else(`Best response` %in%c("CR","PR","SD"),"R","NR"), 
         Response=factor(Response,levels=c("R","NR"))) %>% 
  filter(!is.na(`Best response`))

aab <- rqb %>% as_tibble() %>% 
  dplyr::select(-contains("counts")) %>%
  pivot_longer(cols = -c(Barcode, Response,`Best response`), names_to = "cell_type", values_to = "prop") %>%
  mutate(prop = 100 * as.numeric(prop), prop_other = 100 - prop) %>%
  group_by(cell_type) %>%
  mutate(prop = cbind(prop, prop_other)) %>%
  nest() %>%
  filter(cell_type!="uncharacterized cell") %>%
  rowwise() %>%
  summarise(mod = list(glm(prop ~ Response ,family = quasibinomial(), data = data)))

aa_resb <- aab %>%
  rowwise() %>%
  mutate(coeffs = list(broom::tidy(mod)),
                       em = list(emmeans(mod,  ~Response ,type="response")),
                       p = list(plot(em) + theme_light() + labs(title = cell_type))
                                )

get_pval <- function(mod){
  broom::tidy(mod) %>% 
    filter(term == "ResponseNR") %>% 
    pull(p.value)
}

get_est <- function(mod){
  broom::tidy(mod) %>% 
    filter(term == "ResponseNR") %>% 
    pull(estimate)
}

aa_resb <- aa_resb %>% 
  rowwise() %>% 
  mutate(pval = get_pval(mod), 
         est = get_est(mod),  
         p = list(p + theme(plot.title = element_text(face = "bold")) + labs(title = str_glue("{cell_type} {if_else(pval<=0.05, '*', '')}"), subtitle = if_else(est >= 0, "NR > R", "NR < R" ))+ xlab("Estimated marginal mean"))) 


p_allb <- patchwork::wrap_plots(arrange(aa_resb, est)$p) 

ggsave(glue('{out_dir}/files/quantiseq_crprsd_pd_plot.pdf'),p_allb,width=16,height=12)

aa_dtb <- aa_resb %>%
  rowwise() %>%
  mutate(em = list(broom::tidy(em))) %>%
  unnest(em) %>% 
  mutate(prob=prob*100,cell_type=factor(cell_type,levels=arrange(aa_res, est)$cell_type)) %>% 
  group_by(cell_type) %>% 
  mutate(ratio=case_when(prob[Response=="R"]>prob[Response=="NR"] ~ round(prob[Response=="R"]/prob[Response=="NR"],2),
                         TRUE ~ round(prob[Response=="NR"]/prob[Response=="R"],2)),
         Group = case_when(prob[Response=="R"]>prob[Response=="NR"] ~ "R/NR", TRUE ~ "NR/R"),
         Proba = case_when(prob[Response=="R"]>prob[Response=="NR"] ~ glue::glue('{round(prob[Response=="R"],2)}/{round(prob[Response=="NR"],2)}'), TRUE ~ glue::glue('{round(prob[Response=="NR"],2)}/{round(prob[Response=="R"],2)}')),
         Signif= if_else(pval<=0.05, '*', '')) %>% 
  filter(Response=="R") %>% 
  dplyr::select(cell_type,Proba, Group, ratio,pval,Signif) %>% 
  arrange(cell_type)

write.table(aa_dtb,glue('{out_dir}/files/quantiseq_crprsd_pd_table_prop.txt'),quote=F,col.names=T,row.names=F,sep="\t")    

rqb_all2 <- rqb %>% 
  dplyr::select(-contains("counts")) %>% 
  pivot_longer(cols = -c(Barcode, Response, `Best response`), names_to = "cell_type", values_to = "prop") %>% 
  mutate(prop = 100 * as.numeric(prop), cell_type=factor(cell_type,levels=arrange(aa_res, est)$cell_type)) %>% filter(cell_type!="uncharacterized cell") %>% ggplot(aes(x = prop, y =  `Best response`, fill = Response)) + 
  geom_boxplot() + 
  facet_wrap(~cell_type, scales = "free_x") + 
  theme_light() + 
  scale_fill_brewer(palette= "Set2") + 
  theme(strip.text = element_text(color = "gray10",face="bold"), strip.background = element_rect(fill = "grey95")) + 
  labs(title="Proportion of immune cells") + xlab("Proportion") 


ggsave(glue('{out_dir}/files/quantiseq_crprsd_pd_table_prop_proportion_woUncharacterized_max50.pdf'),rqb_all2 + facet_wrap(~cell_type,scales="free_y")+coord_cartesian(xlim = c(NA, 50)),width=16,height=12)
