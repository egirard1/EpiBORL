Sys.setenv(EXPERIMENT_HUB_CACHE = "~/.ExperimentHub_newR") # for sesame cache 

library(tidyverse)
library(pheatmap)
library(glue)
library(ggrepel)
library(RColorBrewer)
library(limma)
library(doParallel)

library(sesame) 
library(SummarizedExperiment)
library(ChAMP)  
library(ChAMPdata) # dev version installed from github to get AnnoEPICv2 annotation
library(DMRcate)
data(AnnoEPICv2) 
library(BiocParallel)
options(scipen=0)
set.seed(1)

# Need to install preprocessCore with a disable threading
# BiocManager::install( "preprocessCore", configure.args = c(preprocessCore = "--disable-threading"), force = TRUE, update = TRUE, type = "source")

sesameDataCache()

num.cores <- 8
input_dir <- "~/EpiBORL/git/EPIC_data/idat/"
output_dir <- "~/EpiBORL/paper/scripts/files/"

# getting beta values using sesame and their processing QCDPB

addr <- read.table("~/EpiBORL/git/EPIC_data/sesame_EPICV2_hg38_manifest.txt",sep="\t",h=T)
prepv <- "QCDPB" #QualityMask, inferIniniumIChannel, dyeBiasNL, pOOBAH,noob

betas <- openSesame(input_dir, prep=prepv, manifest=addr, BPPARAM = BiocParallel::MulticoreParam(num.cores), collapseToPfx = FALSE) #937690 probes
saveRDS(betas,glue('{output_dir}/epiborl_betas.Rds'))

# impute missing values using CHAMP

betas_imput <- champ.impute(beta = as.matrix(betas), pd = NULL, SampleCutoff=0.75,ProbeCutoff=0.2) # no sample removal, # 683017 probes
saveRDS(betas_imput,glue('{output_dir}/epiborl_betas_impute.Rds'))


## Description file 

in_dir <- "~/EpiBORL/git/scripts/rmd_files"

desc <- as_tibble(read.table(glue('{in_dir}/epiborl_sp_102024.txt'),h=T,sep="\t",check.names=F)) 
desc_sp <- desc %>% dplyr::select(Bloc,Barcode_RNA,EPIC_sample,Array,Slide, `Best response`)  %>% 
  mutate(IDAT=paste(Array,Slide,sep="_")) %>% 
  filter(!is.na(Slide)) %>% 
  filter(`Best response`!="NK") %>% 
  mutate(Status = if_else(`Best response` %in% c("CR","PR"),"R","NR"),
         Status2 = if_else(`Best response` %in% c("CR","PR","SD"),"R","NR"),
         Status=factor(Status, levels=c("R","NR")),
         Status2=factor(Status2, levels=c("R","NR"))) %>%
  as.data.frame()

rownames(desc_sp) <- desc_sp$IDAT

## Annotation using EPICv2 manifest

gene <- read.table("/data/kdi_prod/.kdi/project_workspace_0/1975/acl/01.00/git/EPIC_data/EPICv2.hg38.manifest.gencode.v41.tsv",h=T,sep="\t") %>% 
  as_tibble() %>% 
  mutate(Probe_ID=sapply(strsplit(probeID,"_"),"[",1)) %>% 
  dplyr::select(CpG_chrm:genesUniq,probeID,Probe_ID)

onco <- read.table("~/GIT/snpeff_update/oncokb_genes.txt",h=T)

all_gene <- gene %>% 
  separate_rows(genesUniq,sep=";") %>% 
  left_join(onco,by=c("genesUniq"="Gene")) %>% 
  dplyr::rename("Gene"="genesUniq")

all_gene_onco <- all_gene %>% filter(!is.na(Oncokb)) 

### DMP between R (CR/PR) vs NR (SD/PD)

# params

minProbes <- 7
pval_thresh <- 0.01
est_thresh <- 0.25
logFC_thresh <- 0.25

pval_thresh_dmr <- 0.05

##### using CHAMP

source("~/EpiBORL/git/scripts/champ.DMP.R") # adapt for EPICv2

# CR/PR vs SD/PD 

myDMPchamp <- champ.DMP(beta = betas_imput[,rownames(desc_sp)],
                        pheno = desc_sp[,"Status"], 
                        compare.group=c("R","NR"), 
                        arraytype = "EPICv2",  
                        adjPVal = 1)

saveRDS(myDMPchamp,glue('{output_dir}/epiborl_champ_dmp_crpr_sdpd.Rds'))

# head(myDMPchamp$NR_to_R)

dmpdr_champ <- myDMPchamp$NR_to_R %>% 
  mutate(DEG="NS") %>% 
  mutate(DEG=case_when(round(adj.P.Val,3) <= pval_thresh & logFC >= logFC_thresh ~ "R<NR", 
						  				round(adj.P.Val,3) <= pval_thresh & logFC <= -logFC_thresh ~  "R>NR", 
						  				TRUE ~ DEG))  %>% 
  mutate(DEG=factor(DEG,levels=c("R>NR","R<NR","NS")))

dmpdr_champ$probeID <- rownames(dmpdr_champ)

dmpdr_champ <- dmpdr_champ %>% group_by(DEG) %>% 
  arrange(adj.P.Val,-logFC,.by_group=TRUE) %>% 
  left_join(all_gene[,c("probeID","Probe_ID","Gene","Oncokb")])

write.table(dmpdr_champ,glue('{output_dir}/epiborl_champ_dmp_crpr_sdpd_adjpval0.01_logFC_0.25.txt'),quote=F,col.names=T,row.names=F,sep="\t")


# CR/PR/SD vs PD 

myDMPchamp2 <- champ.DMP(beta = betas_imput[,rownames(desc_sp)],
                        pheno = desc_sp[,"Status2"], 
                        compare.group=c("R","NR"), 
                        arraytype = "EPICv2",  
                        adjPVal = 1)

saveRDS(myDMPchamp2,glue('{output_dir}/epiborl_champ_dmp_crprsd_pd.Rds'))

# head(myDMPchamp2$NR_to_R)

dmpdr_champ2 <- myDMPchamp2$NR_to_R %>% 
  mutate(DEG="NS") %>% 
  mutate(DEG=case_when(round(adj.P.Val,3) <= pval_thresh & logFC >= logFC_thresh ~ "R<NR", 
						  				round(adj.P.Val,3) <= pval_thresh & logFC <= -logFC_thresh ~  "R>NR", 
						  				TRUE ~ DEG))  %>% 
  mutate(DEG=factor(DEG,levels=c("R>NR","R<NR","NS")))

dmpdr_champ2$probeID <- rownames(dmpdr_champ2)

dmpdr_champ2 <- dmpdr_champ2 %>% group_by(DEG) %>% 
  arrange(adj.P.Val,-logFC,.by_group=TRUE) %>% 
  left_join(all_gene[,c("probeID","Probe_ID","Gene","Oncokb")])

write.table(dmpdr_champ2,glue('{output_dir}/epiborl_champ_dmp_crprsd_pd_adjpval0.01_logFC_0.25.txt'),quote=F,col.names=T,row.names=F,sep="\t")

## using sesame

gc()

betas_imput <- betas_imput[,rownames(desc_sp)]

# CR/PR vs SD/PD

se_dup <- SummarizedExperiment(assays=list(betas=as.matrix(betas_imput)), 
                               colData=(desc_sp[colnames(betas_imput),] %>% 
                                          dplyr::select(IDAT,Status) %>% 
                                          mutate(Status=factor(Status,levels=c("R","NR")))))

smry_dup <- DML(se_dup, ~ Status, BPPARAM = BiocParallel::MulticoreParam(num.cores))

saveRDS(smry_dup, glue('{output_dir}/epiborl_sesame_dmp_crpr_sdpd.Rds'))

gc(verbose = FALSE)

dmp_sesame  <-  summaryExtractTest(smry_dup)
colnames(dmp_sesame) 
saveRDS(dmp_sesame, glue('{output_dir}/epiborl_sesame_dmp_crpr_sdpd_extract.Rds')) 

dmp_sesame <- dmp_sesame %>% 
  mutate(AdjPval_StatusNR=p.adjust(Pval_StatusNR,method="BH")) %>%
  mutate(DEGadj="NS") %>% 
  mutate(DEGadj=case_when(round(AdjPval_StatusNR,3) <= pval_thresh & Est_StatusNR >= est_thresh ~ "R<NR",
                          round(AdjPval_StatusNR,3) <= pval_thresh & Est_StatusNR <= -est_thresh ~ "R>NR",
                          TRUE ~ DEGadj)) %>%
  mutate(DEGadj=factor(DEGadj,levels=c("R>NR","R<NR","NS"))) 

dmp_sesame <- dmp_sesame %>% 
  group_by(DEGadj) %>% 
  arrange(Pval_StatusNR,-Est_StatusNR,.by_group=TRUE) %>% 
  # dplyr::rename("probeID"="Probe_ID") %>% 
  left_join(all_gene[,c("probeID","Probe_ID","Gene","Oncokb")])

write.table(dmp_sesame,glue('{output_dir}/epiborl_sesame_dmp_crpr_sdpd_adjpval0.01_est_0.25.txt'),quote=F,col.names=T,row.names=F,sep="\t")

# CR/PR/SD vs PD

se_dup2 <- SummarizedExperiment(assays=list(betas=as.matrix(betas_imput)), 
                               colData=(desc_sp[colnames(betas_imput),] %>% 
                                          dplyr::select(IDAT,Status2) %>% 
                                          mutate(Status=factor(Status2,levels=c("R","NR")))))

smry_dup2 <- DML(se_dup2, ~ Status2, BPPARAM = BiocParallel::MulticoreParam(num.cores))

saveRDS(smry_dup2, glue('{output_dir}/epiborl_sesame_dmp_crprsd_pd.Rds'))

gc(verbose = FALSE)

dmp_sesame2  <-  summaryExtractTest(smry_dup2)
colnames(dmp_sesame2) 
saveRDS(dmp_sesame2, glue('{output_dir}/epiborl_sesame_dmp_crprsd_pd_extract.Rds')) 

dmp_sesame2 <- dmp_sesame2 %>% 
  mutate(AdjPval_Status2NR=p.adjust(Pval_Status2NR,method="BH")) %>%
  mutate(DEGadj="NS") %>% 
  mutate(DEGadj=case_when(round(AdjPval_Status2NR,3) <= pval_thresh & Est_Status2NR >= est_thresh ~ "R<NR",
                          round(AdjPval_Status2NR,3) <= pval_thresh & Est_Status2NR <= -est_thresh ~ "R>NR",
                          TRUE ~ DEGadj)) %>%
  mutate(DEGadj=factor(DEGadj,levels=c("R>NR","R<NR","NS"))) 

dmp_sesame2 <- dmp_sesame2 %>% 
  group_by(DEGadj) %>% 
  arrange(Pval_Status2NR,-Est_Status2NR,.by_group=TRUE) %>% 
  # dplyr::rename("probeID"="Probe_ID") %>% 
  left_join(all_gene[,c("probeID","Probe_ID","Gene","Oncokb")])

write.table(dmp_sesame2,glue('{output_dir}/epiborl_sesame_dmp_crprsd_pd_adjpval0.01_est_0.25.txt'),quote=F,col.names=T,row.names=F,sep="\t")

# Join sesame/champ DMP

# CR/PR vs SD/PD

dm_dup <- full_join(dmpdr_champ %>% dplyr::select(probeID,Probe_ID,logFC,adj.P.Val,DEG,Gene,Oncokb),
                    dmp_sesame %>% dplyr::select(probeID,Probe_ID,Est_StatusNR,Pval_StatusNR,AdjPval_StatusNR,DEGadj:Oncokb),relationship = "many-to-many")

dm_dups <- dm_dup %>% filter(DEG!="NS" | DEGadj !="NS") # 175 probes

dm_dups <- dm_dups %>% mutate("Method"=case_when(DEG!="NS" & DEGadj !="NS"~"ChAMP/Sesame",
								  DEG!="NS" & DEGadj == "NS" ~ "ChAMP",
								  DEG=="NS" & DEGadj!="NS" ~ "Sesame"))

dm_dups %>% filter(DEG!=DEGadj) # 93 probes

write.table(dm_dups,glue('{output_dir}/epiborl_champ_sesame_dmp_crpr_sdpd_adjpval0.01_est_0.25_signif.txt'),quote=F,col.names=T,row.names=F,sep="\t")

## Volcano

gv <- read.table("~/GIT/methylation/EPICv2_manifest_short.txt",h=T,sep="\t") %>% 
  dplyr::select(-CHR,-MAPINFO,-Gene)

cpg_dru <- dm_dups %>% 
  distinct(probeID,.keep_all=TRUE) %>%
  left_join(gv,by="probeID")


all_cpg <- dmp_sesame %>% distinct(probeID,.keep_all=TRUE) %>% as.data.frame()
rownames(all_cpg) <- all_cpg$probeID 

# annocol"DEGadj"=c("NS"="gray55","R>NR"="steelblue","R<NR"="tomato3"))

temp <- all_cpg %>% 
  left_join(cpg_dru[,c("Probe_ID","DEG","Method")]) %>% 
  mutate(DEGadj2 = if_else(!is.na(DEG) & grepl("ChAMP",Method),DEG,DEGadj)) 

p <- ggplot(temp,aes(color=DEGadj2)) + 
  geom_point(aes(Est_StatusNR, -log10(round(AdjPval_StatusNR,3)))) + 
  theme_bw() +  scale_color_manual(values=annocol$DEGadj)  + 
  ggtitle("Best response: R=(CR/PR) vs NR=(SD/PD)") + 
  geom_label_repel(data=cpg_dru %>%  filter(!is.na(Oncokb)), 
                   aes(Est_StatusNR,-log10(round(AdjPval_StatusNR,3)),
                       color=DEG,
                       label=as.character(cpg_dru %>%filter(!is.na(Oncokb)) %>% pull(Gene))),
                   nudge_x=0.05,nudge_y=0.05,show.legend=F)


pdf(glue('{output_dir}/epiborl_volcano_oncokb_crpr_sdpd.pdf'),width=10,height=10)
p
dev.off()


## Heatmap 

bdmr <- betas_imput[dm_dups %>% distinct(probeID,.keep_all=TRUE) %>% pull(probeID),] %>% 
  as.data.frame() %>% 
  dplyr::select(desc_sp$IDAT) %>% 
  distinct() 

anno_col <- desc_sp[,c("Status","Best response")]
anno_row <- dm_dups %>% 
  filter(!is.na(Oncokb)) %>% 
  dplyr::select(probeID,Oncokb,Gene) %>% 
  column_to_rownames("probeID")

br_color_code <- c("CR"="#CA7B3A","PR"="#BE8C64","SD"="#157DC6","PD"="#427231")

g <- RColorBrewer::brewer.pal(nrow(anno_row), "Paired")
names(g) <- anno_row$Gene
 
ann_colors = list(Status = c(R="#66c2a5", NR="#fc8d62"),
                  Oncokb=c("Oncogene"="tomato3","TSG"="steelblue","Unknown"="grey80"), 
                  Gene=g, 
                  "Best response"=br_color_code)

ph <- pheatmap(bdmr,
         annotation_col=anno_col,
         fontsize_row = 6, fontsize_col=6, 
         annotation_colors = ann_colors, 
         annotation_row=anno_row)

pdf(glue('{output_dir}/epiborl_heatmap_dmp_crpr_sdpd.pdf'),width=10,height=12)
ph
dev.off()


### DMR between R (CR/PR) vs NR (SD/PD)

#using ChAMP 

source("~/EpiBORL/git/scripts/champ.DMR.R")


betas_dmr <- betas_imput[,rownames(desc_sp)]
rownames(betas_dmr) <- str_remove(rownames(betas_dmr),"\\_.*")


myDMR_3000 <- champ.DMR(beta = as.matrix(betas_dmr), 
                        pheno = as.character(desc_sp[colnames(betas_imput),"Status"]), 
                        compare.group==c("R","NR"), 
                        method="Bumphunter", 
                        arraytype = "EPICv2",
                        adjPvalDmr = pval_thresh_dmr,
                        cores=1,
                        minProbes=minProbes,
                        maxGap=3000, 
                        resultsDir=output_dir)

# saveRDS(myDMR_3000,glue('{output_dir}/epiborl_champ_dup_bumphunter_maxGap3000_dmr.RDS'))

data("probe.features.epicv2")
probe.features <- probe.features %>% rownames_to_column("probes") %>% 
  filter(probes %in% rownames(betas_imput)) %>% 
  column_to_rownames("probes")
Annot <- probe.features %>% arrange(CHR,MAPINFO) 

Anno <- map_dfr(1:nrow(myDMR_3000[[1]]),function(x){
		probe.features %>% 
    filter(CHR==myDMR_3000[[1]][x,"seqnames"] & MAPINFO >= myDMR_3000[[1]][x,"start"] &  MAPINFO <= myDMR_3000[[1]][x,"end"]) %>% 
    mutate(DMRindex=rownames(myDMR_3000[[1]][x,]),Size=myDMR_3000[[1]][x,"width"])})  %>%
  dplyr::select(DMRindex,everything())
	

DMRgenes <- Anno %>% dplyr::select(DMRindex,gene) %>% 
  distinct()%>% 
  filter(gene!="") %>% 
  arrange(as.numeric(gsub("DMR_","",DMRindex)))

DMRgenes <- DMRgenes %>% left_join(onco,by=c("gene"="Gene")) %>% 
  group_by(DMRindex) %>% 
  mutate(
    Oncokb = Oncokb %>% 
      na.omit() %>% 
      unique() %>% 
      str_flatten(collapse = "/") %>%
      na_if("")
  ) %>% 
  ungroup() %>% 
  distinct(DMRindex,.keep_all=TRUE) 

DMRgenes_onco <-  DMRgenes %>% 
  filter(!is.na(Oncokb))

myDMRt <-  myDMR_3000[[1]] %>% rownames_to_column("DMRindex") %>% dplyr::select(DMRindex,seqnames:width,p.valueArea) %>% left_join(DMRgenes)

write.table(myDMRt,glue('{output_dir}/epiborl_champ_dup_bumphunter_maxGap3000_dmr_anno_crprsd_pd.txt'),quote=F,col.names=T,row.names=F,sep="\t")

# >> head rmd_files_epic/epiborl_champ_dup_bumphunter_maxGap3000_dmr_anno.txt
# DMRindex	seqnames	start	end	width	p.valueArea	gene	Onco
# DMR_1	chr2	176109273	176122497	13224	1.88001729615912e-05	HOXD9	Unknown
# DMR_2	chr12	22331223	22333550	2327	0.0060207553909496	NA	NA


#### Using sesame

merge_sesame_dup <- DMR(se_dup, smry_dup, "StatusNR", platform="EPICv2") # merge CpGs to regions
saveRDS(merge_sesame_dup,glue('{output_dir}/epiborl_sesame_dup_dmr.RDS'))

merge_sesame_dup_signif <- merge_sesame_dup %>% 
  dplyr::filter(Seg_Pval_adj <= pval_thresh_dmr) 

dmr_se <- merge_sesame_dup_signif %>% 
  add_count(Seg_ID) %>% 
  arrange(desc(n)) %>% 
  dplyr::rename("Nb_Probe"="n") %>% 
  mutate(Seg_Size=Seg_End-Seg_Start)  

dmr_sep <- dmr_se %>% filter(Nb_Probe >= minProbes)

dmr_sep <- dmr_sep %>% 
  left_join(all_gene[,c("probeID","Probe_ID","Gene","Oncokb")],by=c("Probe_ID"="probeID"),relationship="many-to-many") %>% 
  arrange(Seg_Pval_adj)

dmr_sep2 <- dmr_sep %>% 
  group_by(Seg_ID) %>%
  mutate(Genes=paste(unique(Gene),collapse="/"),
    Onco = Oncokb %>% 
      na.omit() %>% 
      unique() %>% 
      str_flatten(collapse = "/") %>%
      na_if(""),
    Probes_NB=n(),Genes_NB=length(unique(Gene))
  ) %>% 
  ungroup() %>% 
  dplyr::select(Seg_ID:Seg_End,Seg_Size, Nb_Probe,Genes:Genes_NB) %>% 
  distinct()  #544 dmrs

Anno_sesame <- map_dfr(1:nrow(dmr_sep2),function(x){
	probe.features %>% 
    filter(CHR==unlist(dmr_sep2[x,"Seg_Chrm"]) & MAPINFO >= unlist(dmr_sep2[x,"Seg_Start"]) &  MAPINFO <= unlist(dmr_sep2[x,"Seg_End"])) %>% 
    mutate(DMRindex=unlist(dmr_sep2[x,"Seg_ID"]),Size=unlist(dmr_sep2[x,"Seg_Size"]))})  %>% 
  dplyr::select(DMRindex,everything()) 

write.table(dmr_sep,glue('{output_dir}/epiborl_sesame_dup_dmr_pval0.05_minProbes7_allProbes.txt'),quote=F,col.names=T,row.names=F,sep="\t")
write.table(dmr_sep2,glue('{output_dir}/epiborl_sesame_dup_dmr_pval0.05_minProbes7.txt'),quote=F,col.names=T,row.names=F,sep="\t")
write.table(Anno_sesame,glue('{output_dir}/epiborl_sesame_dup_dmr_pval0.05_minProbes7_annot.txt'),quote=F,col.names=T,row.names=F,sep="\t")


# Overlapping DMR between sesame and champ

dmr_champ <- myDMRt %>% dplyr::rename("Onco"="Oncokb")
dmr_sesame <- dmr_sep2 %>% as.data.frame

dmr_sesame_onco <- dmr_sesame %>% filter(!is.na(Onco))
dmr_champ_onco <- dmr_champ %>% filter(grepl("TSG|Oncogene|Both|Unknown",Onco)) 

dso <- unique(unlist(strsplit(dmr_sesame_onco$Genes,"/")))
dco <- unique(unlist(strsplit(dmr_champ_onco$gene,"/")))
gg <- names(table(c(dso,dco))[which(table(c(dso,dco))!=1)])

c_dmr_onco_sesame <- dmr_sesame_onco %>% filter(grepl(paste(gg[which(gg!="NA")],collapse="|"),Genes))
c_dmr_onco_champ <- dmr_champ_onco %>% filter(grepl(paste(gg[which(gg!="NA")],collapse="|"),gene))

dmr_sesamegr <-  GRanges(seqnames=dmr_sesame[,"Seg_Chrm"],
                         IRanges(start=dmr_sesame[,"Seg_Start"],
                                 end=dmr_sesame[,"Seg_End"]),
                         genes_sesame=dmr_sesame[,"Genes"],
                         onco_sesame=dmr_sesame[,"Onco"])

dmr_champgr <- GRanges(seqnames=dmr_champ[,"seqnames"],
                       IRanges(start=dmr_champ[,"start"],
                               end=dmr_champ[,"end"]),
                       genes_champ=dmr_champ[,"gene"],
                       onco_champ=dmr_champ[,"Onco"])

ov <- suppressMessages(ChIPpeakAnno::findOverlapsOfPeaks(dmr_sesamegr, dmr_champgr))

ovc <- as.data.frame(ov$overlappingPeaks[[1]])
rownames(ovc) <- 1:nrow(ovc)

colnames(ovc)[2:5] <- paste0(colnames(ovc)[2:5],"_sesame")
colnames(ovc)[10:13] <- paste0(colnames(ovc)[10:13],"_ChAMP")

ovc$Genes_sesame <- gsub("/NA/","",gsub("/NA$","",gsub("^NA/","",unlist(sapply(strsplit(ovc$genes_sesame,"/"),function(x){paste(x[grep("^ENSG|^SNORD|^LINC|Metazoa| RTAP",x,invert=T)],collapse="/")})))))

ovc$Genes_Champ <- gsub("/NA/","",gsub("/NA$","",gsub("^NA/","",unlist(sapply(strsplit(ovc$genes_champ,"/"),function(x){paste(x[grep("^ENSG|^SNORD|^LINC|Metazoa|^KRTAP",x,invert=T)],collapse="/")})))))

# Overlapping DMRs between ChAMP & sesame
ovc2 <- ovc %>% as.data.frame() %>% 
  dplyr::select(-peaks1,-strand,-peaks2,-genes_sesame,-genes_champ)

# Overlapping DMRs between ChAMP & sesame with Oncokb genes

ovc2_onco <- ovc2 %>% filter(!is.na(onco_sesame) & !is.na(onco_champ))

write.table(ovc2,glue('{output_dir}/epiborl_sesame_champ_overlap_dmr_pval0.05_minProbes7_annot.txt'),quote=F,col.names=T,row.names=F,sep="\t")




