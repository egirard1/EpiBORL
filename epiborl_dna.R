library(tidyverse)
library(ComplexHeatmap)
library(readxl)
library(RColorBrewer)

options(scipen=0)
set.seed(1)

in_dir <- "~/EpiBORL/git/scripts/rmd_files"
out_dir <- "~/EpiBORL/paper/scripts"

hpvkite_file <- "~/EpiBORL/epibol_rna_wes_summary.txt"

in_dir_files <- "/data/kdi_prod/.kdi/project_workspace_0/1975/acl/01.00/Share/From_Lyon_Wes_Mars_2025/vcf_annotated/reannot/oncodriver"

# Pathway gene file 
pathway_file <- "~/Pathways_DRAGON_07092021_v3.xlsx"

pathway <- read_excel(pathway_file) %>% 
  pivot_longer(cols=everything(),names_to="Pathway",values_to="Gene",values_drop_na = TRUE) %>% arrange(Pathway)

# Description file 

colum_desc <- c("ID","Best response","Primary site of cancer","Alcohols","Tobaccos", "Disease type","HPV_RNAseq_3.5%")

desc <- as_tibble(read.table(glue('{in_dir}/epiborl_sp_102024.txt'),h=T,sep="\t",check.names=F)) 

desc <- desc %>%   mutate("Alcohols"=if_else(Alcohol=="NON CONSUMER" ,"non consumer","consumer"),
                  "Tobaccos"=if_else(Tobacco=="NON CONSUMER","non consumer","consumer")) %>%
                   mutate(`Disease type`=case_when(grepl("Loc",disease_type) ~ "Local",grepl("Meta",disease_type) ~ "Metastatic", disease_type == "Both" ~ "Both", TRUE ~ NA )) 

hpvkite <- read.table(hpvkite_file,sep="\t",h=T,check.names=F) %>% distinct()
colnames(hpvkite) <- gsub("HPV","HPV16",colnames(hpvkite))

hpvkite <- hpvkite %>% mutate("HPV_RNAseq_3.5%"=if_else(round(score1,3)>=0.035 ,"P16pos","P16neg"))

desc <- desc %>% left_join(hpvkite[,c("ID","HPV_RNAseq_3.5%")]) %>% 
  filter(!is.na(`HPV_RNAseq_3.5%`)) # 149 patients

# Variant after oncoDriver 

lfile <- list.files(in_dir_files,pattern=".oncodriver.txt",full.names=TRUE,recursive=TRUE)
lfile <- lfile[grep("CDKN2A",lfile,invert=T)]

res <- bind_rows(map(lfile,function(file){

  name <- stringr::str_remove(basename(file),".oncodriver.txt") 

  tab <- read.table(file,header=TRUE,sep="\t",check.names=FALSE) %>% 
    separate(ONCODRIVER,sep=",",into=c("NM","oncokb","Decision")) %>% 
    dplyr::rename("Gene"="ANN[0].GENE") %>% 
    mutate(ID=name, Gene=as.character(Gene)) %>% 
    as_tibble()
  
})) %>% mutate(Oncogene=if_else(oncokb %in% c("oncogene","both"),"Yes","No"),
               TSG=if_else(oncokb %in% c("tsg","both"),"Yes","No") ) 


onco_gene <- res %>% dplyr::select(Gene,oncokb) 

res_all <- res %>% dplyr::select(ID,Gene,Decision,Oncogene,TSG) %>% 
  mutate(ID=factor(ID,levels=desc$ID))


# Per pathway

res_pathway <- map_dfc(unique(pathway$Pathway),function(pthw){
  res_all %>% filter(Gene %in% (pathway %>% filter(Pathway==pthw) %>%  pull(Gene))) %>% 
  dplyr::count(ID,.drop=FALSE, name=pthw) %>% 
  column_to_rownames("ID")
  })


# CNV after oncoDriver

lfile_cnv <- list.files(in_dir_files,pattern=".oncodriver.cnv.txt",full.names=TRUE,recursive=FALSE)

res_cnv <- bind_rows(map(lfile_cnv,function(file){

  name <- stringr::str_remove(basename(file),".oncodriver.cnv.txt") 

  tab <- read.table(file,header=TRUE,sep=",",check.names=FALSE) %>% 
    as_tibble()  %>% 
    dplyr::rename("Gene"="gene", "oncokb"="gene_type") %>% 
    mutate(ID=name, call=as.character(call), Gene=as.character(Gene), oncokb=as.character(oncokb),Geno=as.character(Geno),driver_status=as.character(driver_status))  

  tab_pik3CA <- tab %>% filter(call == "AMP" & Gene == "PIK3CA" & (loc.end-loc.start)<=10000000) %>% 
    dplyr::rename("Gene_pathway"="Gene") %>% 
    mutate(Gene=Gene_pathway) %>% 
    dplyr::select(ID,Gene,call, LOH, oncokb)
  
  # Remove DEL for chrX
  tab <- tab %>% filter(!(chrom==23 & call=="DEL"))

  tab <- tab %>% dplyr::select(ID,Gene,call, LOH, oncokb)

  # AMP if includes SOX2/PIK3CA/TP63 -> remove focal criteria & name it 3q26-28 (around 15% in HNSCC)

  tab <- tab %>% rename("Gene_pathway"="Gene") %>% 
    mutate(Gene = if_else(Gene_pathway %in% c("SOX2","PIK3CA","TP63") & call=="AMP","3q26-28",Gene_pathway))
  tab <- bind_rows(tab,tab_pik3CA)

  # Keep only AMP for CCND1 & remove FGF3/4/19 & SF3B2/LRP5 & RPS6KB2

  if(tab %>% filter(Gene=="CCND1" & call=="AMP") %>% nrow()==1){
    tab <- tab %>%  filter(!(Gene %in%c("FGF3","FGF4","FGF19","SF3B2","LRP5","RPS6KB2") & call=="AMP"))
  } 
    return(tab)

})) %>% mutate(Oncogene=if_else(oncokb %in% c("oncogene","both"),"Yes","No"),
               TSG=if_else(oncokb %in% c("tsg","both"),"Yes","No") ) %>% 
  mutate(ID=factor(ID,levels=desc$ID))


# Per pathway cnv 

res_cnv_pathway <- map_dfc(unique(pathway$Pathway),function(pthw){
  res_cnv %>% filter(Gene_pathway %in% (pathway %>% filter(Pathway==pthw) %>%  pull(Gene))) %>% 
  dplyr::count(ID,.drop=FALSE, name=pthw) %>% 
  column_to_rownames("ID")
  })


# Per pathway Variants + CNV

res_snv_cnv_pathway <- map_dfc(unique(pathway$Pathway),function(pthw){
  (bind_rows(res_all[,c("ID","Gene")],res_cnv[,c("ID","Gene_pathway")] %>% rename("Gene"="Gene_pathway")) %>% distinct()) %>% filter(Gene %in% (pathway %>% filter(Pathway==pthw) %>%  pull(Gene))) %>% 
  dplyr::count(ID,.drop=FALSE, name=pthw) %>% 
  column_to_rownames("ID")
  })

res_snv_cnv_pathway <- res_snv_cnv_pathway %>% 
  rownames_to_column("ID") %>% 
  pivot_longer(!ID,values_to="Nb",names_to="Pathway")  %>% 
  mutate(Status=if_else(Nb==0,"","Altered")) %>% 
  pivot_wider(!Nb,names_from = Pathway, values_from = Status) %>% 
  column_to_rownames("ID")


# Oncomat 

onco_gene_cnv <- res_cnv %>% dplyr::select(Gene,oncokb) %>% filter(!(Gene=="3q26-28"&oncokb=="both"))
res_cnv <- res_cnv %>% dplyr::select(-oncokb)

onco_res <- full_join(res_all,res_cnv,by=c("ID"="ID","Gene"="Gene","TSG"="TSG","Oncogene"="Oncogene")) %>% 
  group_by(ID,Gene,Decision) %>% 
  add_count(Decision,name="Count") %>% 
  ungroup() %>% 
  distinct() %>%
  mutate(Count = case_when(is.na(Decision) ~ NA_integer_, TRUE ~ Count)) %>% 
  arrange(ID,Gene)

onco_res <- onco_res %>% mutate(Final=case_when(
  Decision %in% c("conservative_inframe_deletion","disruptive_inframe_deletion") ~ "inframe indel",
  Decision %in% c("frameshift_variant") ~ "frameshift indel",
  Decision %in% c("missense_variant") ~ "missense",
  Decision %in% c("splice_acceptor_variant","splice_donor_variant") ~ "splicing",
  Decision %in% c("stop_gained") ~ "stop gain",
  (Decision %in% c("upstream_gene_variant") & Gene == "TERT") ~ "TERT",
  !is.na(call) ~ call,
  TRUE ~ NA
))

# filter on type of oncoKB genes - remove unknown

onco_res <- onco_res %>% filter(Oncogene!="No" | TSG!="No")

oncomat <- onco_res  %>%  
             dplyr::select(ID,Gene, Final) %>% 
             group_by(ID,Gene) %>% 
             mutate(t=paste(unique(Final), collapse=";")) %>% 
             dplyr::select(-Final) %>% 
             distinct()  %>% 
             pivot_wider(names_from = ID, values_from = t ,values_fill="") %>% 
             arrange(Gene) %>%
             column_to_rownames("Gene") %>%  
             ungroup() 


# Fill empty columns for samples without variants

if(ncol(oncomat)!=nrow(desc)){
  onco_empty <- data.frame(matrix(nrow = nrow(oncomat), ncol = length(desc[which(!desc$ID %in% colnames(oncomat)),"ID"] %>% unlist())))  %>% replace(is.na(.), "")
  colnames(onco_empty) <- desc[which(!desc$ID %in% colnames(oncomat)),"ID"] %>% unlist()
  rownames(onco_empty) <- rownames(oncomat)
  oncomat <- bind_cols(oncomat,onco_empty)
}

oncoGenes <- bind_rows(onco_gene,onco_gene_cnv) %>% 
  distinct() %>% 
  column_to_rownames("Gene")

oncoGenes <- oncoGenes  %>% 
  mutate(oncokb=case_when(oncokb=="tsg"~"TSG",
                          oncokb=="unknown" ~ "Unknown", 
                          oncokb == "oncogene" ~ "Oncogene", 
                          oncokb=="both" ~ "Both"),
         oncokb=factor(oncokb,levels=c("TSG","Oncogene","Both","Unknown")))

## Oncoprint ##

column_title <- "Oncoprint EpiBORL"

# Annotations

desc_anno <- desc %>% dplyr::select(colum_desc) %>% 
  filter(ID %in% colnames(oncomat)) %>% 
  distinct() %>% 
   mutate(`Best response`=factor(`Best response`,levels=c("CR","PR","SD","PD","NK")),
          Response=case_when(`Best response` %in% c("CR","PR") ~ "R",
                                           `Best response` %in% c("SD","PD") ~ "NR",
                                           TRUE ~ "NK"), 
                        Response=factor(Response,levels=c("R","NR","NK"))) %>%
  as.data.frame()

rownames(desc_anno) <- desc_anno$ID


pathway_color <- c("Altered"="steelblue")
resp_color_code <- c("R"="#D99E0F", "NR"="#728A69","NK"="#c4c4bf")
br_color_code <- c("CR"="#CA7B3A","PR"="#BE8C64","SD"="#157DC6","PD"="#427231","NK"="#c4c4bf")
HPV_RNAseq_3_color_code <- c("P16neg"="#315a72","P16pos"="#f3b00c","NA"="#c4c4bf")
prim_color_code <- c("Oral cavity"="#88AD60","Oropharynx"="#54950F","Hypopharynx"="#4D7226","Larynx"="#3B6B09","Other"="#2E5307")
alc_color_code <- c("consumer"="#b6ad95","non consumer"="#9e7f32")
tob_color_code <- c("consumer"="#b6ad95","non consumer"="#9e7f32")
disease_color_code <- c("Both"="#E89989","Local"="#E34122","Metastatic"="#9C311C","NA"="#c4c4bf")

annocol <- list("Best_response"=br_color_code,
                "Primary_tumor_site"=prim_color_code, 
                "HPV16"=HPV_RNAseq_3_color_code,
                "Alcohol"=alc_color_code,
                "Tobacco"=tob_color_code,
                "Disease_type"=disease_color_code,
                "Response"=resp_color_code)

nb_pc_oncogene <- 3 # keep only genes altered more than nb_pc%
nb_pc_both <- 3 # keep only genes altered more than nb_pc%
nb_pc_tsg <- 5 # keep only genes altered more than nb_pc%
nb_pc <- glue::glue('oncogene_{nb_pc_oncogene}_both_{nb_pc_both}_tsg_{nb_pc_tsg}')

# Gene order 
gene_order_value <- apply(oncomat,1,function(x){round(length(which(x!=""))*100/ncol(oncomat))})
gene_order <- c(names(which(gene_order_value[oncoGenes %>% filter(oncokb=="Oncogene") %>% rownames_to_column("Gene") %>% pull(Gene)]>=nb_pc_oncogene)),
names(which(gene_order_value[oncoGenes %>% filter(oncokb=="Both") %>% rownames_to_column("Gene") %>% pull(Gene)]>=nb_pc_both)),
names(which(gene_order_value[oncoGenes %>% filter(oncokb=="TSG") %>% rownames_to_column("Gene") %>% pull(Gene)]>=nb_pc_tsg)))

oncomat_op <- oncomat[gene_order,]

# order by response
sample_order <-  desc_anno %>% arrange(Response) %>% pull(ID)
tumor_order <- factor(desc_anno %>% arrange(Response) %>% pull(Response),levels=c("R","NR","NK"))


ha = HeatmapAnnotation(
    "Number alterations" = anno_oncoprint_barplot(),
    "Alcohol"= desc_anno[sample_order,"Alcohols"],
    "Tobacco"= desc_anno[sample_order,"Tobaccos"],
    "Best_response"= desc_anno[sample_order,"Best response"],
    "Primary_tumor_site"= desc_anno[sample_order,"Primary site of cancer"],
    "Disease_type"=desc_anno[sample_order,"Disease type"],
    "HPV16"= desc_anno[sample_order,"HPV_RNAseq_3.5%"],
    "Response"=desc_anno[sample_order,"Response"],
    col = annocol
)

alter_function = list(
                background = function(x, y, w, h)
                grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
                Altered = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Altered"], col = NA)),                 
                AMP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["AMP"], col = NA)),                                                      
                DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["DEL"], col = NA)),
                LOH = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["LOH"], col = NA)),
                `stop gain` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5,gp = gpar(fill = col["stop gain"], col = NA)),                                                                               
                `frameshift indel` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, gp = gpar(fill = col["frameshift indel"], col = NA)),                                                                            
                `inframe indel` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, gp = gpar(fill = col["inframe indel"], col = NA)),                                                                                
                missense = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, gp = gpar(fill = col["missense"], col = NA)),                                                                           
                splicing = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, gp = gpar(fill = col["splicing"], col = NA)),                                                                           
                `stop lost` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, gp = gpar(fill = col["stop_lost"], col = NA)),                                                                             
                `start lost` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, gp = gpar(fill = col["start_lost"], col = NA)),                                                            
                TERT = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, gp = gpar(fill = col["TERT"], col = NA)),
                multihits = function(x, y, w, h) grid.points(x, y, w*0.9, h*0.5, pch = 16)                                                       
                )

col_var = c("background"="#CCCCCC","frameshift indel" ="orange","inframe indel"="yellowgreen", "splicing"="pink", "missense"="forestgreen", "stop gain"="gold", "TERT"="darkblue", "multihits"="#041014")
col_cnv = c("background"="#CCCCCC","DEL" ="deepskyblue2", "AMP"="red", "LOH"="#49be25")
col <- c(col_var,col_cnv,pathway_color)

# For the project remove those genes  
gr <- c("EMSY","FRS2","PGR","FOXA1")
oncomat_op <- oncomat_op[setdiff(rownames(oncomat_op),gr),]

largeur_fixe = unit(35, "cm")

op <- oncoPrint(oncomat_op[,sample_order],
                heatmap_width = largeur_fixe,
   alter_fun = alter_function,
   col = col,
   top_annotation = ha,
   right_annotation = rowAnnotation(
        rbar = anno_oncoprint_barplot(
            axis_param = list(side = "bottom",labels_rot = 0))),
   column_title_gp = gpar(fontsize = 0, col = "transparent"),
   column_names_gp = gpar(fontsize = 6),
   row_names_gp = gpar(fontsize = 9), 
   pct_gp = gpar(fontsize = 8), 
   row_split = oncoGenes[match(rownames(oncomat_op), rownames(oncoGenes)),"oncokb"],
   column_split = tumor_order ,
   remove_empty_columns = TRUE, remove_empty_rows = FALSE, show_column_names = TRUE, alter_fun_is_vectorized = TRUE)

op

pdf(glue::glue('{out_dir}/files/epibORL_oncoKB_oncoprint_sup{nb_pc}_split_response.pdf'), width=16,height=10)
draw(op, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# Oncoprint per pathway

column_title_pathway <- "Oncoprint EpiBORL Pathway (snv, indels & cnv)"

res_sc_pathway <- t(res_snv_cnv_pathway)[,sample_order]

pathway_order_value <- apply(res_sc_pathway,1,function(x){round(length(which(x!=""))*100/ncol(res_sc_pathway))})

freq_onco_pathway <- data.frame("Pathway"=names(pathway_order_value),"Frequency"=pathway_order_value)

op_pathway <- oncoPrint(t(res_snv_cnv_pathway)[,sample_order],
   alter_fun = alter_function,
   heatmap_width = largeur_fixe,
   col = col,
   top_annotation = ha,
   right_annotation = rowAnnotation(
        rbar = anno_oncoprint_barplot(
            axis_param = list(side = "bottom",labels_rot = 0))),
   column_title_gp = gpar(fontsize = 0, col = "transparent"),
   column_names_gp = gpar(fontsize = 6),
   row_names_gp = gpar(fontsize = 10), 
   pct_gp = gpar(fontsize = 8), 
   column_split = tumor_order ,
   remove_empty_columns = TRUE, remove_empty_rows = FALSE, show_column_names = TRUE, alter_fun_is_vectorized = TRUE)

op_pathway


pdf(glue::glue('{out_dir}/files/epiborl_oncoKB_oncoprint_pathway_snv_cnv__split_response.pdf'), width=18,height=8)
draw(op_pathway, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

