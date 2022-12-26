## Load Package

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(BuenColors)
library(RColorBrewer)

## Load data

expr <- read_tsv("/path/GeneExpr.tsv")
modules_info <- read_tsv("/path/GeneExpr.tsv")
Module_selected <- c("tan","darkmagenta","white","sienna3","lightgreen","grey60")
keygene <- read_csv("/path/KeyGenes.tsv")
Module_selected_gene <- modules_info %>% filter(ID %in% Module_selected)

## Select Modules
module_gene_sel_expr_mt <- expr[,-1] %>% as.matrix() %>% apply(1,scale) %>% t()
rownames(module_gene_sel_expr_mt) <- expr$ID
tan_tb <- modules_info %>% filter(module == "tan")
darkmagenta_tb <- modules_info %>% filter(module == "darkmagenta")
white_tb <- modules_info %>% filter(module == "white")
sienna3_tb <- modules_info %>% filter(module == "sienna3")
lightgreen_tb <- modules_info %>% filter(module == "lightgreen")
grey60_tb <- modules_info %>% filter(module == "grey60")

## Plot fig5A

HaFUN <- function(mt){
  pos <- which(rownames(mt) %in% keygene$X3)
  key_gene_pos <- tibble(X3=rownames(mt)[pos],pos =pos) %>% 
    left_join(keygene)
  label=key_gene_pos$X2
  ha = rowAnnotation(
    foo = anno_mark(at = pos, labels = label)
  )
  return(ha)
}

ha1 <- HaFUN(module_gene_sel_expr_mt[tan_tb$ID,])
ha2 <- HaFUN(module_gene_sel_expr_mt[darkmagenta_tb$ID,])
ha3 <- HaFUN(module_gene_sel_expr_mt[white_tb$ID,])
ha4 <- HaFUN(module_gene_sel_expr_mt[sienna3_tb$ID,])
ha5 <- HaFUN(module_gene_sel_expr_mt[lightgreen_tb$ID,])
ha6 <- HaFUN(module_gene_sel_expr_mt[grey60_tb$ID,])

cols <- jdb_palette("wolfgang_basic")
col_fun = colorRamp2(seq(-2:2,length.out=9), cols)

ht1 <- Heatmap(module_gene_sel_expr_mt[tan_tb$ID,],
         cluster_columns = F,right_annotation = ha1,
         show_row_names = F,col = col_fun
 )


ht2 <- Heatmap(module_gene_sel_expr_mt[darkmagenta_tb$ID,],
         cluster_columns = F,right_annotation = ha2,
         show_row_names = F,col = col_fun
 )

ht3 <- Heatmap(module_gene_sel_expr_mt[white_tb$ID,],
         cluster_columns = F,right_annotation = ha3,
         show_row_names = F,col = col_fun
)

ht4 <- Heatmap(module_gene_sel_expr_mt[sienna3_tb$ID,],
         cluster_columns = F,right_annotation = ha4,
         show_row_names = F,col = col_fun
)

ht5 <- Heatmap(module_gene_sel_expr_mt[lightgreen_tb$ID,],
         cluster_columns = F,right_annotation = ha5,
         show_row_names = F,col = col_fun
 )

ht6 <- Heatmap(module_gene_sel_expr_mt[grey60_tb$ID,],
         cluster_columns = F,right_annotation = ha6,
         show_row_names = F,col = col_fun
 )

ht_list <- ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6
draw(ht_list,ht_gap = unit(1, "mm"),show_heatmap_legend = FALSE)

## Plot fig5B

atac_rna <- read_tsv("Path/atac_rna.tsv")

atac_rna_promter_mt <- atac_rna %>% 
  filter(str_detect(annotation,"Promoter")) %>% 
  filter(geneId %in% Module_selected_gene$ID) %>% 
  select(starts_with("DAP"),geneId,seqnames:end)

atac_rna_promter_mt1 <- atac_rna_promter_mt %>% apply(1,scale) %>% t()
pca_atac <- prcomp(t(atac_rna_promter_mt1[,-c(1:2)]),scale. = T,center = T)
pca_atac1 <- pca_atac$x %>%as.data.frame() %>%  rownames_to_column(.,var = "DAP")

pca_rna <- prcomp(t(module_gene_sel_expr_mt[Module_selected_gene,],scale. = T,center = T))
pca_rna1 <- pca_rna$x %>%as.data.frame() %>%  rownames_to_column(.,var = "DAP")

pca_atac1$DAP <- factor(pca_atac1$DAP,levels = pca_atac1$DAP)
pca_rna1$DAP <- factor(pca_rna1$DAP,levels = pca_rna1$DAP)

library(ggrepel)
cols1 <- c("#3288BD","#66C2A5","#ABDDA4","#FFFFBF","#FDAE61","#D53E4F")

atac_pca_p <- ggplot(pca_atac1,aes(PC1,PC2,fill=DAP))+
  geom_point(shape=21,size=6)+
  geom_text_repel(aes(label=DAP))+
  scale_fill_manual(values = cols1)+
  theme_bw(base_size = 16)+
  theme(legend.position = "none")+
  labs(x="PC1 (52.7%)",y="PC2 (19.8%)")+
  coord_cartesian(xlim = c(-230,130),ylim = c(-200,100))

rna_pca_p <- ggplot(pca_rna1,aes(PC1,PC2,fill=DAP))+
  geom_point(shape=21,size=6)+
  geom_text_repel(aes(label=DAP))+
  scale_fill_manual(values = cols1)+
  theme_bw(base_size = 16)+
  theme(legend.position = "none")+
  labs(x="PC1 (56.3%)",y="PC2 (21.2%)")+
  coord_cartesian(xlim = c(-150,100),ylim = c(-80,70))

atac_dist <- stats::dist(as.matrix(pca_atac1[,2:3]))
rna_dist <- stats::dist(as.matrix(pca_rna1[,2:3]))

atac_dist <- c(0,108.79402,175.37901,202.69558,30.02080,42.98150) %>% accumulate(sum)
rna_dist <- c(0,89.34644,106.28854,72.12761,58.29370,43.37830) %>% accumulate(sum)

atac_dist_norm <- (atac_dist-min(atac_dist))/(max(atac_dist)-min(atac_dist))*10
rna_dist_norm <- (rna_dist-min(rna_dist))/(max(rna_dist)-min(rna_dist))*10
pca_atac1$dist <- atac_dist_norm
pca_rna1$dist <- rna_dist_norm

rna_pca_p1 <- ggplot(pca_atac1)+
  geom_hline(yintercept = 0)+
  geom_point(aes(dist,0,fill=DAP),shape=21,size=6)+
  geom_text(aes(dist,-1,label=DAP))+
  theme_classic(base_size = 16)+
  scale_fill_manual(values = cols1)+
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  labs(x="",y="")+
  coord_cartesian(ylim = c(-2,1),xlim = c(-1,10))
  
atac_pca_p1 <- ggplot(pca_atac1)+
  geom_hline(yintercept = 0)+
  geom_point(aes(dist,0,fill=DAP),shape=21,size=6)+
  geom_text(aes(dist,-1,label=DAP))+
  theme_classic(base_size = 16)+
  scale_fill_manual(values = cols1)+
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  labs(x="",y="")+
  coord_cartesian(ylim = c(-2,1),xlim = c(-1,10))
  
## Plot fig5C

library(modelr)

grid <- data.frame(time = seq(0, 10, length.out = 500))

atac_modle_fun <- function(value){
  time = atac_dist_norm
  data = tibble(value = value[-c(1:2)],time =time)
  model = loess(value ~ time, data)
  predict <- grid %>% add_predictions(model)
  return(predict)
}

rna_modle_fun <- function(value){
  time = rna_dist_norm
  data = tibble(value = value[-c(1:2)],time =time)
  model = loess(value ~ time, data)
  predict <- grid %>% add_predictions(model)
  return(predict)
}

atac_modle_list <- apply(atac_data,1,atac_modle_fun)
rna_modle_list <- apply(rna_data,1,rna_modle_fun)
rna_modle_df <- rna_modle_list %>% reduce(bind_rows)
atac_modle_df <- atac_modle_list %>% reduce(bind_rows)

rna_modle_wide <-rna_modle_df %>% 
  mutate(rank = rep(1:500,10426), 
         ID = rep(rna_data_df$ID,each=500)) %>% 
  pivot_wider(names_from = rank, values_from = pred)

atac_modle_wide <-atac_modle_df %>% 
  mutate(rank = rep(1:500,31466), 
         ID = rep(c(atac_rna_promter$ID,atac_rna_distal$ID),each=500),
         seqnames = rep(c(atac_rna_promter$seqnames,atac_rna_distal$seqnames),each=500),
         start = rep(c(atac_rna_promter$start,atac_rna_distal$start),each=500),
         end = rep(c(atac_rna_promter$end,atac_rna_distal$end),each=500)
         ) %>% select(-time) %>% 
  pivot_wider(names_from = rank, values_from = pred)

rna_pca_data_gene <- rna_pca_data_gene %>% mutate(Atan = atan2(PC2,PC1)) %>% arrange(rank)
atac_pca_data_gene <- atac_pca_data_gene %>% mutate(Atan = atan2(PC1,PC2)) %>% arrange(rank)


cols2 <- rev(jdb_palette("algae_earth"))
cols3 <- jdb_palette("ocean_brick")
col_fun1 = colorRamp2(seq(-4,4,length.out=9), cols2)
col_fun2 = colorRamp2(seq(-4,4,length.out=9), cols3)

Heatmap(rna_pca_data_gene[,-rank],cluster_columns = F,
                               cluster_rows = F,show_row_names = F,
                               show_column_names = F,col = col_fun1,
                               border = TRUE)
Heatmap(atac_pca_data_gene[,-rank],cluster_columns = F,
                               cluster_rows = F,show_row_names = F,
                               show_column_names = F,col = col_fun2,
                               border = TRUE)
                               
## Network
library(clusterProfiler)
library(plyranges,include.only = c("as_granges","join_overlap_inner"))

## load footprint result
dap4 <- read_tsv("mpbs_bed/DAP4_mpbs.bed",col_names = F)
dap6 <- read_tsv("mpbs_bed/DAP6_mpbs.bed",col_names = F)
dap8 <- read_tsv("mpbs_bed/DAP8_mpbs.bed",col_names = F)
dap12 <- read_tsv("mpbs_bed/DAP12_mpbs.bed",col_names = F)
dap16 <- read_tsv("mpbs_bed/DAP16_mpbs.bed",col_names = F)
dap22 <- read_tsv("mpbs_bed/DAP22_mpbs.bed",col_names = F)

dap <- list(dap4=dap4,dap6=dap6,dap8=dap8,dap12=dap12,dap16=dap16,dap22=dap22)
dap1 <- dap %>% map(
  function(data){
    res1 <- data %>% filter(str_detect(X4,"MA"))
    res2 <- data %>% filter(X4 %in% motif_int)
    res <- bind_rowsWOX(res1,res2)
    return(res)
  }
)

dapx <- dap1 %>% reduce(bind_rows)
dapx <- dapx %>% select(X1:X4) %>% unique()

dapx_gr <- dapx %>% mutate(ID=paste(X1,X2,X3,sep = "_")) %>% 
  dplyr::rename(seqnames="X1",start="X2",end="X3") %>% 
  as_granges()
  
atac_promter_cluster_gr <- atac_promter_cluster %>% as_granges()
atac_promter_cluster_gr_overlap <- dapx_gr %>% join_overlap_inner(atac_promter_cluster_gr)
atac_promter_cluster_gr_tb <- as_tibble(atac_promter_cluster_gr_overlap) %>% 
  select(ID.x,cluster) %>% distinct()
  
t2g_motif <- dapx %>% transmute(Motif = X4,ID = paste(X1,X2,X3,sep = "_"))

atac_promter_cluster_motif <- compareCluster(
  ID.x~cluster,data=atac_promter_cluster_gr_tb,
  fun="enricher",TERM2GENE= t2g_motif,maxGSSize= 100000,
  minGSSize = 100)

atac_promter_cluster_motif_df <- as_tibble(atac_promter_cluster_motif)

atac_promter_cluster_motif_df_ma <- atac_promter_cluster_motif_df %>% 
  filter(str_detect(ID,"MA")) %>% 
  select(Cluster,ID,geneID) %>% 
  mutate(MotifID = str_sub(ID,1,8)) %>% 
  filter(MotifID %in% Module_group_gene_tf$MotifID)

atac_promter_ID <- atac_promter_cluster_motif_df_ma %>% 
  mutate(ID.x = str_split(geneID,"/")) %>% 
  unnest(ID.x) %>% 
  left_join(as_tibble(atac_promter_cluster_gr_overlap)) %>% 
  select(Cluster,ID,MotifID,ID.x,ID.y)  

atac_promter_net <- atac_promter_ID %>% 
  inner_join(Module_group_gene_tf[,c(1,2,3,6)]) %>% 
  dplyr::rename(TF_name = "Name",TF_class = "Class",
                TF_ID = "X2",TF_Motif= "MotifID") %>% 
  select(Cluster,TF_Motif,TF_name,TF_class,TF_ID,ID.x,X2=ID.y) %>% 
  inner_join(Module_group_gene_tf[,c(1,2,3,6)]) %>% 
  dplyr::rename(Target_name = "Name",Target_class = "Class",
                Target_ID = "X2",Target_Motif= "MotifID") 
  
atac_promter_ID_total <- t2g_motif %>%
  inner_join(atac_promter_cluster_gr_tb,by=c("ID"="ID.x"))%>% 
  dplyr::rename(ID.x="ID") %>% 
  inner_join(as_tibble(atac_promter_cluster_gr_overlap)) %>% 
  select(Cluster=cluster,ID=Motig,MotifID=Motig,ID.x,ID.y)
  
atac_promter_ID_total_ma <- atac_promter_ID_total %>% 
  filter(str_detect(ID,"MA")) %>% 
  mutate(MotifID = str_sub(ID,1,8))
  
atac_promter_net_total <- atac_promter_ID_total_ma %>% 
  inner_join(Module_group_gene_tf[,c(1,2,3,6)]) %>% 
  dplyr::rename(TF_name = "Name",TF_class = "Class",
                TF_ID = "X2",TF_Motif= "MotifID") %>% 
  select(Cluster,TF_Motif,TF_name,TF_class,TF_ID,ID.x,X2=ID.y) %>% 
  left_join(Module_group_gene_tf[,c(1,2,3,6)]) %>% 
  dplyr::rename(Target_name = "Name",Target_class = "Class",
                Target_ID = "X2",Target_Motif= "MotifID") 

library(qvalue)
expr_422 <- expr_mt %>% select(DAP4:DAP22) %>% as.matrix()
rownames(expr_422) <- rownames(expr_mt1)

corFun <- function(x,y){
  x_expr <- expr_422[x,]
  y_expr <- expr_422[y,]
  cor <- cor.test(x_expr ,y_expr,method = "spearman")$estimate
  return(cor)
}

network_total_ditinct <- network_total_ditinct %>% filter(TF_ID %in% rownames(expr_422) & Target_ID %in% rownames(expr_422))
for (i in 1:nrow(network_total_ditinct)){cor[i] <- corFun(network_total_ditinct$TF_ID[i], network_total_ditinct$Target_ID[i])}
network_total_ditinct$cor <- cor

network_total_ditinct_list <- network_total_ditinct %>% group_by(TF_ID) %>% group_split()

Getpvalue <- function(data){
  TF <- data$TF_ID[1]
  Gene0 <- setdiff(rownames(expr_422),data$Target_ID) %>% as.list()
  Cor0 <- map2(TF,Gene0,corFun) %>% unlist()
  Cor <- data$cor
  data$pvalue <- empPvals(stat = Cor, stat0 = Cor0)
  return(data)
}
Getfdr <- function(data){
  qobj <- qvalue(p =data$pvalue,pi0 = 1)
  lfdr <- qobj$lfdr
  data$fdr <- lfdr
  res <- data %>% filter(fdr <= 0.05)
  return(res)
}

network_total_ditinct_list1 <- network_total_ditinct_list %>% map(Getpvalue)
network_total_ditinct_list2 <- network_total_ditinct_list1 %>% 
  map(Getfdr)

network_new <- map(network_total_ditinct_list2,filter,fdr <= 0.05) %>% reduce(bind_rows)
network_total_new <- network_total %>% inner_join(network_new)

