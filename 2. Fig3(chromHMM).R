library(tidyverse)
library(plyranges)

## Load Data
cellFile <- dir(pattern = "*_15_segments.bed")
cell <- map(cellFile,read_tsv,col_names=F)
names(cell) <- c("cell1","cell2","cell3","cell4")
 colnames(cell$cell1) <- c("seqnames","start","end","state")
 colnames(cell$cell2) <- c("seqnames","start","end","state")
 colnames(cell$cell3) <- c("seqnames","start","end","state")
 colnames(cell$cell4) <- c("seqnames","start","end","state")
 
## Convert Data to Granges
cell_gr <- cell %>% map(as_granges)
## Intersect Granges
cell12 <- join_overlap_intersect(cell_gr$cell1, cell_gr$cell2, minoverlap = 2) %>% as_tibble()
cell13 <- join_overlap_intersect(cell_gr$cell1, cell_gr$cell3, minoverlap = 2) %>% as_tibble()
cell14 <- join_overlap_intersect(cell_gr$cell1, cell_gr$cell4, minoverlap = 2) %>% as_tibble()
cell23 <- join_overlap_intersect(cell_gr$cell2, cell_gr$cell3, minoverlap = 2) %>% as_tibble()
cell24 <- join_overlap_intersect(cell_gr$cell2, cell_gr$cell4, minoverlap = 2) %>% as_tibble()
cell34 <- join_overlap_intersect(cell_gr$cell3, cell_gr$cell4, minoverlap = 2) %>% as_tibble()

cell_change <- list(cell12=cell12,cell13=cell13,cell14=cell14,cell23=cell23,cell24=cell24,cell34=cell34)

statFun <- function(data){
  data_stat <- data %>% group_by(state.x,state.y) %>% summarise(base=sum(width-1))
  data_stat1 <- data_stat %>% group_by(state.x)%>% summarise(base_total=sum(base))
  data_stat2 <- data_stat %>% group_by(state.x) %>% filter(state.x != state.y)%>% summarise(base_var=sum(base))
  res1 <- tibble(state=data_stat1$state.x ,base_var=data_stat2$base_var,base_total=data_stat1$base_total) %>% mutate(ratio=base_var/base_total)
  res2 <- sum(data_stat2$base_var)/sum(data_stat1$base_total)
  res <- list(res1=res1,res2=res2)
  return(res)
}

## Calculate percentage
cell_change_stat_data <- bind_rows(cell_change_stat$cell12$res1,
  cell_change_stat$cell13$res1,cell_change_stat$cell14$res1,
  cell_change_stat$cell23$res1,cell_change_stat$cell24$res1,
  cell_change_stat$cell34$res1) %>% 
  mutate(change=rep(c("cell12","cell13","cell14","cell23","cell24","cell34"),each=15)) %>% 
  inner_join(state_iden)

cell_change_stat_data$state <- factor(cell_change_stat_data$state,
  levels=c("E13","E9","E7","E12","E11","E10","E6","E14","E15","E8","E4","E5","E3","E1","E2"))
cell_change_stat_data$state <- factor(
  cell_change_stat_data$state,
  levels=c("E5","E4","E7","E9",
           "E6","E8","E10",
           "E2","E3","E1",
           "E11","E12",
           "E13","E15","E14"))  


cell_change_stat_data$change <- factor(cell_change_stat_data$change,
  levels=c("cell34","cell24","cell23","cell14","cell13","cell12"))

## plot
cols <- c(jdb_palette("Darjeeling")[1:3],jdb_color_map(c("NK")),jdb_palette("Darjeeling")[5])
cols <- c(c("#be6f7c", "#fac699", "#94ad90","#447593","#66446f"))
  

p1 <- ggplot(cell_change_stat_data)+
  geom_tile(aes(state,change,alpha=ratio,fill=group),color="black")+
  theme_void()+scale_fill_manual(values = cols)

p2 <- ggplot(cell_change_stat_data)+
  geom_boxplot(aes(state,ratio,fill=group))+
  theme_classic()+scale_fill_manual(values = cols)+
  theme(legend.position = "none")

p3 <- ggplot(cell_change_stat_summary)+
  geom_col(aes(change, percentage),width = 0.7,fill=jdb_color_map(c("GMP-A")),color="black")+
  theme_classic()+coord_flip()
