## Packages

library(RColorBrewer)
library(ggplot2)
library(ape)
library(ggtree)
library(ggstance)
library(tidyverse)
library(ggrepel)


## Tree

em <- read_tsv("expr_em.txt")
head(em)
data <- em[,-1] %>% as.matrix()
colnames(data)<- paste0("DAP",c(0,2,4,6,8,12,16,22))
rownames(data) <- em$gene
data <- t(data)
hc <- hclust(dist(data))
den <- as.dendrogram(hc)

clus <- cutree(hc, 5)
d = data.frame(label=names(clus), member=factor(clus))
head(d)

p <- ggtree(as.phylo(hc),color = "#487AA1",size=1.5) %<+% d +
  geom_tiplab(aes(color = member), angle=0, hjust=-0.2,size=6,vjust=0.4)+
  geom_tippoint(aes(fill = member),size=8,shape=21)+
  scale_fill_manual(values = brewer.pal(6,"Set3")[-2])+
  scale_color_manual(values = brewer.pal(6,"Set3")[-2])+
  theme_tree2()+
  xlim(0,300)+
  theme(legend.position = "none")

p1 + geom_text(aes(label=node))
p1 <- p %>% flip(12,13) %>% flip(1,2) %>% flip(3,4)


## PCA

pca_res <- prcomp(data,center=TRUE)
summary(pca_res)

pca_res1 <- pca_res$x
pca_res1 <- as.data.frame(pca_res1)
pca_res1$ID <- rownames(pca_res1)
pca_res1$member = d$member

ggplot(pca_res1,aes(PC1,PC2))+
  geom_point(aes(fill=member),size=8,shape=21)+
  geom_text_repel(aes(label=ID,color=member),size=6)+
  scale_fill_manual(values = brewer.pal(6,"Set3")[-2])+
  scale_color_manual(values = brewer.pal(6,"Set3")[-2])+
  theme_classic(base_size = 16)+
  labs(x="PC1 (46.9%)",y="PC2 (24.6%)")+
  theme(legend.position = "none")

