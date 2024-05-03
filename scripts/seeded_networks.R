library(dplyr)
library(WGCNA)
library(biomaRt)


mart <- useMart(
  biomart = 'ENSEMBL_MART_ENSEMBL', 
  host    = 'grch37.ensembl.org',
  path    = '/biomart/martservice',
  dataset = 'hsapiens_gene_ensembl')

getinfo <- c("ensembl_gene_id","chromosome_name","start_position",
             "end_position","strand","band","gene_biotype", "external_gene_name")


# Visualize seeded networks 

seed_gene = "CCBL1" #KAT1
kyna.cor <- Kyna.genes_corr_list$KAT1 #list of correlated genes
datExpr <- datExpr.reg # processed bulk expression data
#sum(colnames(datExpr) == datMeta$X) 
#pos.cor1<- pos.cor ## save in another variable for now
#pos.cor <- pos.cor[-which(is.na(pos.cor$ID)),]

idx <- abs(kyna.cor$`R`) > 0.5 & kyna.cor$`FDR` < 0.05 
sum(idx) 

idx_genes <-data.frame(row.names =  kyna.cor$Gene[idx])
#rownames(idx_genes) <- UpdateSymbolList(symbols = rownames(idx_genes), verbose = T)
temp<- getBM(attributes = getinfo,filters=c("external_gene_name"),values= rownames(idx_genes),mart=mart)
idx_genes$ID <- temp$ensembl_gene_id[match(rownames(idx_genes), temp$external_gene_name)]

# Set up coordinates that the seeded networks will share
mat<- as.matrix(filter(as.data.frame(datExpr), rownames(datExpr)%in%idx_genes$ID))
ADJ1 <- adjacency(t(mat), type = "unsigned", power = 1)
network <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values= rownames(ADJ1),mart=mart)
colnames(ADJ1) = rownames(ADJ1) = network$external_gene_name[match(rownames(ADJ1), network$ensembl_gene_id)]
ADJ1[ADJ1 < 0.5] <- 0
graph.gene1 <- graph.adjacency(as.matrix(ADJ1), mode = "undirected", 
                               weighted = TRUE, diag = FALSE)
#unconnected_vertices <- which(degree(graph.gene1) == 0)

# Remove unconnected vertices
#graph.gene1 <- delete.vertices(graph.gene1, unconnected_vertices)

coords <- data.frame(layout.fruchterman.reingold(graph.gene1), 
                     gene = colnames(ADJ1))


## Plot the seeded network

ADJ1 <- adjacency(t(mat), type = "unsigned", power = 1)
ADJ1.neg <- cor(t(mat))
network <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values= rownames(ADJ1),mart=mart)
network.neg <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values= rownames(ADJ1.neg),mart=mart)
colnames(ADJ1) = rownames(ADJ1) = network$external_gene_name[match(rownames(ADJ1), network$ensembl_gene_id)]  #changed ensembl ids to gene names and added them here if necessary
colnames(ADJ1.neg) = rownames(ADJ1.neg) = network.neg$external_gene_name[match(rownames(ADJ1.neg), network.neg$ensembl_gene_id)] 
#ADJ1 <- ADJ1[-which(rownames(ADJ1)==""),-which(colnames(ADJ1)=="")]
#ADJ1.neg <- ADJ1.neg[-which(rownames(ADJ1.neg)==""),-which(colnames(ADJ1.neg)=="")]
#ADJ1 <- ADJ1[-which(is.na(rownames(ADJ1))),-which(is.na(colnames(ADJ1)))]
#ADJ1.neg <- ADJ1.neg[-which(is.na(rownames(ADJ1.neg))),-which(is.na(colnames(ADJ1.neg)))]

ADJ1[ADJ1 < 0.5] <- 0
ADJ1.neg[abs(ADJ1.neg) < 0.5] <- 0
graph.gene1 <- graph.adjacency(as.matrix(ADJ1), mode = "undirected", 
                               weighted = TRUE, diag = FALSE)
#unconnected_vertices <- which(degree(graph.gene1) == 0)

# Remove unconnected vertices
#graph.gene1 <- delete.vertices(graph.gene1, unconnected_vertices)

edgelist1 <- as_edgelist(graph.gene1)
edgelist.temp <- as_edgelist(graph.gene1)
edgelist1[, 1] <- match(edgelist1[, 1], coords$gene)
edgelist1[, 2] <- match(edgelist1[, 2], coords$gene)
edges1 <- data.frame(coords[edgelist1[, 1], 1:2], 
                     coords[edgelist1[, 2], 1:2], 
                     Network = "Gene", 
                     PCC = NA)
for (i in 1:nrow(edges1)) {
  edges1$PCC[i] <- ADJ1[edgelist.temp[i, 1], edgelist.temp[i, 2]]
  edges1$PCC.neg[i] <- ADJ1.neg[edgelist.temp[i, 1], edgelist.temp[i, 2]]
}
edges1$PCC.neg[edges1$PCC.neg > 0] <- 0
edges1$PCC[edges1$PCC.neg < 0] <- 0


#add groups
coords$group <- ifelse(coords$gene %in% kyna.cor$Gene[which(kyna.cor$R>0)], "POS", "NEG")
coords$keep <- ifelse(coords$gene %in% kyna.cor$Gene[order(abs(kyna.cor$R), decreasing = TRUE)[1:50]], "in","notin")

##### MAIN PLOT #######
ggplot() + 
  geom_segment(aes(x = X1, y = X2, xend = X1.1, yend = X2.1), color = 'darkred', 
               data = edges1, size = (edges1$PCC)^5, alpha = 0.1) +
  geom_segment(aes(x = X1, y = X2, xend = X1.1, yend = X2.1), color = 'darkblue', 
               data = edges1, size = abs(edges1$PCC.neg)^5, alpha = 0.1) +
  geom_point(aes(X1, X2), data = coords) +
  geom_point(aes(X1, X2), data = coords[which(coords$group=="POS"),][match(colnames(ADJ1), coords[which(coords$group=="POS"),]$gene), ], color = 'red', size=3) + 
  geom_point(aes(X1, X2), data = coords[which(coords$group=="NEG"),][match(colnames(ADJ1), coords[which(coords$group=="NEG"),]$gene), ], color = 'blue', size=3) +
  geom_point(aes(X1, X2), data = coords[coords$gene == "GOT2", ], alpha = 1, 
             color = 'black', size = 5.5) + 
  theme_classic() + labs(x = "", y = "") + 
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(1,1,1,1), "cm"), 
        rect = element_rect(fill = "transparent", colour = NA)) +
  #geom_text_repel(aes(x = X1, y = X2 - 2.0, label = gene), 
  #         data = coords[coords$gene == seed_gene, ], size = 5.0) + 
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[which(coords$keep=="in"),][match(rownames(idx_genes), coords[which(coords$keep=="in"),]$gene), ], 
                  size = 3.0, color = "black", 
                  box.padding = 0.5, force=2, max.overlaps = 30) + ggtitle(label = "KAT I-seeded co-expression network")
                  
#coords$keep <- ifelse(coords$gene %in% kyna.cor$Gene[order(kyna.cor$R, decreasing = TRUE)[1:10]], "in","notin")
### Got a beautiful graph. Change colours and visuals later. 
## Dont save EPS. It loses the edges. 

#saveRDS(list("edges"=edges1, "coords"= coords, "ADJ1"= ADJ1, "ADJ.neg"= ADJ1.neg, "idx_genes"=idx_genes), "KAT1_plotting_graph_objects.rds")
#save as pdf size 10x7
## Lowered threshold to 0.8 for KAT4. 

dev.off()










