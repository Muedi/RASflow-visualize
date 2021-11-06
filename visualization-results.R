library(tidyverse)
library(yaml)
library(DESeq2)
library(biomaRt)
library(tximport)
library(mygene)
library(hash)
  ### heatmap of topgenes
library("genefilter")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("pheatmap")
library("RColorBrewer")
library(viridis)
library("ggbeeswarm")
library(EnhancedVolcano)
library(openxlsx)
norm.path <- "/home/max/projects/NGS/neutrophiles/RASflow/output/neutrophiles_combined/trans/dea/countGroup"
dea.path <- "/home/max/projects/NGS/neutrophiles/RASflow/output/neutrophiles_combined/trans/dea/DEA/gene-level"
out.path <- "/home/max/projects/NGS/neutrophiles/RASflow/output/neutrophiles_combined/trans/dea/visualization"

yaml.file <- yaml.load_file('configs/config_main_combined.yaml')
yaml.file.old <- yaml.load_file('configs/config_main_neutro_old.yaml')


# extract the information from the yaml file
project <- yaml.file$PROJECT  # project name
project.old <- yaml.file.old$PROJECT  # project name
dea.tool <- yaml.file$DEATOOL  # tool used for DEA
quant.path <- file.path(yaml.file$FINALOUTPUT, project, "trans/quant")
quant.path.old <- file.path(yaml.file.old$FINALOUTPUT, project.old, "trans/quant")
gene.level <- yaml.file$GENE_LEVEL  # whether to do gene-level analysis
controls <- yaml.file$CONTROL  # all groups used as control
treats <- yaml.file$TREAT  # all groups used as treat, should correspond to control
filter.need <- yaml.file$FILTER$yesOrNo
pair.test <- yaml.file$PAIR
meta.file <- "/home/max/projects/NGS/neutrophiles/RASflow/configs/metadata_combined_sets.tsv"
meta.file.old <- "/home/max/projects/NGS/neutrophiles/RASflow/configs/metadata_old.tsv"
ENSEMBL <- yaml.file$ENSEMBL
dataset <- yaml.file$EnsemblDataSet
dea.path <- file.path(yaml.file$FINALOUTPUT, project, "trans/dea")
fdr.thr <- yaml.file$FDR  # threshold of FDR/adjusted P-value for significantlly differentially expressed genes
num.control <- length(controls)  # number of comparisons that the user wants to do
num.treat <- length(treats)  # should equals to num.control


if (num.control != num.treat) {
  message("Error: Control groups don't match with treat groups!")
  message("Please check config_dea.yaml")
  quit(save = 'no')
}

num.comparison <- num.control

convert.id2symbol <- function(gene.id) {
  gene.symbol <- gene.id  # initialize the gene symbol with the gene id
  
  # it may happen that no symbol can be found for any id. In that case, "queryMany" will throw an error
  # so "try" is used here to take care of that error
  try({
    gene.symbol.all <- queryMany(gene.id, scopes = 'ensembl.gene', fields = 'symbol')
    
    h <- hash()
    for (i in 1:nrow(gene.symbol.all)) {
      query <- gene.symbol.all$query[i]
      symbol <- gene.symbol.all$symbol[i]
      if (has.key(query, h)) {  # if there's duplicate for the same query
        h[[query]] <- paste(hash::values(h, keys = query), symbol, sep = ', ')
      } else {
        if (is.na(symbol)) {  # if there's no hit for the query, keep the original id
          h[[query]] <- query
        } else {
          h[[query]] <- symbol
        }
      }
    }
    
    for (i in c(1:length(gene.symbol))) {
      gene.symbol[i] <- h[[gene.id[i]]]
    }
  })
  
  return(gene.symbol)
}



### load salomns quants as deseq2 object.
### the original quant files from Salmon
meta.data <- read.csv(meta.file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
meta.data.old <- read.csv(meta.file.old, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
samples.all <- meta.data$sample
group.all <- meta.data$group
subject.all <- meta.data$subject
samples.all.old <- meta.data.old$sample
group.all.old <- meta.data.old$group
subject.all.old <- meta.data.old$subject
# get complete counts
samples <- factor(samples.all) 
samples.old <- factor(samples.all.old) 
  ### import quantification as txi
  # files <- file.path(quant.path, samples, "quant.sf")
  # names(files) <- samples
  # load tx2gene
  output.path <- file.path(yaml.file$FINALOUTPUT, project, "trans/dea")
  load(file.path(output.path, "countGroup", 'tx2gene.RData'))
  # noVersion because ensemble is used 
  files.noVersion <- file.path(quant.path, samples, "quant_noVersion.sf")
  names(files.noVersion) <- samples
  files.noVersion.old <- file.path(quant.path.old, samples.old, "quant_noVersion.sf")
  names(files.noVersion.old) <- samples.old
  # load data fro salom, gene level
  txi_all <- tximport(files.noVersion, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
  txi_all_old <- tximport(files.noVersion.old, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")


# all.dea.tibble <- tibble()

# joined_table will hold TPM vals and the padj/lfc and an indicator if diff expressed of all comparisons by the end. 
# before the loop only tpms are in there and it is used for plotting of the heatmap and PCA of all samples by expression value.
joined_table <- as.data.frame(txi_all$counts) # length sclaed tpm since indicated in tximport function
joined_table <- tibble(rownames_to_column(joined_table, "ensembl"))   
# filter rows, with only counts below 2
joined_table <- joined_table %>% filter_at(vars(-ensembl), any_vars(. > 2))
joined_table.2 <- joined_table
# old 
joined_table.old <- as.data.frame(txi_all_old$counts) # length sclaed tpm since indicated in tximport function
joined_table.old <- tibble(rownames_to_column(joined_table.old, "ensembl"))   
# filter rows, with only counts below 2
joined_table.old <- joined_table.old %>% filter_at(vars(-ensembl), any_vars(. > 2))

# deseq2 normalization 

colData = data.frame(samples.all, subject.all, group.all)
design <- model.matrix(~group.all)
dds <- DESeqDataSetFromTximport(txi_all, colData = colData, design = design)
dds <- estimateSizeFactors(dds)
normalized_counts_deseq2 <- counts(dds, normalized=TRUE)
normalized_counts_deseq2 <- as_tibble(normalized_counts_deseq2, rownames="ensembl")



# pca of all samples
library(gridExtra)
library(grid)
library(sjmisc)
library(org.Mm.eg.db)

df_pca <- prcomp(joined_table %>% column_to_rownames("ensembl") %>% 
                  dplyr::select( c("45__PU1_Neu_1712_S45_R1_001",
                  "46__PU1_Neu_1715_S46_R1_001",
                  "48__PU1_Neu_1718_S48_R1_001",
                  "41__PU1_WT_1731_S41_R1_001",
                  "42__PU1_WT_1691_S42_R1_001",
                  "43__PU1_WT_1735_S43_R1_001",
                  "44__PU1_WT_1734_S44_R1_001")) %>%
                  replace(is.na(.), 0) %>% rotate_df()) # make ensembl the rownames again, replkace all nas with 0 and finally transpose matrix
df_out <- as.data.frame(df_pca$x)
# head(df_out)
# add group
df_out$group <- (meta.data %>% as_tibble() %>% column_to_rownames("sample") )[rownames(df_out),"group"]
# plt
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,label=substring(rownames(df_out), 1, 10) )) + 
      geom_point() + 
      geom_label_repel(hjust="inward", nudge_y = - 20000, max.overlaps=60) +
      xlab(percentage[1]) + ylab(percentage[2])
ggsave(filename = file.path(out.path, 'PCA_Prim_length-scaled-tpm.png'))


# heatmap of gene expression
topVarGenes <- head(order(rowVars(as.matrix(joined_table %>% column_to_rownames("ensembl"))), decreasing = TRUE), 1000)
  mat  <- joined_table[ topVarGenes, ] %>% column_to_rownames("ensembl")
  #mat  <- mat - rowMeans(mat)
  # annotate mat
  ens.str <- rownames(mat)
  # sym.str <- mapIds(org.Mm.eg.db,
  #                   keys=ens.str,
  #                   column="SYMBOL",
  #                   keytype="ENSEMBL",
  #                   multivals="first")
  sym.str <- convert.id2symbol(ens.str)
  sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
  rownames(mat) <- sym.str
  group.anno <- meta.data %>% column_to_rownames("sample") %>% dplyr::select("group")

# mat_breaks <- seq(min(mat), max(mat), length.out = 10)

png(file = file.path(out.path, 'top-gene-heatmap_all_samples.png'), width=3300, height=6000, res=300, title = 'top genes by variance')
pheatmap(log2(mat + 1), 
          annotation_col = group.anno, 
          color=inferno(20), 
          main="Top 1000 Genes, log transformed lengthScaledTPM") 
dev.off() 

# heatmap normalized by row means and color scheme by quantile
# heatmap of gene expression
topVarGenes <- head(order(rowVars(as.matrix(joined_table %>% column_to_rownames("ensembl"))), decreasing = TRUE), 1000)
  mat  <- joined_table[ topVarGenes, ] %>% column_to_rownames("ensembl")
  mat  <- mat - rowMeans(mat)
  # annotate mat
  ens.str <- rownames(mat)
  # sym.str <- mapIds(org.Mm.eg.db,
  #                   keys=ens.str,
  #                   column="SYMBOL",
  #                   keytype="ENSEMBL",
  #                   multivals="first")
  sym.str <- convert.id2symbol(ens.str)
  sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
  rownames(mat) <- sym.str
  group.anno <- meta.data %>% column_to_rownames("sample") %>% dplyr::select("group")

# mat_breaks <- seq(min(mat), max(mat), length.out = 10)


# qualntile breaks for equally distributed color
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

colnames(mat) <- paste("n" , colnames(mat), sep="")
mat_breaks <- quantile_breaks(unlist(mat), n = 11)

rownames(group.anno) <- colnames(mat)

png(file = file.path(out.path, 'top-gene-heatmap_all_samples_quantile_breaks.png'),width=3300, height=6000, res=300,  title = 'top genes by variance')
pheatmap(mat, 
          annotation_col = group.anno, 
          breaks = mat_breaks,
          color=inferno(length(mat_breaks) - 1), 
          main="Top 1000 Genes, lengthScaledTPM normalized by rowMeans, quantile breaks in color scheme") 
dev.off() 

old_expression_clusters <- read.xlsx("/mnt/c/Users/masprang/Desktop/Projects/Neutrophil-PU.1-project/mRNA_Hoxis_naiive_Ca_clusters_genelist_JF.xlsx", 1, startRow = 2)

clust_genes <- old_expression_clusters[c("Cluster", "Name")]
clust_genes <- clust_genes[clust_genes$Cluster != "none",]

# # avergae expression of all samples/groups
# exp.plot <- joined_table %>% gather(key="sample", value="value", -ensembl) %>% 
#     left_join(meta.data %>% dplyr::select(sample, group), by="sample") %>% 
#     filter(value>0) %>%
#     mutate(value=log2(value)) %>% 
#     mutate(sample=as.factor(sample)) %>% 
#     ggplot(aes(x=value, y=group, color=group)) + 
#     geom_boxplot() +
#     labs(title="Avg. Gene Expression in groups",
#          x="Log2 Gene Expression", y="Group"
#     ) +
#     coord_flip() +
#     theme_minimal()
# ggsave(filename = file.path(out.path, 'avg_expression_group.png'))


## TODO: 
# get TPMs from RPKMs of the old data to compare directly. 

# get log normalized TPMs:
logTPM <- joined_table %>% mutate(across(-ensembl, ~ log(.x + 1)))
logTPM.old <- joined_table.old %>% mutate(across(-ensembl, ~ log(.x + 1)))

# compute z-scores
zscore<- function(x){
    z<- (x - mean(x)) / sd(x)
    return(z)
}
# older data zscores driectly from RPKMs
# zscores_old_rpkm <- as_tibble(old_expression_clusters) %>% mutate(across(c(-Cluster, -Name), ~ zscore(.x)))
zscores <- logTPM %>% mutate(across(-ensembl, ~ zscore(.x)))
zscores$Symbol <- convert.id2symbol(zscores$ensembl)
# old data
zscores.old <- logTPM.old %>% mutate(across(-ensembl, ~ zscore(.x)))
zscores.old$Symbol <- convert.id2symbol(zscores.old$ensembl)
# # avg expr over all
# exp.plot <- zscores %>% gather(key="sample", value="value",  c(-Symbol, -ensembl)) %>% 
#     left_join(meta.data %>% dplyr::select(sample, group), by="sample") %>% 
#     mutate(sample=as.factor(sample)) %>% 
#     ggplot(aes(x=value, y=sample, color=group)) + 
#     geom_boxplot() +
#     labs(title="Avg. Gene Expression in groups",
#          x="z-score", y="sample"
#     ) +
#     coord_flip() +
#     theme_minimal()
# ggsave(filename = file.path(out.path, 'avg_expression_sample.png'))
# Cluster I
exp.plot <- zscores %>% filter(Symbol %in% clust_genes[clust_genes$Cluster == "I","Name"]) %>%
    gather(key="sample", value="value", c(-Symbol, -ensembl)) %>% 
    left_join(meta.data %>% dplyr::select(sample, group), by="sample") %>% 
    mutate(sample=as.factor(sample)) %>% 
    ggplot(aes(x=value, y=group, color=group)) + 
    geom_boxplot() +
    labs(title="Avg. Gene Expression in groups",
         x="Z-score", y="Group"
    ) +
    coord_flip() +
    theme_minimal()
ggsave(filename = file.path(out.path, 'avg_expression_score_clustI.png'))

# Cluster II
exp.plot <- zscores %>% filter(Symbol %in% clust_genes[clust_genes$Cluster == "II","Name"]) %>%
    gather(key="sample", value="value", c(-Symbol, -ensembl)) %>% 
    left_join(meta.data %>% dplyr::select(sample, group), by="sample") %>% 
    mutate(sample=as.factor(sample)) %>% 
    ggplot(aes(x=value, y=group, color=group)) + 
    geom_boxplot() +
    labs(title="Avg. Gene Expression in groups",
         x="Z-score", y="Group"
    ) +
    coord_flip() +
    theme_minimal()
ggsave(filename = file.path(out.path, 'avg_expression_score_clustII.png'))

# Cluster III
exp.plot <- zscores %>% filter(Symbol %in% clust_genes[clust_genes$Cluster == "III","Name"]) %>%
    gather(key="sample", value="value", c(-Symbol, -ensembl)) %>% 
    left_join(meta.data %>% dplyr::select(sample, group), by="sample") %>% 
    mutate(sample=as.factor(sample)) %>% 
    ggplot(aes(x=value, y=group, color=group)) + 
    geom_boxplot() +
    labs(title="Avg. Gene Expression in groups",
         x="Z-score", y="Group"
    ) +
    coord_flip() +
    theme_minimal()
ggsave(filename = file.path(out.path, 'avg_expression_score_clustIII.png'))

# Cluster IV
exp.plot <- zscores %>% filter(Symbol %in% clust_genes[clust_genes$Cluster == "IV","Name"]) %>%
    gather(key="sample", value="value", c(-Symbol, -ensembl)) %>% 
    left_join(meta.data %>% dplyr::select(sample, group), by="sample") %>% 
    mutate(sample=as.factor(sample)) %>% 
    ggplot(aes(x=value, y=group, color=group)) + 
    geom_boxplot() +
    labs(title="Avg. Gene Expression in groups",
         x="Z-score", y="Group"
    ) +
    coord_flip() +
    theme_minimal()
ggsave(filename = file.path(out.path, 'avg_expression_score_clustIV.png'))
# old data plotted 
# Cluster I
exp.plot <- zscores.old %>% filter(Symbol %in% clust_genes[clust_genes$Cluster == "I","Name"]) %>%
    gather(key="sample", value="value", c(-Symbol, -ensembl)) %>% 
    left_join(meta.data.old %>% dplyr::select(sample, group, subject), by="sample") %>% 
    mutate(sample=as.factor(sample)) %>% 
    ggplot(aes(x=value, y=group, color=subject)) + 
    geom_boxplot() +
    labs(title="Avg. Gene Expression in groups",
         x="Z-score", y="Group"
    ) +
    coord_flip() +
    theme_minimal()
ggsave(filename = file.path(out.path, 'avg_expression_score_clustI.old.png'))

# Cluster II
exp.plot <- zscores.old %>% filter(Symbol %in% clust_genes[clust_genes$Cluster == "II","Name"]) %>%
    gather(key="sample", value="value", c(-Symbol, -ensembl)) %>% 
    left_join(meta.data.old %>% dplyr::select(sample, group, subject), by="sample") %>% 
    mutate(sample=as.factor(sample)) %>% 
    ggplot(aes(x=value, y=group, color=subject)) + 
    geom_boxplot() +
    labs(title="Avg. Gene Expression in groups",
         x="Z-score", y="Group"
    ) +
    coord_flip() +
    theme_minimal()
ggsave(filename = file.path(out.path, 'avg_expression_score_clustII.old.png'))

# Cluster III
exp.plot <- zscores.old %>% filter(Symbol %in% clust_genes[clust_genes$Cluster == "III","Name"]) %>%
    gather(key="sample", value="value", c(-Symbol, -ensembl)) %>% 
    left_join(meta.data.old %>% dplyr::select(sample, group, subject), by="sample") %>% 
    mutate(sample=as.factor(sample)) %>% 
    ggplot(aes(x=value, y=group, color=subject)) + 
    geom_boxplot() +
    labs(title="Avg. Gene Expression in groups",
         x="Z-score", y="Group"
    ) +
    coord_flip() +
    theme_minimal()
ggsave(filename = file.path(out.path, 'avg_expression_score_clustIII.old.png'))

# Cluster IV
exp.plot <- zscores.old %>% filter(Symbol %in% clust_genes[clust_genes$Cluster == "IV","Name"]) %>%
    gather(key="sample", value="value", c(-Symbol, -ensembl)) %>% 
    left_join(meta.data.old %>% dplyr::select(sample, group, subject), by="sample") %>% 
    mutate(sample=as.factor(sample)) %>% 
    ggplot(aes(x=value, y=group, color=subject)) + 
    geom_boxplot() +
    labs(title="Avg. Gene Expression in groups",
         x="Z-score", y="Group"
    ) +
    coord_flip() +
    theme_minimal()
ggsave(filename = file.path(out.path, 'avg_expression_score_clustIV.old.png'))

# comparison to proteomics expression
diff_express_protein_old <- read.xlsx("/mnt/c/Users/masprang/Desktop/Projects/Neutrophil-PU.1-project/DEPs_pathway_PU1_BM_neutros_HA.xlsx", 1)
diff_express_protein_old <- as_tibble(diff_express_protein_old)
diff_express_protein_old <- diff_express_protein_old %>% separate_rows(Gene_name)
#diff_express_protein_old <- diff_express_protein_old %>% mutate()


# combined_degs holds top X of all comparisons 
combined_degs = tibble()
combined_degs_1000 = tibble()
combined_degs_100 = tibble()
plot_spi1 = tibble()
plot_ets2 = tibble()
plot_junb = tibble()
plot_ccl3 = tibble()

## create the DESeqDataSet
# for subset samples (comparisons)
for (i in 1:length(controls)) {
#for (i in 2) {

  control <- controls[i]
  treat <- treats[i]
  samples <- factor(samples.all[c(which(group.all == control), which(group.all == treat))]) 
  ### import quantification as txi
  # files <- file.path(quant.path, samples, "quant.sf")
  # names(files) <- samples
  # load tx2gene
  output.path <- file.path(yaml.file$FINALOUTPUT, project, "trans/dea")
  load(file.path(output.path, "countGroup", 'tx2gene.RData'))
  # noVersion because ensemble is used 
  files.noVersion <- file.path(quant.path, samples, "quant_noVersion.sf")
  names(files.noVersion) <- samples
  # load data fro salom, gene level
  txi <- tximport(files.noVersion, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")

  subject <- factor(subject.all[c(which(group.all == control), which(group.all == treat))])
  group <- relevel(factor(group.all[c(which(group.all == control), which(group.all == treat))]), ref = control)
  colData = data.frame(samples, subject, group)
  design <- model.matrix(~group)
  dds <- DESeqDataSetFromTximport(txi, colData = colData, design = design)

  # Filtering
  if (filter.need) {
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
  }
  ## specify the control group
  dds$group <- relevel(dds$group, ref = control)
  colnames(dds) <- dds@colData$samples
  # count data
  countdata <- round(assays(dds)[["counts"]])

  # compare transformations

  vsd <- vst(dds, blind = FALSE)
  rld <- rlog(dds, blind = FALSE)

  dds <- estimateSizeFactors(dds)

  df <- bind_rows(
    as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
          mutate(transformation = "log2(x + 1)"),
    as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
    as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
    
  colnames(df)[1:2] <- c("x", "y")  
  # lvls <- c("log2(x + 1)", "vst", "rlog")
  # df$transformation <- factor(df$transformation, levels=lvls)
  #######
  # png(file = file.path(out.path, paste('nornalization-comparison', control, '_', treat, '.png', sep = '')), title = 'comparison of transformations for counts')
  # ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  #   coord_fixed() + facet_grid( . ~ transformation)  
  # dev.off()


  # # distances 
  # sampleDists <- dist(t(assay(vsd)))
  # sampleDistMatrix <- as.matrix( sampleDists )
  # rownames(sampleDistMatrix) <- vsd@colData@rownames

  # colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

  # png(file = file.path(out.path, paste('sample-distances', control, treat, '.png', sep = '_')), width=750, height=500,  title = 'Sample distances of variance stabilizing transformed counts')
  # pheatmap(sampleDistMatrix,
  #         labels_row = rownames(sampleDistMatrix),
  #         show_colnames = F,
  #         clustering_distance_rows = sampleDists,
  #         clustering_distance_cols = sampleDists,
  #         col = colors)
  # dev.off()

  # # pca
  # pcaData <- plotPCA(vsd, intgroup = "group", returnData = T)
  # percentVar <- round(100 * attr(pcaData, "percentVar"))

  # ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
  #   geom_point(size =3) +
  #   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  #   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  #   coord_fixed() +
  #   ggtitle("PCA with VST data")
  # ggsave(filename = file.path(out.path, paste('PCA_', control, '_', treat, '.png', sep = '')))


  # topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
  # mat  <- assay(vsd)[ topVarGenes, ]
  # mat  <- mat - rowMeans(mat)
  # # annotate mat
  # ens.str <- rownames(mat)
  # sym.str <- mapIds(org.Mm.eg.db,
  #                   keys=ens.str,
  #                   column="SYMBOL",
  #                   keytype="ENSEMBL",
  #                   multivals="first")
  # sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
  # rownames(mat) <- sym.str
  # group.anno <- as.data.frame(colData(vsd)[, c("samples", "group")])
  # group.anno <- group.anno["group"]

  # library(grid)
  # ## Edit body of pheatmap:::draw_colnames, customizing it to your liking
  # # draw_colnames_45 <- function (coln, ...) {
  # #   m = length(coln)
  # #   x = (1:m)/m - 1/2/m
  # #   grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
  # #       hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
  # # }

  # # ## 'Overwrite' default draw_colnames with your own version 
  # # assignInNamespace(x="draw_colnames", value="draw_colnames_45",
  # # ns=asNamespace("pheatmap"))

  # png(file = file.path(out.path, paste('top-gene-heatmap_', control, '_', treat, '.png', sep = '')),width=750, height=750,  title = 'Top genes by variance')
  # pheatmap(mat, annotation_col = group.anno) 
  # dev.off()

  ## perform DEA
  dds <- DESeq(dds)
  res <- results(dds)

  # enrichment/ annotation
  topGene <- rownames(res)[which.min(res$padj)]

  # png(file = file.path(out.path, paste('countplot_topgene_', control, '_', treat, '.png', sep = '')), title = 'PCA of variance stabilizing transformed counts')
  # plotCounts(dds, gene = topGene, intgroup=c("group"))
  # dev.off()

  # annotate results
  ens.str <- rownames(res)
  # sym.str <- mapIds(org.Mm.eg.db,
  #                     keys=ens.str,
  #                     column="SYMBOL",
  #                     keytype="ENSEMBL",
  #                     multiVals="first")
  sym.str <- convert.id2symbol(ens.str)
  # sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
  res$symbol <- sym.str

  # vulcano <- 
  #png(file = file.path(out.path, paste('volcano_plot_', control, '_', treat, '.png', sep = '')), width=750, height=500,  title = paste('Volcano plot: ', control, treat, sep = ' '))
  EnhancedVolcano(res,
    lab = res$symbol,
    x = 'log2FoldChange',
    title = paste('Volcano plot: ', control, treat, sep = ' '),
    y = 'pvalue')
  #dev.off()
  ggsave(filename = file.path(out.path, paste('volcano_plot_', control, '_', treat, '.png', sep = '')), width=7, height=7)
  #as.pdf(volcano, width = 9, height = 6, scaled = TRUE, file = file.path(out.path, paste('volcano_plot_', name.control, '_', name.treat, '.png', sep = '')))

  # ordered after pval
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  
  top20 <- as_tibble(resOrdered[1:20,] , rownames="ensembl")
  if ( plyr::empty(combined_degs)   ) {
    combined_degs <- top20 %>% mutate(combination=paste(control, treat, sep = "_"))
  } else {
    combined_degs <- combined_degs %>%
      bind_rows(top20 %>% mutate(combination=paste(control, treat, sep = "_")))
  }
  # top 1000 degs
  top1000 <- as_tibble(resOrdered[1:1500,] , rownames="ensembl")
  if ( plyr::empty(combined_degs_1000)   ) {
    combined_degs_1000 <- top1000 %>% mutate(combination=paste(control, treat, sep = "_"))
  } else {
    combined_degs_1000 <- combined_degs_1000 %>%
      bind_rows(top1000 %>% mutate(combination=paste(control, treat, sep = "_")))
  }
  # top 100 degs
  top100 <- as_tibble(resOrdered[1:100,] , rownames="ensembl")
  if ( plyr::empty(combined_degs_100)   ) {
    combined_degs_100 <- top100 %>% mutate(combination=paste(control, treat, sep = "_"))
  } else {
    combined_degs_100 <- combined_degs_100 %>%
      bind_rows(top100 %>% mutate(combination=paste(control, treat, sep = "_")))
  }

  # heatmap 1 vs 1
  mat  <- joined_table.2 %>% 
                    filter(ensembl %in% top100$ensembl) %>% 
                    column_to_rownames("ensembl") %>%
                    dplyr::select( (meta.data %>% filter(group %in% c(control, treat)) )$sample )
  #mat  <- mat - rowMeans(mat)
  # annotate mat
  ens.str <- rownames(mat)
  # sym.str <- mapIds(org.Mm.eg.db,
  #                   keys=ens.str,
  #                   column="SYMBOL",
  #                   keytype="ENSEMBL",
  #                   multivals="first")
  sym.str <- convert.id2symbol(ens.str)
  sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
  rownames(mat) <- sym.str
  group.anno <- meta.data %>% filter(group %in% c(control, treat)) %>% column_to_rownames("sample") %>% dplyr::select("group")
  # mat_breaks <- seq(min(mat), max(mat), length.out = 10)

  png(file = file.path(out.path, paste('top100_degs',control, treat ,'heatmap.png', sep='_')),width=3300, height=3600, res=300)
  pheatmap(log2(mat + 1), 
          annotation_col = group.anno, 
          color=inferno(20), 
          cluster_cols=F,
          main=paste(control, treat, "Top 100 differential genes, log transformed lengthScaledTPM", sep=" ") 
          )
  dev.off() 

  tibb.dea <- as_tibble(resOrdered , rownames="ensembl")
  tibb.dea <- tibb.dea[!is.na(tibb.dea$padj),]
  # tibb.dea <- tibb.dea[tibb.dea$padj < 0.05,]
  # tibb.dea <- tibb.dea[abs(tibb.dea$log2FoldChange) > 2,]
  to_join <- tibb.dea %>% dplyr::select("ensembl", "log2FoldChange", "padj")
  to_join <- to_join %>% mutate("DEG_lfc>2" = padj < 0.05 & (log2FoldChange>2 | log2FoldChange < -2) )
  to_join <- to_join %>% rename_with(~ paste(., control, treat, sep="_"), -"ensembl")
  #tibb.dea <- tibb.dea %>% column_to_rownames("symbol")

  # colnames(tibb.dea)
  # if (is_empty(all.dea.tibble)) {
  #   all.dea.tibble <- tibb.dea
  #   # all.dea.tibble <- all.dea.tibble %>% rename_with(~ paste(., control, treat, sep="_"), -"symbol")
  #   }
  # else {
  #   all.dea.tibble  <- all.dea.tibble %>%
  #   full_join(tibb.dea, by = "symbol")
  # }
  joined_table <- joined_table %>%
      full_join(to_join, by = "ensembl")
  
  # plot fold changes in genomic space
  # takes quite long
  # library(apeglm)
  # resGR <- lfcShrink(dds, coef=c("groupKO.PU"), type="apeglm", format="GRanges")
}
library(openxlsx)

ens.str <- joined_table$ensembl
sym.str <- convert.id2symbol(ens.str)

joined_table <- joined_table %>% mutate(symbol = sym.str) 

write.xlsx(joined_table, file=file.path(out.path, "master.table.xlsx"), row.names=F, overwrite=T)

tpm <- joined_table[c(samples.all, "ensembl", "symbol")]

heatmap_input <- tpm %>% filter(ensembl %in% combined_degs$ensembl)
combis <- combined_degs %>% filter(ensembl %in% heatmap_input$ensembl) %>% dplyr::select(combination)
combined_degs$symbol <- convert.id2symbol(combined_degs$ensembl)
# get rid of duplications and (hopefully retain info on both combis
test <- combined_degs %>% group_by(ensembl) %>% summarise(combination_summed=paste0(combination, sep = "|"))
combined_degs$combination_summed <- test$combination_summed
distinct_ensembl <-  combined_degs %>% 
     distinct(ensembl, .keep_all=T)
# add toi heatmap input
x <- distinct_ensembl %>% dplyr::select("combination", "ensembl")
heatmap_input <- heatmap_input %>% left_join(x, by="ensembl")
heatmap_input <- heatmap_input[order(heatmap_input$combination),]
# heatmap of gene expression
  mat  <- heatmap_input %>% column_to_rownames("ensembl") %>% dplyr::select(-symbol, -combination) 
  #mat  <- mat - rowMeans(mat)
  # annotate mat
  ens.str <- rownames(mat)
  # sym.str <- mapIds(org.Mm.eg.db,
  #                   keys=ens.str,
  #                   column="SYMBOL",
  #                   keytype="ENSEMBL",
  #                   multivals="first")
  sym.str <- convert.id2symbol(ens.str)
  sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
  rownames(mat) <- sym.str
  group.anno <- meta.data %>% column_to_rownames("sample") %>% dplyr::select("group")
  anno.genes <- distinct_ensembl %>% column_to_rownames("symbol") %>% dplyr::select("combination")
# mat_breaks <- seq(min(mat), max(mat), length.out = 10)

png(file = file.path(out.path, 'degs-all-combis-heatmap_all_samples.png'),width=3300, height=3600, res=300)
pheatmap(log2(mat + 1), 
          annotation_col = group.anno, 
          annotation_row = anno.genes,
          color=inferno(20), 
          cluster_cols=F,
          main="Top 20 differential genes, log transformed lengthScaledTPM") 
dev.off() 

# without clustering
# heatmap of gene expression
  mat  <- heatmap_input %>% column_to_rownames("ensembl") %>% dplyr::select(-symbol, -combination) 
  #mat  <- mat - rowMeans(mat)
  # annotate mat
  ens.str <- rownames(mat)
  # sym.str <- mapIds(org.Mm.eg.db,
  #                   keys=ens.str,
  #                   column="SYMBOL",
  #                   keytype="ENSEMBL",
  #                   multivals="first")
  sym.str <- convert.id2symbol(ens.str)
  sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
  rownames(mat) <- sym.str
  group.anno <- meta.data %>% column_to_rownames("sample") %>% dplyr::select("group")
  anno.genes <- distinct_ensembl %>% column_to_rownames("symbol") %>% dplyr::select("combination")
# mat_breaks <- seq(min(mat), max(mat), length.out = 10)

png(file = file.path(out.path, 'degs-all-combis-noclust-heatmap_all_samples.png'),width=3300, height=3600, res=300)
pheatmap(log2(mat + 1), 
          annotation_col = group.anno, 
          annotation_row = anno.genes,
          color=inferno(20), 
          cluster_rows=F,
          cluster_cols=F,
          main="Top 20 differential genes, non-clustered, log transformed lengthScaledTPM") 
dev.off() 


# only Hox
# combined_degs_wt <- combined_degs %>% filter( grepl("WT_", combined_degs$combination)) 
combined_degs_wt <- combined_degs %>% filter(!grepl("Prim", combined_degs$combination)) 
# combined_degs_wt <- combined_degs_100 %>% filter(grepl("Hox-KO_v2_Hox-2KOs", combined_degs_100$combination)) 
heatmap_input <- normalized_counts_deseq2 %>% 
                      filter(ensembl %in% combined_degs_wt$ensembl) %>%
                      # column_to_rownames("ensembl") %>% 
                      dplyr::select(-c("45__PU1_Neu_1712_S45_R1_001",
                      "46__PU1_Neu_1715_S46_R1_001",
                      "48__PU1_Neu_1718_S48_R1_001",
                      "41__PU1_WT_1731_S41_R1_001",
                      "42__PU1_WT_1691_S42_R1_001",
                      "43__PU1_WT_1735_S43_R1_001",
                      "44__PU1_WT_1734_S44_R1_001") )
combis <- combined_degs_wt %>% filter(ensembl %in% heatmap_input$ensembl) %>% dplyr::select(combination)
combined_degs_wt$symbol <- convert.id2symbol(combined_degs_wt$ensembl)
# get rid of duplications and (hopefully retain info on both combis
test <- combined_degs_wt %>% group_by(ensembl) %>% summarise(combination_summed=paste0(combination, sep = "|"))
combined_degs_wt$combination_summed <- test$combination_summed
distinct_ensembl <-  combined_degs_wt %>% 
     distinct(ensembl, .keep_all=T)
# add toi heatmap input
x <- distinct_ensembl %>% dplyr::select("combination", "ensembl")
heatmap_input <- heatmap_input %>% left_join(x, by="ensembl")
heatmap_input <- heatmap_input[order(heatmap_input$combination),]
# heatmap of gene expression
  mat  <- heatmap_input %>% column_to_rownames("ensembl") %>% dplyr::select( -combination) 
  #mat  <- mat - rowMeans(mat)
  # annotate mat
  ens.str <- rownames(mat)
  # sym.str <- mapIds(org.Mm.eg.db,
  #                   keys=ens.str,
  #                   column="SYMBOL",
  #                   keytype="ENSEMBL",
  #                   multivals="first")
  sym.str <- convert.id2symbol(ens.str)
  sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
  rownames(mat) <- sym.str
  group.anno <- meta.data %>% column_to_rownames("sample") %>% dplyr::select("group")
  anno.genes <- distinct_ensembl %>% column_to_rownames("symbol") %>% dplyr::select("combination")
# mat_breaks <- seq(min(mat), max(mat), length.out = 10)

png(file = file.path(out.path, 'degs-heatmap_Hox_deseq2.png'),width=3300, height=3600, res=300,  title = 'top genes by variance')
pheatmap(log(mat +1) , 
          annotation_col = group.anno, 
          annotation_row = anno.genes,
          color=inferno(20), 
          cluster_cols=T,
          cluster_rows=T,
          main="Top 20 differential genes, normalized DESEQ2") 
dev.off() 


# get clusters print map of 1000 degs for the wt combis
# only wt vs everything
combined_degs_wt <- combined_degs_1000 %>% filter( grepl("WT_", combined_degs_1000$combination)) 
combined_degs_wt <- combined_degs_wt %>% distinct(ensembl, .keep_all=T)
heatmap_input <- tpm %>% filter(ensembl %in% combined_degs_wt$ensembl) 
combis <- combined_degs_wt %>% filter(ensembl %in% heatmap_input$ensembl) %>% dplyr::select(combination)
combined_degs_wt$symbol <- convert.id2symbol(combined_degs_wt$ensembl)
# get rid of duplications and (hopefully retain info on both combis
test <- combined_degs_wt %>% group_by(ensembl) %>% summarise(combination_summed=paste0(combination, sep = "|"))
combined_degs_wt$combination_summed <- test$combination_summed
# distinct_ensembl <-  combined_degs_wt %>% distinct(symbol, .keep_all=T)
     #distinct(ensembl, .keep_all=T)
# add toi heatmap input
x <- combined_degs_wt %>% dplyr::select("combination", "ensembl")
heatmap_input <- heatmap_input %>% left_join(x, by="ensembl")
heatmap_input <- heatmap_input[order(heatmap_input$combination),]
# heatmap of gene expression
  mat  <- heatmap_input %>% column_to_rownames("ensembl") %>% dplyr::select(-symbol, -combination) 
  #mat  <- mat - rowMeans(mat)
  # annotate mat
  # ens.str <- rownames(mat)
  # # sym.str <- mapIds(org.Mm.eg.db,
  # #                   keys=ens.str,
  # #                   column="SYMBOL",
  # #                   keytype="ENSEMBL",
  # #                   multivals="first")
  # sym.str <- convert.id2symbol(ens.str)
  # sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
sym.str <- heatmap_input$symbol
sym.str[is.na(sym.str)] <- heatmap_input$ensembl
sym.str <- make.unique(sym.str)
rownames(mat) <- sym.str
group.anno <- meta.data %>% column_to_rownames("sample") %>% dplyr::select("group")
anno.genes <- distinct_ensembl %>% column_to_rownames("symbol") %>% dplyr::select("combination")
# mat_breaks <- seq(min(mat), max(mat), length.out = 10)
png(file = file.path(out.path, 'degs-wt-combis-heatmap_all_samples_top1500.png'),width=3300, height=3600, res=300,  title = 'top genes by variance')
wt_clust <- pheatmap(log2(mat + 1), 
          annotation_col = group.anno, 
          annotation_row = anno.genes,
          color=inferno(20), 
          cluster_cols=F,
          main="Top 1000 differential genes, log transformed lengthScaledTPM") 
wt_clust
dev.off() 
clusters_wt_combis <- cutree(wt_clust$tree_row, k = 3)

# Venn diagramm of old an New clusters? Dot plot better? 
clusts_list_new <- list(
  "1" = names(clusters_wt_combis[clusters_wt_combis == 1]),
  "2" = names(clusters_wt_combis[clusters_wt_combis == 2]),
  "3" = names(clusters_wt_combis[clusters_wt_combis == 3])#,
  # "4" = names(clusters_wt_combis[clusters_wt_combis == 4]),
  # "5" = names(clusters_wt_combis[clusters_wt_combis == 5]),
  # "6" = names(clusters_wt_combis[clusters_wt_combis == 6])
)

clusts_list_old <- list(
  "I" = clust_genes[clust_genes$Cluster == "I", "Name"],
  "II" = clust_genes[clust_genes$Cluster == "II", "Name"],
  "III" = clust_genes[clust_genes$Cluster == "III", "Name"],
  "IV" = clust_genes[clust_genes$Cluster == "IV", "Name"]
)

# "35__HoxPU1_Neu_v2_1_S35_R1_001"             "36__HoxPU1_Neu_v2_2_S36_R1_001"             "37__HoxPU1_Neu_v2_3_S37_R1_001"             "38__HoxPU1_Neu_creEts2_c3_KO1_1_S38_R1_001" "39__HoxPU1_Neu_creEts2_c3_KO1_2_S39_R1_001" "40__HoxPU1_Neu_creEts2_c3_KO1_3_S40_R1_001"
#  "41__PU1_WT_1731_S41_R1_001"                 "42__PU1_WT_1691_S42_R1_001"                 "43__PU1_WT_1735_S43_R1_001"                 "44__PU1_WT_1734_S44_R1_001"                 "45__PU1_Neu_1712_S45_R1_001"                "46__PU1_Neu_1715_S46_R1_001"               
#  "48__PU1_Neu_1718_S48_R1_001"                "ensembl"                                    "symbol"      

# heatmap of old Clusters in the new data.
heatmap_input <- tpm %>% filter(symbol %in% clust_genes$Name) %>% dplyr::select(c(-"38__HoxPU1_Neu_creEts2_c3_KO1_1_S38_R1_001",
                                                                                  -"39__HoxPU1_Neu_creEts2_c3_KO1_2_S39_R1_001",
                                                                                  -"40__HoxPU1_Neu_creEts2_c3_KO1_3_S40_R1_001",
                                                                                  -"45__PU1_Neu_1712_S45_R1_001",
                                                                                  -"46__PU1_Neu_1715_S46_R1_001",
                                                                                  -"48__PU1_Neu_1718_S48_R1_001"))
  mat  <- heatmap_input %>% column_to_rownames("ensembl") %>% dplyr::select(-symbol) 
  #mat  <- mat - rowMeans(mat)
  # annotate mat
  # ens.str <- rownames(mat)
  # # sym.str <- mapIds(org.Mm.eg.db,
  # #                   keys=ens.str,
  # #                   column="SYMBOL",
  # #                   keytype="ENSEMBL",
  # #                   multivals="first")
  # sym.str <- convert.id2symbol(ens.str)
  # sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
sym.str <- heatmap_input$symbol
sym.str[is.na(sym.str)] <- heatmap_input$ensembl
sym.str <- make.unique(sym.str)
rownames(mat) <- sym.str
group.anno <- meta.data %>% column_to_rownames("sample") %>% dplyr::select("group")
anno.genes <- clust_genes %>% as_tibble() %>% column_to_rownames("Name") 
# mat_breaks <- seq(min(mat), max(mat), length.out = 10)
png(file = file.path(out.path, 'degs-old-clusters-heatmap_all_samples_top1500.png'),width=3300, height=3600, res=300,  title = 'Old cluster genes')
wt_clust <- pheatmap(log2(mat + 1), 
          annotation_col = group.anno, 
          annotation_row = anno.genes,
          color=inferno(20), 
          cluster_cols=F,
          main="Old Cluster genes, log transformed lengthScaledTPM") 
wt_clust
dev.off() 
# library(ggvenn)
# save_venn <- function(clust_old) {
#   x <- c(clusts_list_old[clust_old], clusts_list_new)
#   ggvenn(
#     x, 
#     fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#     stroke_size = 0.5, set_name_size = 4
#     )
#   ggsave(filename = file.path(out.path, paste('Venn_', "CLust-",clust_old,"-vs-1-2-3", '.png', sep = '')))
# }
# x <- c(clusts_list_old["III"], clusts_list_new[4:6])
#   ggvenn(
#     x, 
#     fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#     stroke_size = 0.5, set_name_size = 4
#     )
#   ggsave(filename = file.path(out.path, paste('Venn_', "CLust-","III","-vs-1-2-3", '.png', sep = '')))


# library(eulerr)
lists <- c(clusts_list_new, clusts_list_old)
# lists_eulerr <- lapply(lists, length)

# # Get the combinations of names of list elements
# nms <- combn( names(lists) , 2 , FUN = paste0 , collapse = "&" , simplify = FALSE )
# # Make the combinations of list elements
# ll <- combn( lists , 2 , simplify = FALSE )
# # Intersect the list elements
# list.intersections <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )
# # Output with names
# list.intersections <- setNames( list.intersections , nms )
# lists_eulerr <- append(lists_eulerr , list.intersections)
# fit_euler <- euler(unlist(lists_eulerr), shape = "ellipse")
# png(file = file.path(out.path, 'euler-plot.png'),width=2000, height=2000, res=300,  title = 'Euler diagram of all clusters (old: roman numbers)')
# plot(fit_euler, quantities = list(type = "counts"))
# dev.off() 

library(UpSetR)

upset_plot <- upset(fromList(lists), 
      nsets = length(lists), nintersects = 60,
      mainbar.y.label = "Intersections", sets.x.label = "DE Genes in Cluster X",
      mb.ratio = c(0.6, 0.4),
      text.scale = c(1.3, 1.3, 1, 1, 1, 0.75)
)
# ggsave(file.path("output","all proj UpSet plot.png"))
png(file = file.path(out.path, 'upset-plot.png'),width=2000, height=2000, res=300,  title = 'Upset plot of all clusters (old: roman numbers)')
upset_plot
dev.off()

# function for countplots
save_combined_countplots <- function(tpm, sym) {
  countplot <- tpm[tpm$symbol == sym,] %>% dplyr::select(-ensembl, -symbol)
  # countplot <- countplot + 0.00001
  countplot <- countplot %>% rotate_df() %>% dplyr::rename(count = V1) %>% rownames_to_column("sample")
  countplot <- countplot %>% left_join(meta.data[c("sample", "group")], by="sample")

  # combined countplots
  ggplot(countplot, aes(x=group, y=count, color=sample)) +
      #scale_y_log10() +  
      geom_beeswarm(size = 3,cex = 3) + ggtitle(paste0("Countplot: ", sym)) + theme(axis.text.x = element_text(angle=45, hjust=1))
  ggsave(filename = file.path(out.path, paste('countplot_', sym, '_all_samples_lsTPM', '.png', sep = '')))
}

# function for countplots with subgroups
save_combined_countplots_subset <- function(tpm, sym, list_of_samples, samples_text) {
  countplot <- tpm[tpm$symbol == sym,] %>% dplyr::select(all_of(list_of_samples), -ensembl, -symbol)
  # countplot <- countplot + 0.00001
  countplot <- countplot %>% rotate_df() %>% dplyr::rename(count = V1) %>% rownames_to_column("sample")
  countplot <- countplot %>% left_join(meta.data[c("sample", "group")], by="sample")

  # combined countplots
  ggplot(countplot, aes(x=group, y=count, color=sample)) +
      #scale_y_log10() +  
      geom_beeswarm(size = 3,cex = 3) + ggtitle(paste0("Countplot: ", sym)) + theme(axis.text.x = element_text(angle=45, hjust=1))
  ggsave(filename = file.path(out.path, paste('countplot_', sym, '_' , samples_text, '_lsTPM', '.png', sep = '')))
}

# Rnf213 WT vs KO-celline
samples_subgroup_WT_KOcellline = c("35__HoxPU1_Neu_v2_1_S35_R1_001",
"36__HoxPU1_Neu_v2_2_S36_R1_001",
"37__HoxPU1_Neu_v2_3_S37_R1_001",
"41__PU1_WT_1731_S41_R1_001",
"42__PU1_WT_1691_S42_R1_001",
"43__PU1_WT_1735_S43_R1_001",
"44__PU1_WT_1734_S44_R1_001")
save_combined_countplots_subset(tpm, "Rnf213", samples_subgroup_WT_KOcellline, "WT_KO-PU-cellline")

# Rnf213 WT vs KO
samples_subgroup_WT_KO = c("45__PU1_Neu_1712_S45_R1_001",
"46__PU1_Neu_1715_S46_R1_001",
"48__PU1_Neu_1718_S48_R1_001",
"41__PU1_WT_1731_S41_R1_001",
"42__PU1_WT_1691_S42_R1_001",
"43__PU1_WT_1735_S43_R1_001",
"44__PU1_WT_1734_S44_R1_001")
save_combined_countplots_subset(tpm, "Rnf213", samples_subgroup_WT_KO, "WT_KO-PU")
save_combined_countplots(tpm, "Rnf213")


# Ets2
save_combined_countplots(tpm, "Ets2")
# Ccl3
save_combined_countplots(tpm, "Ccl3")
# Junb
save_combined_countplots(tpm, "Junb")
# pu.1
save_combined_countplots(tpm, "Spi1")
# Cxcl3
save_combined_countplots(tpm, "Cxcl3")
# Csf1
save_combined_countplots(tpm, "Csf1")
# Trib1
save_combined_countplots(tpm, "Trib1")
# Spp1
save_combined_countplots(tpm, "Spp1")
# Gpr68
save_combined_countplots(tpm, "Gpr68")
# Thbs1
save_combined_countplots(tpm, "Thbs1")
# Cd63
save_combined_countplots(tpm, "Cd63")
# C3
save_combined_countplots(tpm, "C3")
# Slpi
save_combined_countplots(tpm, "Slpi")
# S100a4
save_combined_countplots(tpm, "S100a4")


# Venn diagrams 
# o   HoxWT (alt, denn wir haben auch nur die WT aus dem Paper) vs. HoxKO_v2
# o   HoxKO_v2 vs. Hox2KO
# =>
#   o   HoxWT down mit Hox2Ko down
#   o   Overlap von denen für GO BP und MF

# genes down in hox WT in comparison to hoxKO_v2
hoxWT_down <- joined_table %>% 
                          filter(`padj_Hox-WT_Hox-KO_v2` < 0.05 ) %>% 
                          filter(`log2FoldChange_Hox-WT_Hox-KO_v2` < 0 ) %>% 
                          dplyr::select("log2FoldChange_Hox-WT_Hox-KO_v2", "ensembl")
# genes down in HoxKO_v2 in comparison to Hox_2KOs
hoxKO_v2_up <- joined_table %>% 
                          filter(`padj_Hox-KO_v2_Hox-2KOs` < 0.05 ) %>% 
                          filter(`log2FoldChange_Hox-KO_v2_Hox-2KOs` > 0 ) %>% 
                          dplyr::select("log2FoldChange_Hox-KO_v2_Hox-2KOs", "ensembl")
library(ggvenn)


x <- list("Hox-WT vs HoxKO_v2 down" = hoxWT_down[["ensembl"]], 
          "Hox-KO_v2 vs Hox-2KOs up" = hoxKO_v2_up[["ensembl"]])
ggvenn(
    x, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4
    )
ggsave(filename = file.path(out.path, paste('Venn_', 'HoxWT-down_Hox-2KO-down', '.png', sep = '')))

# take intersect and do pathway analysis
intersect_ <- intersect(hoxWT_down[["ensembl"]], hoxKO_v2_up[["ensembl"]])

joined_dfs <- hoxWT_down %>% left_join(hoxKO_v2_up, by="ensembl") %>% filter(ensembl %in% intersect_) %>% dplyr::arrange(desc( `log2FoldChange_Hox-WT_Hox-KO_v2`))

# for the universe take relevant samples and compute rowsums, remove rowsums < 10 
universe = joined_table.2 %>% column_to_rownames("ensembl") %>% dplyr::select(c("35__HoxPU1_Neu_v2_1_S35_R1_001",
                                        "36__HoxPU1_Neu_v2_2_S36_R1_001",
                                        "37__HoxPU1_Neu_v2_3_S37_R1_001",
                                        "41__PU1_WT_1731_S41_R1_001",
                                        "42__PU1_WT_1691_S42_R1_001",
                                        "43__PU1_WT_1735_S43_R1_001",
                                        "44__PU1_WT_1734_S44_R1_001",
                                        "38__HoxPU1_Neu_creEts2_c3_KO1_1_S38_R1_001",
                                        "39__HoxPU1_Neu_creEts2_c3_KO1_2_S39_R1_001",
                                        "40__HoxPU1_Neu_creEts2_c3_KO1_3_S40_R1_001"
                                        ))
keep <- rowSums(universe) >= 3
universe <- universe[keep, ]
universe <- rownames(universe)
library(clusterProfiler)
ego <- enrichGO(gene          = joined_dfs[["ensembl"]],
                  universe      = universe,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "Bp", # CC: cellular compartment, MF: molecul. function, BP: biol. process
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(ego, 10)
ego <- as_tibble(ego)
write.xlsx(ego, file=file.path(out.path, paste("intersection_HoxWT-down_Hoix2KO-down_GO-enrichment.BP","genes.rowsums.gt.10.xlsx", sep=".")), row.names=F, overwrite=T)
