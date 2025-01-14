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
library(fgsea)
library(data.table)                 
library(ComplexHeatmap)

norm.path <- "/home/max/projects/NGS/neutrophiles/RASflow/output/neutrophiles_primary/trans/dea/countGroup"
dea.path <- "/home/max/projects/NGS/neutrophiles/RASflow/output/neutrophiles_primary/trans/dea/DEA/gene-level"
out.path <- "/home/max/projects/NGS/neutrophiles/RASflow/output/neutrophiles_primary/trans/dea/visualization"
dir.create(file.path(out.path), showWarnings = FALSE)

yaml.file <- yaml.load_file('configs/config_main.yaml')
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
meta.file <- "/home/max/projects/NGS/neutrophiles/RASflow/configs/metadata_primary.tsv"
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
design <- model.matrix(~0+group.all)
dds_all <- DESeqDataSetFromTximport(txi_all, colData = colData, design = design)
dds_all <- estimateSizeFactors(dds_all)
normalized_counts_deseq2 <- counts(dds_all, normalized=TRUE)
normalized_counts_deseq2 <- as_tibble(normalized_counts_deseq2, rownames="ensembl")



# pca of all samples
library(gridExtra)
library(grid)
library(sjmisc)
library(org.Mm.eg.db)

df_pca <- prcomp(joined_table %>% column_to_rownames("ensembl") %>% 
                  replace(is.na(.), 0) %>% rotate_df()) # make ensembl the rownames again, replkace all nas with 0 and finally transpose matrix
df_out <- as.data.frame(df_pca$x)
# head(df_out)
# add group
df_out$group <- (meta.data %>% as_tibble() %>% column_to_rownames("sample") )[rownames(df_out),"group"]
df_out$treat <- (meta.data %>% as_tibble() %>% column_to_rownames("sample") )[rownames(df_out),"subject"]
# plt
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
p<-ggplot(df_out,aes(x=PC1,
                    y=PC2,
                    color=group,
                    # shape=treat,
                    label=substring(rownames(df_out), 1, 10) )) + 
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
  group.anno <- meta.data %>% column_to_rownames("sample") %>% dplyr::select("group") #, "subject")

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
  group.anno <- meta.data %>% column_to_rownames("sample") %>% dplyr::select("group") #, "subject")

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

# old_expression_clusters <- read.xlsx("/mnt/c/Users/masprang/Desktop/Projects/Neutrophil-PU.1-project/mRNA_Hoxis_naiive_Ca_clusters_genelist_JF.xlsx", 1, startRow = 2)

# clust_genes <- old_expression_clusters[c("Cluster", "Name")]
# clust_genes <- clust_genes[clust_genes$Cluster != "none",]

ets2_expression_clusters <- read.xlsx("/mnt/c/Users/masprang/Desktop/Projects/Neutrophil-PU.1-project/WT_crETS2_Down_HA.xlsx", 1, startRow = 1)

clust_genes <- ets2_expression_clusters[["symbol"]]


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
exp.plot <- zscores %>% filter(Symbol %in% clust_genes) %>%
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
ggsave(filename = file.path(out.path, 'avg_expression_score_clust_Hoxis.png'))


# comparison to proteomics expression
diff_express_protein_prot <- read.xlsx("/mnt/c/Users/masprang/Desktop/Projects/Neutrophil-PU.1-project/DEPs_pathway_PU1_BM_neutros_HA.xlsx", 1)
diff_express_protein_prot <- as_tibble(diff_express_protein_prot)
diff_express_protein_prot <- diff_express_protein_prot %>% separate_rows(Gene_name)
#diff_express_protein_prot <- diff_express_protein_prot %>% mutate()


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
  # load(file.path(output.path, "countGroup", 'tx2gene.RData'))
  # # noVersion because ensemble is used 
  # files.noVersion <- file.path(quant.path, samples, "quant_noVersion.sf")
  # names(files.noVersion) <- samples
  # # load data fro salom, gene level
  # txi <- tximport(files.noVersion, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")

  # subject <- factor(subject.all[c(which(group.all == control), which(group.all == treat))])
  # group <- relevel(factor(group.all[c(which(group.all == control), which(group.all == treat))]), ref = control)
  # colData = data.frame(samples, group)
  # design <- model.matrix(~group)
  # dds <- DESeqDataSetFromTximport(txi_all, colData = colData, design = design)
  dds <- dds_all[,samples]
  # new design, as some were dropped
  new_design <- design[c(which(group.all == control), which(group.all == treat)),c(paste0("group.all", control), paste0("group.all", treat))]
  dds@design <- new_design
  # Filtering
  if (filter.need) {
    keep <- rowSums(counts(dds)) >= 20
    dds <- dds[keep,]
  }
  ## specify the control group
  dds$group.all <- factor(dds$group.all)
  dds$group.all <- relevel(factor(dds$group.all), ref = control)
  colnames(dds) <- dds@colData$samples.all
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
  # fgsea immune sets
  
  ranks_input <- as_tibble(res)  %>% 
    dplyr::select("symbol", "stat") %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(symbol) %>% 
    summarize(stat=mean(stat))
  ranks_input <- deframe(ranks_input)
  names(ranks_input) <- toupper(names(ranks_input))
  immune_sets <- gmtPathways("/home/max/projects/NGS/neutrophiles/genesets_immune_regulation_MM.gmt")
  fgseaRes <- fgsea(pathways = immune_sets, 
                    stats    = ranks_input, nperm=1000)

  fwrite(fgseaRes, file = file.path(out.path, paste('fgsea_results_immunesets_', control, '_', treat, '.tsv', sep = '')), sep = "\t", sep2=c("", " ", ""))
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
  group.anno <- meta.data %>% filter(group %in% c(control, treat)) %>% column_to_rownames("sample") %>% dplyr::select("group") #, "subject")
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
  to_join <- tibb.dea %>% dplyr::select("ensembl", "log2FoldChange", "padj", "stat")
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
  
  # save table to check for log_fold_changes 
  write.xlsx(tibb.dea %>% full_join(as_tibble(countdata, rownames="ensembl")  , by="ensembl"), file=file.path(out.path, paste("table_degs", control, treat, ".xlsx", sep='_')), row.names=F, overwrite=T)

  # plot fold changes in genomic space
  # takes quite long
  # library(apeglm)
  # resGR <- lfcShrink(dds, coef=c("groupKO.PU"), type="apeglm", format="GRanges")
}

ens.str <- joined_table$ensembl
sym.str <- convert.id2symbol(ens.str)
joined_table <- joined_table %>% mutate(symbol = sym.str) 
# write the master table
write.xlsx(joined_table, file=file.path(out.path, "master.table.xlsx"), row.names=F, overwrite=T)
### heatmap combi
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

# complex heatmap with topp 20 DEGs?
ens.str <- normalized_counts_deseq2$ensembl
sym.str <- convert.id2symbol(ens.str)
normalized_counts_deseq2 <- normalized_counts_deseq2 %>% mutate(symbol = sym.str) 
# input for WT vs WTCA
samples_WT_WTCa <- factor(samples.all[c(which(group.all == "Prim-WT-U"), which(group.all == "Prim-WT-Ca"))]) 
tpm_WT_WTCa <- normalized_counts_deseq2 %>% dplyr::select(all_of(samples_WT_WTCa), "ensembl", "symbol")
heatmap_input_WT_WTCa <- tpm_WT_WTCa %>% 
  filter(ensembl %in% combined_degs$ensembl)  %>%
  column_to_rownames("symbol") %>%
  dplyr::select(-ensembl) %>%
  as.matrix()

# input for WT vs KO
samples_WT_KO <- factor(samples.all[c(which(group.all == "Prim-WT-U0"), which(group.all == "Prim-KO-U0"))]) 
tpm_WT_KO <- normalized_counts_deseq2 %>% dplyr::select(all_of(samples_WT_KO),"ensembl", "symbol")
heatmap_input_WT_KO <- tpm_WT_KO %>% 
  filter(ensembl %in% combined_degs$ensembl)  %>%
  column_to_rownames("symbol") %>%
  dplyr::select(-ensembl) %>%
  as.matrix()

h_WT_WTCa <- Heatmap(log2(heatmap_input_WT_WTCa + 1), name = "WTu vs WTCa")
h_WT_KO <- Heatmap(log2(heatmap_input_WT_KO + 1), name = "WTu vs KOu")

png(file = file.path(out.path, 'comparison-heatmaps-top20.png'),width=3300, height=3600, res=300)
h_WT_WTCa + h_WT_KO
dev.off()


# complex heatmap with topp 100 DEGs?
# input for WT vs WTCA
samples_WT_WTCa <- factor(samples.all[c(which(group.all == "Prim-WT-U"), which(group.all == "Prim-WT-Ca"))]) 
tpm_WT_WTCa <- normalized_counts_deseq2 %>% dplyr::select(all_of(samples_WT_WTCa), "ensembl", "symbol")
heatmap_input_WT_WTCa <- tpm_WT_WTCa %>% 
  filter(ensembl %in% combined_degs_100$ensembl)  %>%
  column_to_rownames("symbol") %>%
  dplyr::select(-ensembl) %>%
  as.matrix()

# input for WT vs KO
samples_WT_KO <- factor(samples.all[c(which(group.all == "Prim-WT-U0"), which(group.all == "Prim-KO-U0"))]) 
tpm_WT_KO <- normalized_counts_deseq2 %>% dplyr::select(all_of(samples_WT_KO),"ensembl", "symbol")
heatmap_input_WT_KO <- tpm_WT_KO %>% 
  filter(ensembl %in% combined_degs_100$ensembl)  %>%
  column_to_rownames("symbol") %>%
  dplyr::select(-ensembl) %>%
  as.matrix()

h_WT_WTCa <- Heatmap(log2(heatmap_input_WT_WTCa + 1), name = "WTu vs WTCa")
h_WT_KO <- Heatmap(log2(heatmap_input_WT_KO + 1), name = "WTu vs KOu")

png(file = file.path(out.path, 'comparison-heatmaps-top100.png'),width=3300, height=6600, res=300)
h_WT_WTCa + h_WT_KO
dev.off()
