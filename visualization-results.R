library(tidyverse)
library(yaml)
library(DESeq2)
library(biomaRt)
library(tximport)
library(mygene)
library(hash)

norm.path <- "/home/max/projects/NGS/neutrophiles/RASflow/output/neutrophiles/trans/dea/countGroup"
dea.path <- "/home/max/projects/NGS/neutrophiles/RASflow/output/neutrophiles/trans/dea/DEA/gene-level"
out.path <- "/home/max/projects/NGS/neutrophiles/RASflow/output/neutrophiles/trans/dea/visualization"

yaml.file <- yaml.load_file('configs/config_main_neutrophiles.yaml')


# extract the information from the yaml file
project <- yaml.file$PROJECT  # project name
dea.tool <- yaml.file$DEATOOL  # tool used for DEA
quant.path <- file.path(yaml.file$FINALOUTPUT, project, "trans/quant")
gene.level <- yaml.file$GENE_LEVEL  # whether to do gene-level analysis
controls <- yaml.file$CONTROL  # all groups used as control
treats <- yaml.file$TREAT  # all groups used as treat, should correspond to control
filter.need <- yaml.file$FILTER$yesOrNo
pair.test <- yaml.file$PAIR
meta.file <- "/home/max/projects/NGS/neutrophiles/RASflow/configs/metadata_copy.tsv"
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
samples.all <- meta.data$sample
group.all <- meta.data$group
subject.all <- meta.data$subject
# get complete counts
samples <- factor(samples.all) 
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
  txi_all <- tximport(files.noVersion, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")


all.dea.tibble <- tibble()
joined_table <- as.data.frame(txi_all$counts)
joined_table <- tibble(rownames_to_column(joined_table, "ensembl"))   
# filter rows, with only counts below 2
joined_table <- joined_table %>% filter_at(vars(-ensembl), any_vars(. > 2))

## create the DESeqDataSet
# for subset samples (comparisons)
for (i in 1:length(controls)) {

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

  png(file = file.path(out.path, paste('nornalization-comparison', control, '_', treat, '.png', sep = '')), title = 'comparison of transformations for counts')
  ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)  
  dev.off()


  # distances 
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- vsd@colData@rownames

  library("pheatmap")
  library("RColorBrewer")

  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

  png(file = file.path(out.path, paste('sample-distances', control, '_', treat, '.png', sep = '')), width=750, height=500,  title = 'Sample distances of variance stabilizing transformed counts')
  pheatmap(sampleDistMatrix,
          labels_row = rownames(sampleDistMatrix),
          show_colnames = F,
          clustering_distance_rows = sampleDists,
          clustering_distance_cols = sampleDists,
          col = colors)
  dev.off()

  # pca
  pcaData <- plotPCA(vsd, intgroup = "group", returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size =3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle("PCA with VST data")
  ggsave(filename = file.path(out.path, paste('PCA_', control, '_', treat, '.png', sep = '')))


  ### heatmap of topgenes
  library("genefilter")
  library("AnnotationDbi")
  library("org.Mm.eg.db")

  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
  mat  <- assay(vsd)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  # annotate mat
  ens.str <- rownames(mat)
  sym.str <- mapIds(org.Mm.eg.db,
                    keys=ens.str,
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multivals="first")
  sym.str[is.na(sym.str)] <- names(sym.str[is.na(sym.str)])  
  rownames(mat) <- sym.str
  group.anno <- as.data.frame(colData(vsd)[, c("samples", "group")])
  group.anno <- group.anno["group"]

  library(grid)
  ## Edit body of pheatmap:::draw_colnames, customizing it to your liking
  draw_colnames_45 <- function (coln, ...) {
    m = length(coln)
    x = (1:m)/m - 1/2/m
    grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
        hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
  }

  ## 'Overwrite' default draw_colnames with your own version 
  assignInNamespace(x="draw_colnames", value="draw_colnames_45",
  ns=asNamespace("pheatmap"))

  png(file = file.path(out.path, paste('top-gene-heatmap_', control, '_', treat, '.png', sep = '')),width=750, height=750,  title = 'PCA of variance stabilizing transformed counts')
  pheatmap(mat, annotation_col = group.anno) 
  dev.off()

  ## perform DEA
  dds <- DESeq(dds)
  res <- results(dds)

  # enrichment/ annotation
  topGene <- rownames(res)[which.min(res$padj)]

  png(file = file.path(out.path, paste('countplot_topgene_', control, '_', treat, '.png', sep = '')), title = 'PCA of variance stabilizing transformed counts')
  plotCounts(dds, gene = topGene, intgroup=c("group"))
  dev.off()

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


  library("ggbeeswarm")
  # count plot pu.1
  pu.1 <- rownames(res[res$symbol == "Spi1",] )
  plot <- plotCounts(dds, gene = pu.1 , intgroup=c("group", "samples"), returnData=T)
  ggplot(plot, aes(x=group, y=count, color=samples)) +
    scale_y_log10() +  geom_beeswarm(size = 3, cex = 3)
  ggsave(filename = file.path(out.path, paste('countplot_PU.1_', control, '_', treat, '.png', sep = '')))
  

  # count plot ets2
  ets2 <- rownames(res[res$symbol == "Ets2",] )
  plot <- plotCounts(dds, gene = ets2, intgroup=c("group", "samples"), returnData=T)
  ggplot(plot, aes(x=group, y=count, color=samples)) +
    scale_y_log10() +  geom_beeswarm(size = 3,cex = 3)
  
  ggsave(filename = file.path(out.path, paste('countplot_Ets2_', control, '_', treat, '.png', sep = '')))

  # count plot junb
  Junb <- rownames(res[res$symbol == "Junb",] )
  plot <- plotCounts(dds, gene = Junb, intgroup=c("group", "samples"), returnData=T)
  ggplot(plot, aes(x=group, y=count, color=samples)) +
    scale_y_log10() +  geom_beeswarm(size = 3,cex = 3)
  
  ggsave(filename = file.path(out.path, paste('countplot_Junb_', control, '_', treat, '.png', sep = '')))

  # count ccl3
  Ccl3 <- rownames(res[res$symbol == "Ccl3",] )
  plot <- plotCounts(dds, gene = Ccl3, intgroup=c("group", "samples"), returnData=T)
  ggplot(plot, aes(x=group, y=count, color=samples)) +
    scale_y_log10() +  geom_beeswarm(size = 3,cex = 3)
  ggsave(filename = file.path(out.path, paste('countplot_Ccl3_', control, '_', treat, '.png', sep = '')))
  

  # ordered after pval
  resOrdered <- res[order(res$log2FoldChange),]
  head(resOrdered)
  
  
  tibb.dea <- as_tibble(resOrdered , rownames="ensembl")
  tibb.dea <- tibb.dea[!is.na(tibb.dea$padj),]
  tibb.dea <- tibb.dea[tibb.dea$padj < 0.05,]
  # tibb.dea <- tibb.dea[abs(tibb.dea$log2FoldChange) > 2,]
  to_join <- tibb.dea %>% dplyr::select("ensembl", "log2FoldChange", "padj")
  to_join <- to_join %>% mutate("lfc>2" = as.numeric(abs(log2FoldChange)>2))
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

write.table(joined_table, file=file.path(out.path, "master.table.csv"), sep=",", row.names=F)
