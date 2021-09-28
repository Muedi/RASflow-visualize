## pathway analysis (kegg and go)
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
meta.file <- "/home/max/projects/NGS/neutrophiles/RASflow/configs/metadata_47removed.tsv"
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

  control <- controls[5] # KO-PU-Cellline
  treat <- treats[5] # KO-PU+KO-Ets2
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
  
  dds <- DESeq(dds)
  res <- results(dds, alpha=0.05)
  res <- res[order(res$pvalue),]
  summary(res)
  #packages for annotation
  library("genefilter")
  library("AnnotationDbi")
  library("org.Mm.eg.db")

  # use mapids 
  res$symbol = mapIds(org.Mm.eg.db,
                     keys=row.names(res), 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez = mapIds(org.Mm.eg.db,
                     keys=row.names(res), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$name =   mapIds(org.Mm.eg.db,
                     keys=row.names(res), 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

# packages for pathways:
# library(pathview)
# library(gage)
# library(gageData)
# data(kegg.sets.mm)
# data(sigmet.idx.mm)
# kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
# head(kegg.sets.mm, 3)
# # gage required fold changes vector named with entrez ids.
# foldchanges = res$log2FoldChange
# names(foldchanges) = res$entrez
# head(foldchanges)
# # the actual pathway analysis
# # same.dir =TRUE gives us pathways in both up- and down-regulated directions
# keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)
# # Get the upregulated (greater) pathways
# keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
#   tibble::as_tibble() %>% 
#   filter(row_number()<=5) %>% 
#   .$id %>% 
#   as.character()

# # Get the IDs.
# keggresids = substr(keggrespathways, start=1, stop=8)
# # plot pathways as graphs
# # Define plotting function for applying later
# #plot_pathway = function(pid) pathview(gene.data=foldchanges, kegg.dir=paste(out.path,"/pathways", sep=''), pathway.id=pid, species="mmu", new.signature=FALSE)

# dir.create(file.path(paste(out.path,"/pathways", sep='')), showWarnings = FALSE)
# normal_wd <- getwd()
# setwd(file.path(paste(out.path,"/pathways", sep='')))
# # plot multiple pathways (plots saved to disk and returns a throwaway list object)
# tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges,  pathway.id=pid, species="mmu"))
# setwd(normal_wd)

# ### GO gene ontology enrichemnt
# data(go.sets.mm)
# data(go.subs.mm)
# gobpsets = go.sets.mm[go.subs.mm$BP]
# gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
# lapply(gobpres, head)

# pathway analysis with clusterprofiler
library(clusterProfiler)
# similar to gage requires fold changes vector named with ids.
geneList = res$log2FoldChange
names(geneList) = rownames(res)
# but it additionally needs to be sorted
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)

gene <- names(geneList) #[abs(geneList) > 2] )
degs <- names(geneList[abs(geneList) > 2] )
# gene.df <- bitr(gene, fromType = "ENSEMBL",
#         toType = c("ENTREZID", "SYMBOL"),
#         OrgDb = org.Mm.eg.db)
# head(gene.df)

ego <- enrichGO(gene          = degs,
                universe      = gene,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "MF", # CC: cellular compartment, MF: molecul. function, BP: biol. process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego, 10)


ggo <- groupGO(gene     = degs,
               OrgDb    = org.Mm.eg.db,
                keyType       = 'ENSEMBL',

               ont      = "MF",
               level    = 4,
               readable = TRUE)
head(ggo, 10)




# add the proteomics data of the old project. 
# check how data is comparable to the new knockout (cellline and double KO)
# read protein data
library(openxlsx)
old_expression_clusters <- read.xlsx("/mnt/c/Users/masprang/Desktop/Projects/Neutrophil-PU.1-project/mRNA_Hoxis_naiive_Ca_clusters_genelist_JF.xlsx", 1, startRow = 2)
diff_express_protein_old <- read.xlsx("/mnt/c/Users/masprang/Desktop/Projects/Neutrophil-PU.1-project/DEPs_pathway_PU1_BM_neutros_HA.xlsx", 1)
pathways_old <- read.xlsx("/mnt/c/Users/masprang/Desktop/Projects/Neutrophil-PU.1-project/DEPs_pathway_PU1_BM_neutros_HA.xlsx", 2)

# gsea for Clusters
library(fgsea)
# get list of clusters and genes, to use as pathways for fgsea
clust_genes <- old_expression_clusters[c("Cluster", "Name")]
clust_genes$entrez = mapIds(org.Mm.eg.db,
                     keys=clust_genes$Name, 
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")

clusters <-unique(clust_genes[["Cluster"]]) 
clust_pathways_symbols <- list(
  clust_genes[clust_genes$Cluster == clusters[1],"Name"],
  clust_genes[clust_genes$Cluster == clusters[2],"Name"],
  clust_genes[clust_genes$Cluster == clusters[3],"Name"],
  clust_genes[clust_genes$Cluster == clusters[4],"Name"],
  clust_genes[clust_genes$Cluster == clusters[5],"Name"]
)
clust_genes <- clust_genes[!is.na(clust_genes$entrez),]
clust_pathways_entrez <- list(
  clust_genes[clust_genes$Cluster == clusters[1],"entrez"],
  clust_genes[clust_genes$Cluster == clusters[2],"entrez"],
  clust_genes[clust_genes$Cluster == clusters[3],"entrez"],
  clust_genes[clust_genes$Cluster == clusters[4],"entrez"],
  clust_genes[clust_genes$Cluster == clusters[5],"entrez"]
)


names(clust_pathways) <- clusters

