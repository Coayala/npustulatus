# ---- Packages needed ----
library(edgeR)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(ggsci)
library(pheatmap)
library(clusterProfiler)
library(topGO)
library(tidyverse)
# ---- Count load and normalization ----

# Create a variable storing the paths of kallisto counts

setwd("C:/Users/Christian/Desktop/Christian/Hagen_projects/beetle")

samples <- read_tsv("Nicropus_inputs/samples_switched.txt")

counts_paths <- file.path('kallisto_Npus-counts', samples$sample_name)
kallisto_counts <- catchKallisto(counts_paths)
colnames(kallisto_counts$counts) <- samples$sample_name

# create a grouping vector based on file information

group <- factor(samples$Condition)
group <- relevel(group, ref = 'Starved')

# Create a DGEList object with the counts and grouping information

dge <- DGEList(counts = kallisto_counts$counts/kallisto_counts$annotation$Overdispersion, group = group)

# Filter lowly expressed genes and calculate normalization factors

isexpr <- filterByExpr(dge, group = group)
dge <- dge[isexpr, keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# Obtain normalized counts

scale <-  dge$samples$lib.size*dge$samples$norm.factors
norm.counts <- round(t(t(kallisto_counts$counts)/scale)*mean(scale))

# Transform normalized counts into a tibble to write the table

norm.counts_table <- as_tibble(norm.counts) %>% 
  add_column(target_id = rownames(norm.counts), .before = samples$sample_name[1])

write_tsv(norm.counts_table, 'Nicropus_switched_outputs/normalized_counts.tsv')

# ---- Sample QC ----

# Explore data with MDS plots

tiff('Nicropus_switched_outputs/MDS_plot_kallisto.tiff', width = 900)

my_colors <- rep(c('black', 'red'), 3)
plotMDS(dge, col = my_colors[group], 
        main = "MDS plot")
legend('top', legend = levels(group), col = my_colors, pch = 16)

dev.off()

# Log2 transform read counts

log.nc <- log2(norm.counts + 1)

# PCA data

pca_data <- prcomp(t(log.nc))

pca_data4plot <- data.frame(PC1 = pca_data$x[,1], 
                            PC2 = pca_data$x[,2],
                            sample = colnames(log.nc), 
                            condition = rep(c("Fed", "Starved"), each =3),
                            sex = rep(c("Female", "Male"), each = 6))

pca_plot <- ggplot(pca_data4plot,
                   aes(x = PC1,
                       y = PC2,
                       color = condition,
                       shape = sex)) +
  geom_point(size = 8)

pca_plot

# Explore read counts

## Box plot of transformed and unstransformed read counts

boxplot1 <- ggplot(stack(as.data.frame(log.nc)), 
                   aes(x=ind, y = values)) +
  theme_bw() +
  geom_boxplot(notch = TRUE) +
  labs (x = 'Samples', y = 'log2(read counts', title = "Log2-transformed read counts" ) +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 6))

boxplot2 <- ggplot(stack(as.data.frame(norm.counts)), 
                   aes(x=ind, y = values)) +
  theme_bw() +
  geom_boxplot(notch = TRUE) +
  labs (x = 'Samples', y = 'read counts', title = "Untransformed read counts" ) +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 6))

boxplot.log2_vs_unt <- plot_grid(boxplot1, boxplot2, labels = 'AUTO', rel_widths = c(1,1))
ggsave('Nicropus_switched_outputs/box_vs_box_kallisto.jpg', boxplot.log2_vs_unt)

## MA plots

avalues <- (log.nc[,1] + log.nc[,2])/2
mvalues <- (log.nc[,1] - log.nc[,2])

MA_data <- tibble(avalues, mvalues)

MA1 <- ggplot(MA_data) +
  geom_point(aes(x = avalues, y = mvalues)) +
  labs(x = 'A', y = 'M')

## Pearson Correlation heatmap

sampleCor <- (1- cor(log.nc, method = 'pearson'))

mapcolor <- colorRampPalette(rev(brewer.pal(9, 'GnBu')))(16)

col_annot <- column_to_rownames(samples, var = 'sample_name')
side_colors <- list(
  Sex = c(Male = '#5977ff', Female = '#f74747'),
  Condition = c( Fed = '#82ed82', Starved = '#e89829')
)

pdf('Nicropus_switched_outputs/pearson_heatmap_kallisto.pdf')

pheatmap(sampleCor,
         annotation_row = col_annot,
         annotation_col = col_annot,
         color = mapcolor,
         annotation_colors = side_colors,
         main = 'Pearson correlation heatmap (N. pustulatus)')
dev.off()


### Distance heatmap

sampleDists <- as.matrix(dist(t(log.nc)))
sampleDists <- sampleDists/max(sampleDists)

pdf('Nicropus_switched_outputs/distance_heatmap_kallisto.pdf')

pheatmap(sampleDists,
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         cutree_cols = 5,
         annotation_row = col_annot,
         annotation_col = col_annot,
         color = mapcolor,
         annotation_colors = side_colors,
         main = 'Sample distance (N. pustulatus)')
dev.off()

# CPM Heatmap

log.cpm <- cpm(dge, log = TRUE, prior.count = 1)

pheatmap(log.cpm,
         show_rownames = FALSE,
         annotation_col = col_annot,
         cutree_cols = 4)

# ---- Differential expression analysis ----

# Perform differential expression analysis with edgeR

## Create a matrix design based on sample information

design <- model.matrix(~group, dge$samples)

## Estimate common and tagwise dispersion

dge <- estimateDisp(dge, design = design, robust = TRUE)

plotBCV(dge)

dgeTest <- exactTest(dge)
dgeTest

## Inspect results of DGE

summary(de <- decideTests(dgeTest, adjust.method = 'fdr', p.value = 0.05))

top_results <- topTags(dgeTest, n=nrow(dgeTest$table))
top_sig_results <- topTags(dgeTest, n=nrow(dgeTest$table), p.value = 0.05)

head(top_results)
sum(top_results$table$FDR < 0.05)

# ---- Results annotation ----

## Load the results of Blast search of transcriptome to annotate the results and curate the gene names and accesion numbers

blast_table <- read_tsv('Nicropus_inputs/Nicropus_results.blastx', col_names = FALSE)
blast_table <- blast_table %>% 
  separate(X2, c('delete', 'refseq_acc', 'delete2'), 
           sep = "[|]") %>%
  select(-c(delete, delete2)) %>% 
  separate(X6, c('unite', 'name'), 
           sep = 'PREDICTED:') %>% 
  unite('name', unite:name, sep = '') %>% 
  separate(name, c('name', 'delete'), sep = '\\[') %>% 
  select(-delete)


colnames(blast_table) <- c('transcript_id', 'refseq_acc', 'p.ident', 'Evalue', 'bitscore', 'name', 'taxid')
blast_table$name <- str_replace(blast_table$name, 'uncharacterized protein', 'u.c')
blast_table$name <- str_remove(blast_table$name, 'LOW QUALITY PROTEIN:')
blast_table$name <- str_trim(blast_table$name, side = 'left')

taxid_names <- read_tsv('Nicropus_inputs/tax_report (1).txt')
taxid_names <- distinct(taxid_names, taxid, taxname)
blast_table <- left_join(blast_table, taxid_names, by = 'taxid' )

final_table <- blast_table

## Eliminate the transcripts that match twice to the same refseq gene

final_table <- distinct(final_table, transcript_id, refseq_acc, .keep_all = TRUE)

## merge Blast annotations with Swissprot annotations and GO terms

### Read the GO annotations obtained through BLAST of the Uniprot_Sprot database

uniprot_table <- read_tsv('Nicropus_inputs/Nicropus.t-id.uni-acc.txt', col_names = FALSE)
uniprot_table <- uniprot_table %>% 
  separate(X2, c('delete', 'uniprot_acc', 'delete2'),
           sep = "[|]") %>% 
  select(-delete, -delete2)
colnames(uniprot_table) <- c('transcript_id', 'uniprot_acc')

uniprot_accesions <- tibble(uniprot_table$uniprot_acc)
uniprot_accesions <- distinct(uniprot_accesions, uniprot_table$uniprot_acc)
write_tsv(uniprot_accesions, 'Nicropus_inputs/Nicropus.uniq_uniprot_acc.txt', col_names = FALSE)

GO_ids <- read_tsv('Nicropus_inputs/Nicropus.uni_acc.go_id.txt')

uniprot_table <- left_join(uniprot_table, GO_ids, by = c('uniprot_acc' = 'Entry'))

### Read KEGG annotation obtained form BLAST of the Uniprot_Sprot database

KEGG_ids <- read_tsv('Nicropus_ko_insects.terms')

uniprot_table <- left_join(uniprot_table, KEGG_ids, by = 'transcript_id')
  
uniprot_table <- left_join(uniprot_table, path_final, by = c('KO' = 'koID'))

uniprot_table <- distinct(uniprot_table, transcript_id, uniprot_acc, .keep_all = TRUE)

final_table <- left_join(final_table, uniprot_table, by = 'transcript_id')

### write the annotation table

write_tsv(final_table, 'Nicropus_switched_outputs/Nicropus.transcriptome_annotation.txt')

## Annotate results with BLAST refseq IDS, Swissprot IDS and GO IDS

top_results_annotated <- top_results$table
top_results_annotated <- rownames_to_column(top_results$table, var = 'transcript_id')
top_results_annotated <- left_join(top_results_annotated, final_table, by = 'transcript_id')
top_results_annotated$`Gene ontology IDs` <- str_replace(top_results_annotated$`Gene ontology IDs`, ';', ',')
write_tsv(top_results_annotated, 'Nicropus_switched_outputs/Nicropus.all_results.annotated.tsv')

top_sig_results_annotated <- top_sig_results$table
top_sig_results_annotated <- rownames_to_column(top_sig_results$table, var = 'transcript_id')
top_sig_results_annotated <- left_join(top_sig_results_annotated, final_table, by = 'transcript_id')
top_sig_results_annotated$`Gene ontology IDs` <- str_replace(top_sig_results_annotated$`Gene ontology IDs`, ';', ',')
write_tsv(top_sig_results_annotated, 'Nicropus_switched_outputs/Nicropus.sig_DE.annotated.tsv')


# ---- Plotting of annotated results ----

## MA and MD plots from the edgeR package

pdf('Nicropus_switched_outputs/MA_plot.pdf')

plotSmear(dgeTest, 
          de.tags = rownames(top_sig_results$table),
          main = 'MA plot (edgeR)')
top4plot <- head(top_sig_results_annotated, 50)

text(top4plot$logCPM, top4plot$logFC, 
     labels = top4plot$name,
     col = 'blue',
     cex = 0.75)

dev.off()

pdf('Nicropus_switched_outputs/MD_plot.pdf')

plotMD(dgeTest)
abline(h = c(-2,2), col = 'blue')

dev.off()


## Volcano plot

## pvalue and lfc cutoffs

lfc <- 2
pval <- 0.05

### results to annotate

top4plot <- top_sig_results_annotated[abs(top_results_annotated$logFC) > 1.25*lfc & 
                                -log10(top_results_annotated$PValue) > 3*-log10(pval),]
top4plot <- top4plot %>% 
  separate(`Protein names`, c('prot_names', 'delete'),
           sep = '[(]') %>% 
  select(-delete)

volcanoplot <- ggplot(top_results_annotated, 
                      aes(x = logFC, 
                          y = -log10(PValue))) +
  theme_bw() +
  geom_point(color = ifelse(abs(top_results_annotated$logFC) > lfc & 
                              -log10(top_results_annotated$PValue) > -log10(pval), "#FF0000", "#000000")) +
  geom_vline(xintercept = c(-lfc, lfc),
             linetype = 'dotted',
             size = 1,
             color = 'blue') +
  geom_hline(yintercept = -log10(pval),
             linetype = 'dotted',
             size = 1,
             color = 'blue') +
  labs(title = 'Volcano plot', 
       x = expression("Log"[2]*" Fold Change"), 
       y = expression("-Log"[10]*" pvalue")) +
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = 'bold', 
                                  size = 18)) +
  geom_text_repel(data = top4plot,
                  aes(x = logFC,
                      y= -log10(PValue)),
                  label = top4plot$prot_names,
                  size = 2)

volcanoplot

ggsave('Nicropus_switched_outputs/volcano_plot_kallisto.png', volcanoplot)

# ---- Heatmap and table of most significant genes ----

## Heatmap of the most significant genes

### Obtain the log-counts-per-million log(cpm)

log.cpm <- cpm(dge, log = TRUE, prior.count = 1)

### Select 1% of differentially expressed genes to produce a heatmap

selected_log.cpm <- log.cpm[rownames(top_results$table)[top_results$table$FDR<0.05 & abs(top_results$table$logFC)>2],]

names4heatmap <- rownames(selected_log.cpm)

names4heatmap <- as_tibble(names4heatmap)

names4heatmap <- left_join(names4heatmap, final_table, by  = c('value' = 'transcript_id'))
names4heatmap <- names4heatmap %>% 
  separate(`Protein names`, c('prot_names', 'delete'),
           sep = '[(]') %>% 
  select(-delete)

for(i in 1:nrow(names4heatmap)){
  if (is.na(names4heatmap$prot_names[i])){
    names4heatmap$prot_names[i] <- names4heatmap$value[i]
  }
}

rownames(selected_log.cpm) <- names4heatmap$prot_names



### Draw the heatmap

mapcolor <- colorRampPalette(c('green', 'green2', 'green4', 'darkgreen', 'black', 'red'))(100)[100:1]

tiff('Nicropus_switched_outputs/DGE_heatmap_kallisto.tiff', height = 2000, width = 1500)

setEPS()
postscript('Nicropus_switched_outputs/DGE_heatmap.eps', height = 7.5, width = 7)

pheatmap(selected_log.cpm,
         cutree_cols = 2,
         cutree_rows = 4,
         annotation_col = col_annot,
         color = mapcolor,
         annotation_colors = side_colors,
         main = 'Differentially expressed genes',
         #fontsize = 20,
         fontsize_row = 4,
         )

dev.off()

## Extract and sort downregulated and upregulated expressed genes

sigDownReg <- top_results_annotated[top_results_annotated$FDR< 0.05 & top_results_annotated$logFC < 0,]

sigUpReg <- top_results_annotated[top_results_annotated$FDR< 0.05 & top_results_annotated$logFC > 0,]


## Writing results in csv files

write_tsv(sigDownReg, 'Nicropus_switched_outputs/downreg_genes.tsv')
write_tsv(sigUpReg, 'Nicropus_switched_outputs/upreg_genes.tsv')

# ---- GO enrichment analysis (All DE Genes) ----

### Save a table with Uniprot/Sprot gene accesion numbers and their respective GO annotations

acc_N_GO <- top_sig_results_annotated %>% 
  select(uniprot_acc, `Gene ontology IDs`) %>% 
  drop_na()

write_tsv(acc_N_GO, 'Nicropus_inputs/acc_N_GO_.tsv', col_names = FALSE)

## Create a gene2GO mapping

geneID2GO <- readMappings(file = 'Nicropus_inputs/acc_N_GO.tsv' )

## Create a named vector which contain the information of all genes

geneNames <- names(geneID2GO)
geneListTemp <- drop_na(top_sig_results_annotated, uniprot_acc)
geneList <- geneListTemp$FDR
names(geneList) <- geneListTemp$uniprot_acc

## Create a function to retrieve the differentially expressed genes based in a determined criteria 
##  (function is predefined to find genes with Pvalue < 0.01 in this case)

topDiffGenes <- function(allScore) {
  return(allScore < 0.01) 
}

##create the topGO object (must be done for each ontology)

#### Biological Process

GOdataBP <- new("topGOdata",
              description = "GO analysis",
              ontology = "BP", 
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO,
              geneSel = topDiffGenes)

### Run enrichment test using runTest functions, with different statistics tests

#### Fisher test
resultFisBP <- runTest(GOdataBP, algorithm = "classic", statistic = "fisher")
resultweightBP <- runTest(GOdataBP, algorithm = 'weight', statistic = 'fisher')
resultKSBP <- runTest(GOdataBP, algorithm = 'classic', statistic = 'ks')
resultweight2BP <- runTest(GOdataBP, algorithm = 'weight01', statistic = 'fisher')


geneData(resultweight2BP)
show(resultweight2BP)

allResGOBP <- GenTable(GOdataBP, Fis = resultFisBP, KS = resultKSBP, weight = resultweightBP, 
                                     ranksOf = 'Fis', topNodes = 20)

test.stat <- new("weight01", testStatistic = GOFisherTest, name = "Fisher test", cutOff = 0.05)
resultFisher <- getSigGroups(GOdataBP, test.stat)

write_tsv(allResGOBP, 'Nicropus_switched_outputs/all_GOBP.tsv', col_names = TRUE)

#### Molecular Function

GOdataMF <- new("topGOdata",
                description = "GO analysis",
                ontology = "MF", 
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO,
                geneSel = topDiffGenes)

### RUn enrichment test using runTest functions, with different statistics tests

#### Fisher test
resultFisMF <- runTest(GOdataMF, algorithm = "classic", statistic = "fisher")
resultweightMF <- runTest(GOdataMF, algorithm = 'weight', statistic = 'fisher')
resultKSMF <- runTest(GOdataMF, algorithm = 'classic', statistic = 'ks')


geneData(resultFisMF)

allResGOMF <- GenTable(GOdataMF, Fis = resultFisMF, KS = resultKSMF, weight = resultweightMF, 
                       ranksOf = 'Fis', topNodes = 20)

write_tsv(allResGOMF, 'Nicropus_switched_outputs/all_GOMF.tsv', col_names = TRUE)

#### Celullar component

GOdataCC <- new("topGOdata",
                description = "GO analysis",
                ontology = "CC", 
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO,
                geneSel = topDiffGenes)

### RUn enrichment test using runTest functions, with different statistics tests

#### Fisher test
resultFisCC <- runTest(GOdataCC, algorithm = "classic", statistic = "fisher")
resultweightCC <- runTest(GOdataCC, algorithm = 'weight', statistic = 'fisher')
resultKSCC <- runTest(GOdataCC, algorithm = 'classic', statistic = 'ks')


geneData(resultFisCC)

allResGOCC <- GenTable(GOdataCC, Fis = resultFisCC, KS = resultKSCC, weight = resultweightCC, 
                       ranksOf = 'Fis', topNodes = 20)

write_tsv(allResGOCC, 'Nicropus_switched_outputs/all_GOCC.tsv', col_names = TRUE)

###### Plot GO enrichment results

allResGOBP <- allResGOBP %>% 
  mutate(Domain = 'Biological Process')

allResGOMF <- allResGOMF %>% 
  mutate(Domain = 'Molecular Function')

allResGOCC <- allResGOCC %>% 
  mutate(Domain = 'Cellular component')

allGO4plot <- rbind(allResGOBP, allResGOMF, allResGOCC) %>% 
  #slice_max(order_by = Significant, n = 30, with_ties = FALSE) %>%
  arrange(Domain, Significant) %>% 
  mutate(Position = n():1) 

GOgraph <- ggplot(allGO4plot, aes(x = fct_reorder(Term, dplyr::desc(Position)), y = Significant, fill = Domain)) +
  #plot the number of times each term occur, normalized
  geom_col() +
  coord_flip() +
  theme_bw() +
  #Aesthetic modifications
  theme(axis.text.x = element_text(hjust = 1, color = "black"), axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 6), legend.title = element_text(size = 8),
        legend.key.size = unit(0.2, "in"), plot.title = element_text(size = 11, hjust = 0.5),
        plot.margin = margin(10, 50, 10, 50)) +
  #Change name of the legend and legend labels
  scale_fill_jama(name = "Domain") +
  #add a plot title in bolds and italic
  labs(x = NULL, y = "Number of genes (significant)", title = "Enriched GO terms") +
  theme(plot.title = element_text(face = "bold"))
GOgraph

ggsave("Nicropus_switched_outputs/enriched_GO.png", GOgraph, height = 10)


# ---- GO enrichment analysis (Upregulated Genes) ----

### Save a table with Uniprot/Sprot gene accesion numbers and their respective GO annotations

UP_acc_N_GO <- sigUpReg %>% 
  select(uniprot_acc, `Gene ontology IDs`) %>% 
  drop_na()

write_tsv(UP_acc_N_GO, 'Nicropus_inputs/UP_acc_N_GO.tsv', col_names = FALSE)

## Create a gene2GO mapping

#geneID2GO <- readMappings(file = 'Nicropus_inputs/UP_acc_N_GO.tsv' )

## Create a named vector which contain the information of all genes

geneNames <- names(geneID2GO)
geneListTemp <- drop_na(sigUpReg, uniprot_acc)
geneList <- geneListTemp$FDR
names(geneList) <- geneListTemp$uniprot_acc

##create the topGO object (must be done for each ontology)

#### Biological Process

GOdataBP <- new("topGOdata",
                description = "GO analysis",
                ontology = "BP", 
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO,
                geneSel = topDiffGenes)


### RUn enrichment test using runTest functions, with different statistics tests

#### Fisher test
resultFisBP <- runTest(GOdataBP, algorithm = "classic", statistic = "fisher")
resultweightBP <- runTest(GOdataBP, algorithm = 'weight', statistic = 'fisher')
resultKSBP <- runTest(GOdataBP, algorithm = 'classic', statistic = 'ks')


geneData(resultFisBP)

UPResGOBP <- GenTable(GOdataBP, Fis = resultFisBP, KS = resultKSBP, weight = resultweightBP, 
                       ranksOf = 'Fis', topNodes = 20)

write_tsv(UPResGOBP, 'Nicropus_switched_outputs/UP_GOCC.tsv', col_names = TRUE)

#### Molecular Function

GOdataMF <- new("topGOdata",
                description = "GO analysis",
                ontology = "MF", 
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO,
                geneSel = topDiffGenes)

### RUn enrichment test using runTest functions, with different statistics tests

#### Fisher test
resultFisMF <- runTest(GOdataMF, algorithm = "classic", statistic = "fisher")
resultweightMF <- runTest(GOdataMF, algorithm = 'weight', statistic = 'fisher')
resultKSMF <- runTest(GOdataMF, algorithm = 'classic', statistic = 'ks')


geneData(resultFisMF)

UPResGOMF <- GenTable(GOdataMF, Fis = resultFisMF, KS = resultKSMF, weight = resultweightMF, 
                       ranksOf = 'Fis', topNodes = 20)

write_tsv(UPResGOMF, 'Nicropus_switched_outputs/UP_GOMF.tsv', col_names = TRUE)

#### Celullar component

GOdataCC <- new("topGOdata",
                description = "GO analysis",
                ontology = "CC", 
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO,
                geneSel = topDiffGenes)

### RUn enrichment test using runTest functions, with different statistics tests

#### Fisher test
resultFisCC <- runTest(GOdataCC, algorithm = "classic", statistic = "fisher")
resultweightCC <- runTest(GOdataCC, algorithm = 'weight', statistic = 'fisher')
resultKSCC <- runTest(GOdataCC, algorithm = 'classic', statistic = 'ks')


geneData(resultFisCC)

UPResGOCC <- GenTable(GOdataCC, Fis = resultFisCC, KS = resultKSCC, weight = resultweightCC, 
                       ranksOf = 'Fis', topNodes = 20)

write_tsv(UPResGOCC, 'Nicropus_switched_outputs/UP_GOCC.tsv', col_names = TRUE)

###### Plot GO enrichment results

UPResGOBP <- UPResGOBP %>% 
  mutate(Domain = 'Biological Process')

UPResGOMF <- UPResGOMF %>% 
  mutate(Domain = 'Molecular Function')

UPResGOCC <- UPResGOCC %>% 
  mutate(Domain = 'Cellular component')

allGO4plot_UP <- rbind(UPResGOBP, UPResGOMF, UPResGOCC) %>% 
  slice_max(order_by = Significant, n = 30, with_ties = FALSE) %>%
  arrange(Domain, Significant) %>% 
  mutate(Position = n():1) %>% 
  mutate(Regulation = 'Upregulated genes')


UP_GOgraph <- ggplot(allGO4plot_UP, aes(x = fct_reorder(Term, dplyr::desc(Position)), y = Significant, fill = Domain)) +
  #plot the number of times each term occur, normalized
  geom_col() +
  coord_flip() +
  theme_bw() +
  #Aesthetic modifications
  theme(axis.text.x = element_text(hjust = 1, color = "black"), axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 6), legend.title = element_text(size = 8),
        legend.key.size = unit(0.2, "in"), plot.title = element_text(size = 11, hjust = 0.5),
        plot.margin = margin(10, 50, 10, 50)) +
  #Change name of the legend and legend labels
  scale_fill_jama() +
  #add a plot title in bolds and italic
  labs(x = NULL, y = "Number of genes (significant)", title = "Upregulated genes") +
  theme(plot.title = element_text(face = "bold"))
UP_GOgraph

ggsave("Nicropus_switched_outputs/UP_enriched_GO.png", UP_GOgraph, dpi = 320)

# ---- GO enrichment analysis (Downregulated Genes) ----

### Save a table with Uniprot/Sprot gene accesion numbers and their respective GO annotations

DOWN_acc_N_GO <- sigDownReg %>% 
  select(uniprot_acc, `Gene ontology IDs`) %>% 
  drop_na()

write_tsv(DOWN_acc_N_GO, 'Nicropus_inputs/DOWN_acc_N_GO.tsv', col_names = FALSE)

## Create a gene2GO mapping

#geneID2GO <- readMappings(file = 'Nicropus_inputs/DOWN_acc_N_GO.tsv' )

## Create a named vector which contain the information of all genes

geneNames <- names(geneID2GO)
geneListTemp <- drop_na(sigDownReg, uniprot_acc)
geneList <- geneListTemp$FDR
names(geneList) <- geneListTemp$uniprot_acc

##create the topGO object (must be done for each ontology)

#### Biological Process

GOdataBP <- new("topGOdata",
                description = "GO analysis",
                ontology = "BP", 
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO,
                geneSel = topDiffGenes)


### RUn enrichment test using runTest functions, with different statistics tests

#### Fisher test
resultFisBP <- runTest(GOdataBP, algorithm = "classic", statistic = "fisher")
resultweightBP <- runTest(GOdataBP, algorithm = 'weight', statistic = 'fisher')
resultKSBP <- runTest(GOdataBP, algorithm = 'classic', statistic = 'ks')

geneData(resultFisBP)

DOWNResGOBP <- GenTable(GOdataBP, Fis = resultFisBP, KS = resultKSBP, weight = resultweightBP, 
                       ranksOf = 'Fis', topNodes = 20)

write_tsv(DOWNResGOBP, 'Nicropus_switched_outputs/DOWN_GOBP.tsv', col_names = TRUE)

#### Molecular Function

GOdataMF <- new("topGOdata",
                description = "GO analysis",
                ontology = "MF", 
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO,
                geneSel = topDiffGenes)

### RUn enrichment test using runTest functions, with different statistics tests

#### Fisher test
resultFisMF <- runTest(GOdataMF, algorithm = "classic", statistic = "fisher")
resultweightMF <- runTest(GOdataMF, algorithm = 'weight', statistic = 'fisher')
resultKSMF <- runTest(GOdataMF, algorithm = 'classic', statistic = 'ks')

geneData(resultFisMF)

DOWNResGOMF <- GenTable(GOdataMF, Fis = resultFisMF, KS = resultKSMF, weight = resultweightMF, 
                       ranksOf = 'Fis', topNodes = 20)

write_tsv(DOWNResGOMF, 'Nicropus_switched_outputs/DOWN_GOMF.tsv', col_names = TRUE)

#### Celullar component

GOdataCC <- new("topGOdata",
                description = "GO analysis",
                ontology = "CC", 
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO,
                geneSel = topDiffGenes)

### RUn enrichment test using runTest functions, with different statistics tests

#### Fisher test
resultFisCC <- runTest(GOdataCC, algorithm = "classic", statistic = "fisher")
resultweightCC <- runTest(GOdataCC, algorithm = 'weight', statistic = 'fisher')
resultKSCC <- runTest(GOdataCC, algorithm = 'classic', statistic = 'ks')

geneData(resultFisCC)

DOWNResGOCC <- GenTable(GOdataCC, Fis = resultFisCC, KS = resultKSCC, weight = resultweightCC, 
                       ranksOf = 'Fis', topNodes = 20)

write_tsv(DOWNResGOCC, 'Nicropus_switched_outputs/DOWN_GOCC.tsv', col_names = TRUE)

###### Plot GO enrichment results

DOWNResGOBP <- DOWNResGOBP %>% 
  mutate(Domain = 'Biological Process')

DOWNResGOMF <- DOWNResGOMF %>% 
  mutate(Domain = 'Molecular Function')

DOWNResGOCC <- DOWNResGOCC %>% 
  mutate(Domain = 'Cellular component')

allGO4plot_DOWN <- rbind(DOWNResGOBP, DOWNResGOMF, DOWNResGOCC) %>% 
  # take the most represented GO terms
  slice_max(order_by = Significant, n = 30, with_ties = FALSE) %>%
  arrange(Domain, Significant) %>% 
  mutate(Position = n():1) %>% 
  mutate(Regulation = 'Downregulated genes')

DOWN_GOgraph <- ggplot(allGO4plot_DOWN, aes(x = fct_reorder(Term, dplyr::desc(Position)), y = Significant, fill = Domain)) +
  #plot the number of times each term occur, normalized
  geom_col() +
  coord_flip() +
  theme_bw() +
  #Aesthetic modifications
  theme(axis.text.x = element_text(hjust = 1, color = "black"), axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 6), legend.title = element_text(size = 8),
        legend.key.size = unit(0.2, "in"), plot.title = element_text(size = 11, hjust = 0.5),
        plot.margin = margin(10, 50, 10, 50)) +
  #Change name of the legend and legend labels
  scale_fill_jama(name = "Domain") +
  #add a plot title in bolds and italic
  labs(x = NULL, y = "Number of genes (significant)", title = "Downregulated genes") +
  theme(plot.title = element_text(face = "bold"))
DOWN_GOgraph

ggsave("Nicropus_switched_outputs/DOWN_enriched_GO.png", DOWN_GOgraph, dpi = 320)

source('../../get_legend.R')

legend <- get_legend(DOWN_GOgraph)

GOGO <- plot_grid(UP_GOgraph + 
                    theme(legend.position = 'none',
                          axis.title.x = element_blank()), 
                  DOWN_GOgraph + 
                    theme(legend.position = 'none',
                          plot.background = ),
                  nrow = 2,
                  align = TRUE,
                  labels = 'AUTO')
GOGO2 <- plot_grid(GOGO + 
                     theme(plot.margin = unit(c(0,-1,0,0), 'cm')),
                   legend,
                   rel_widths = c(3,0.4)) +
  theme(plot.margin = unit(c(0,1,0,0), 'cm'))

title <- ggdraw() +
  draw_label('Enriched GO terms',
             fontface = 'bold',
             size = 16)

GOGO2 <- plot_grid(title,
                   GOGO2,
                   ncol = 1,
                   rel_heights = c(0.07,1))
GOGO2

ggsave('Nicropus_switched_outputs/GOGO2.tiff', GOGO2, width = 7.5, height = 10, 
       units = 'in', dpi = 300, device = 'tiff')

# ---- KEGG Pathway analysis ----

## Pathways enriched in differentially expressed genes

geneListTemp <- top_sig_results_annotated %>% 
  select(transcript_id, KO, logFC, FDR) %>% 
  drop_na()

transcript2ko <- read_tsv('Nicropus_ko_insects.terms') %>% 
  select(KO, transcript_id) %>% 
  drop_na()

geneList <- geneListTemp$logFC
names(geneList) <- geneListTemp$transcript_id
geneList <- sort(geneList, decreasing = TRUE)
gene4enrich <- names(geneList)

kegg_enriched <- enricher(gene4enrich,
                              TERM2NAME = path_final,
                              TERM2GENE = transcript2ko,
                              minGSSize = 3)


Kegg_all_graph <- dotplot(kegg_enriched,
                              showCategory = 10) +
  labs(title = 'Overrepresented KEGG Pathways') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

Kegg_all_graph
ggsave('Nicropus_switched_outputs/KEGG_Pathways.png', Kegg_all_graph)


geneListUP <- geneList[geneList>0]
geneListUP <- sort(geneListUP, decreasing = TRUE)

# kegg_gsea <- GSEA(geneList,
#                   TERM2NAME = path_final,
#                   TERM2GENE = transcript2ko,
#                   minGSSize = 3)
# 
# Kegg_gsea_graph <- dotplot(kegg_gsea,
#                            showCategory = 10) +
#   labs(title = 'Enriched KEGG Pathways') +
#   theme(plot.title = element_text(face = 'bold', hjust = 0.5))

gene4enrichUP <- names(geneListUP)

kegg_enrichedUP <- enricher(gene4enrichUP,
                          TERM2NAME = path_final,
                          TERM2GENE = transcript2ko,
                          minGSSize = 3)


Kegg_all_graphUP <- dotplot(kegg_enrichedUP,
                          showCategory = 10) +
  labs(title = 'Overrepresented KEGG Pathways (upregulated genes)') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

Kegg_all_graphUP
ggsave('Nicropus_switched_outputs/KEGG_PathwaysUP.png', Kegg_all_graphUP)


geneListDOWN <- geneList[geneList<0]
geneListDOWN <- sort(geneListDOWN, decreasing = TRUE)

gene4enrichDOWN <- names(geneListDOWN)

kegg_enrichedDOWN <- enricher(gene4enrichDOWN,
                              TERM2NAME = path_final,
                              TERM2GENE = transcript2ko,
                              minGSSize = 3)


Kegg_all_graphDOWN <- dotplot(kegg_enrichedDOWN,
                              showCategory = 10) +
  labs(title = 'Overrepresented KEGG Pathways (downregulated genes)') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

Kegg_all_graphDOWN
ggsave('Nicropus_switched_outputs/KEGG_PathwaysDOWN.png', Kegg_all_graphDOWN)

KEGG_2_graph <- plot_grid(Kegg_all_graphUP,
                          Kegg_all_graphDOWN,
                          ncol = 1,
                          align = TRUE)
KEGG_2_graph
ggsave('Nicropus_switched_outputs/KEGG_PATH2.tiff', KEGG_2_graph,width = 10, height = 11, 
       units = 'in', dpi = 300, device = 'tiff')
