# ClusterProfiler

Link: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html

Let's use the DESeq2 output file as an example. In this file, we observe the following columns:

![image](https://github.com/user-attachments/assets/c5ee91e3-6c98-4225-b5a9-1a6f0a0cfa5e)

**gene_id:** Gene identification

**baseMean:** Average of normalized counts across all samples

**log2FoldChange:** Relative expression change, expressed in base 2 logarithmic scale

**IfcSE:** Standard error associated with log2FoldChange

**stat:** Test statistic applied

**pvalue:** P-value indicating statistical significance

**padj:** Adjusted p-value for multiple testing correction

**OBS:** It is possible to organize the output data based on adjusted p-values, commonly using a threshold of p-value > 0.05. 

Furthermore, the data can be separated to analyze differentially expressed genes (DEGs) that are either more or less expressed according to logfoldchange. 
This allows for the identification of enriched pathways that include both upregulated and downregulated genes. 
The analysis can also focus solely on upregulated DEGs and their associated enriched biological pathways, or solely on downregulated DEGs and their corresponding enriched pathways

# Script for KEGG (Kyoto Encyclopedia of Genes and Genomes)

library(clusterProfiler)

library(org.Hs.eg.db)

d <- read.csv("name_of_the_file.csv")

geneList <-d [,2]

names(geneList) <- as.character(d[,1])

gene <- names(geneList)

gene2 = bitr(gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

gene3 <- gene2[,2]

kk <- enrichKEGG(gene         = gene3,

organism     = 'hsa',
pvalueCutoff = 0.05)

y <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.csv(y, "results.csv")

dotplot(kk, showCategory=30)

# Script explanation

**library(clusterProfiler) and library(org.Hs.eg.db) -->** load the necessary packages for enrichment analysis and for the human genome database. It is important to mention that the genome database can change according to the organism type

**d <- read.csv("name_of_the_file.csv") -->** reads a CSV file containing gene data. The file should include two columns: one for gene identifiers and another for associated values

**geneList <-d [,2] -->** extracts the values from the second column (e.g., expression values) and **names(geneList) <- as.character(d[,1])** assigns the gene identifiers from the first column as names for the values in **geneList**

**gene <- names(geneList) -->** extracts the gene names from **geneList**

**bitr(gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") -->** converts the gene identifiers from ENSEMBL format to ENTREZ IDs using the human genome database **org.Hs.eg.db**

**gene3 <- gene2[,2] -->** extracts the ENTREZ IDs from the conversion result for further analysis

**kk <- enrichKEGG(gene         = gene3, organism     = 'hsa', pvalueCutoff = 0.05) -->** performs enrichment analysis on the gene list, identifying pathways in the KEGG database ( **hsa** is the KEGG code for human). The p-value cutoff filters pathways based on statistical significance

**y <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID") -->** converts the enrichment results into a human-readable format by mapping ENTREZ IDs back to gene symbols

**write.csv(y, "results.csv") -->** saves the formatted enrichment results into a CSV file named **results.csv**

**dotplot(kk, showCategory=30) -->** generates a dot plot to visualize the top 30 enriched pathways based on their significance

# Script for GO (Gene Onthology)

**OBS:** The beggining of the script is the same as KEGG analysis

library(clusterProfiler)

library(org.Hs.eg.db)

d <- read.csv("name_of_the_file.csv")

geneList <-d [,2]

names(geneList) <- as.character(d[,1])

gene <- names(geneList)

ego <- enrichGO(gene          = gene,
                keyType = "ENSEMBL", 
                OrgDb         = org.Hs.eg.db, 
                ont           = "BP", 
                pAdjustMethod = "BH", 
                pvalueCutoff  = 0.05, 
                qvalueCutoff  = 0.05, 
        readable      = TRUE) 

write.csv(ego, "results.csv")

dotplot(ego, showCategory=30)

# Script explanation

**enrichGO -->** The function name

**gene -->** The list of genes to be analyzed

**keyType = "ENSEMBL" -->** Specifies that the gene identifiers are in ENSEMBL format

**OrgDb= org.Hs.eg.db -->** The organism database to be used for gene annotation (in this case, human)

**ont= "BP" -->** Defines the ontology category for analysis (Biological Processes- BP)

**pAdjustMethod = "BH" -->** Adjustment method for p-values to correct for multiple testing (Benjamini-Hochberg method)

**pvalueCutoff  = 0.05 -->** Sets the threshold for p-values. Only terms with p-values below this limit are considered significant

**qvalueCutoff  = 0.05 -->** Sets the threshold for q-values after p-value adjustment

** readable      = TRUE -->** Converts the output into a human-readable format, mapping ENSEMBL IDs to gene names

**dotplot(ego, showCategory=30) -->** generates a dot plot to visualize the top 30 enriched pathways based on their significance

**write.csv(ego, "results.csv") -->** saves the formatted enrichment results into a CSV file named **results.csv**





