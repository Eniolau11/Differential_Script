# Load necessary libraries
library(DESeq2)  
library(tidyverse)  
library(pheatmap)  
library(RColorBrewer)  
library(ggrepel)  

# Load data
count_data <- read.csv("output.csv", header = TRUE, row.names = 1)  # Read count data from CSV file
sample_info <- read.csv("colData.csv", header = TRUE, row.names = 1) %>%
  mutate(across(c(Sex, Age), as.factor))  # Read sample info and convert Sex and Age columns to factors

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info,
                              design = ~ Sex + Age + Sex:Age)  

# Perform DESeq2 analysis
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ Sex + Age)  # Run DESeq2 using likelihood ratio test

# Convert results to dataframe and add row names as Gene_name
res_lrt <- results(dds_lrt) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene_name") %>%  # Convert results to a dataframe and create a column for gene names
  arrange(pvalue)  # Arrange results by p-value

# Filter based on p-value
filtered_0.01 <- res_lrt %>% filter(pvalue < 0.01)  # Filter results to include only genes with p-value < 0.01

# Plotting histogram for log2fold change
filtered_0.01_results <- filtered_0.01[!is.na(filtered_0.01$log2FoldChange), ]  # Remove rows with NA in log2FoldChange
ggplot(filtered_0.01_results, aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +  # Plot histogram of log2FoldChanges
  labs(title = "Histogram of log2 Fold Changes",
       x = "log2 Fold Change",
       y = "Frequency") +
  theme_minimal() 

# Filter for significant log2 fold change
filtered_Log <- filtered_0.01 %>% filter(abs(log2FoldChange) > 2)  # Filter results for absolute log2FoldChange > 2

# Save DESeq2 results
#write.csv(res_lrt, 'DESeq_results.csv')  # Save the DESeq2 results to a CSV file

# Normalise counts and add Gene_name
normalised_counts <- counts(dds_lrt, normalized = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene_name")  # Normalise counts and convert to dataframe with gene names

# Select top 10 genes based on adjusted p-value and a log2 fold change of 2
top_genes <- res_lrt %>%
  filter(pvalue < 0.1, abs(log2FoldChange) >= 2) %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  pull(Gene_name)  # Identify top 10 significant genes

# Filter normalised counts for the top genes
top_gene_data <- normalised_counts %>%
  filter(Gene_name %in% top_genes)  # Filter normalised count data for top genes

# Prepare data for plotting
combined_data <- top_gene_data %>%
  pivot_longer(cols = -Gene_name, names_to = "Sample", values_to = "Expression") %>%
  mutate(
    Age = as.numeric(gsub("X(\\d+).*", "\\1", Sample)),
    Sex = ifelse(grepl("F", Sample), "Female", "Male")
  ) %>%
  group_by(Gene_name, Age, Sex) %>%
  summarize(
    Mean_Expression = mean(Expression, na.rm = TRUE),
    Standard_Deviation = sd(Expression, na.rm = TRUE),
    .groups = 'drop'
  )  # Combine data by gene, age, and sex for plotting

# Save combined data
#write.csv(combined_data, "Combined_gene_data.csv")  
combined_data$Age <- factor(combined_data$Age, levels = c(0, 8, 16))  # Factorise age for plotting

# For loop to Iterate through unique gene names
for (gene_name in unique(combined_data$Gene_name)) {
  gene_data <- subset(combined_data, Gene_name == gene_name)  # Subset data for current gene
  
  # Plotting (line plot)
  Gene_plots <- ggplot(gene_data, aes(x = Age, y = Mean_Expression, group = Sex, color = Sex)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = Mean_Expression - Standard_Deviation, ymax = Mean_Expression + Standard_Deviation),
                  width = 0.2) +
    labs(x = "Age", y = "Mean Expression", title = paste("Mean Expression vs. Age -", gene_name)) +
    theme_minimal()  # Create a line plot with points and error bars
  
  print(Gene_plots)
}

# PCA plot
vsd <- vst(dds_lrt, blind = FALSE)  # Variance stabilising transformation
plotPCA(vsd, intgroup = c("Sex", "Age"))  # PCA plot colored by Sex and Age

# HeatMap
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)  # Set color scheme
top_hits <- filtered_Log[order(filtered_Log$padj), ][1:10,]  # Select top 10 genes based on adjusted p-value
top_hits <- top_hits[["Gene_name"]]  # Extract gene names
rld <- rlog(dds_lrt, blind = FALSE)  # Regularised log transformation
annot_info <- as.data.frame(colData(dds_lrt)[, c("Age", "Sex")])  # Extract sample annotation information
pheatmap(assay(rld)[top_hits,], cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE,
         annotation_col = annot_info, col = colors)  # Plot heatmap with clustering and annotations

# MA plot
plotMA(dds_lrt, ylim = c(-10, 10))  # MA plot of log fold changes and mean expression
# Remove noise
resLFC <- lfcShrink(dds_lrt, coef = 2, type = "apeglm")  # Shrink log fold changes using apeglm method
plotMA(resLFC, ylim = c(-10, 10))  # Plot shrunk MA values




library(readr)
list_of_locs<-substr(filtered_Log$Gene_name,4,12)
library(mygene)
fulldata<-getGenes(list_of_locs)

write.csv(fulldata, file="differentially_expressed_genes_Age:Sex.csv")


