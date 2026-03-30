# Read your count matrix
# Read count matrix
countdata <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))

# Clean rownames
clean_names <- gsub(".*\\|", "", rownames(countdata))

# Your genes
my_genes <- readLines("geneIDS.txt")
head(my_genes)
# Extract just your genes - raw counts, no filtering
my_raw_counts <- countdata[clean_names %in% my_genes, ]

# See them all
print(my_raw_counts)

# Check which ones have zero across all samples
row_totals <- rowSums(my_raw_counts)
print(row_totals)

# Which of your genes would be removed by the >20 filter?
print("Genes that would be REMOVED by filtering:")
print(my_raw_counts[row_totals <= 20, ])

print("Genes that would be KEPT:")
print(my_raw_counts[row_totals > 20, ])

# Save all raw counts for your genes
write.csv(my_raw_counts, "//wsl.localhost/Ubuntu-24.04/home/spurthi/Quercussuber_scripts/my_LBD_raw_counts.csv")


##partB: creates qslbd to the following LoC gene IDs
library(readxl)
lbd_info <- read_excel("C:/Users/spurt/OneDrive - Auburn University/Desktop/Research/1_Genesequence_oak/oak/geneexcel.xlsx")

# Keep only QsLBD rows
lbd_names <- data.frame(gene_id = lbd_info$geneId, Name = lbd_info$Name, stringsAsFactors = FALSE)
lbd_names <- lbd_names[grepl("^QsLBD", lbd_names$Name), ]

# Read raw counts
my_raw <- read.csv("//wsl.localhost/Ubuntu-24.04/home/spurthi/Quercussuber_scripts/my_LBD_raw_counts.csv", row.names = 1)
my_raw$gene_id <- gsub(".*\\|", "", rownames(my_raw))

# Merge
final <- merge(lbd_names, my_raw, by = "gene_id")

# Sort by QsLBD number
final <- final[order(as.numeric(gsub("QsLBD", "", final$Name))), ]

print(final)


# Save
write.csv(final, "//wsl.localhost/Ubuntu-24.04/home/spurthi/Quercussuber_scripts/my_LBD_raw_counts_named.csv", row.names = FALSE)

#partC-createheatmap
library(pheatmap)
library(RColorBrewer)

# Read your edited file
lbd <- read.csv("my_LBD_raw_counts_named.csv")

# Set up count matrix
counts <- as.matrix(lbd[, 3:ncol(lbd)])
rownames(counts) <- lbd$Name

# Column annotation
anno_col <- data.frame(
  Tissue = c(rep("Root_Tip", 4), rep("Root_Segment", 4)),
  row.names = colnames(counts)
)

# Define colors for tissue
ann_colors <- list(Tissue = c(Root_Tip = "#E41A1C", Root_Segment = "#377EB8"))

# Log2 transform
log_counts <- log2(counts + 1)

library(pheatmap)
library(RColorBrewer)
# ---- HEATMAP 1: All genes, row-scaled (shows pattern) ----
pheatmap(log_counts,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         scale = "none",
         main = "QsLBD All Genes Expression (log2)",
         fontsize_row = 10,
         fontsize_col = 10,
         cellwidth = 40,
         cellheight = 15,
         gaps_col = 4,
         color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100))
# ---- HEATMAP 2: Only expressed genes (remove all-zero rows) ----
expressed <- log_counts[rowSums(counts) > 0, ]

ann_colors <- list(Tissue = c(Root_Tip = "#2ECC71", Root_Segment = "#9B59B6"))

pheatmap(expressed,
         annotation_col = anno_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         scale = "row",
         main = "QsLBD Expressed Genes (row-scaled)",
         fontsize_row = 9,
         cellwidth = 40,
         cellheight = 15,
         color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100))

# ---- HEATMAP 3: Averaged by tissue ----
avg <- data.frame(
  Root_Tip = rowMeans(counts[, 1:4]),
  Root_Segment = rowMeans(counts[, 5:8])
)
rownames(avg) <- lbd$Name
avg <- avg[rowSums(avg) > 0, ]
log_avg <- log2(as.matrix(avg) + 1)

pheatmap(log_avg,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         main = "QsLBD Average Expression by Tissue (log2)",
         fontsize_row = 9,
         color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100))