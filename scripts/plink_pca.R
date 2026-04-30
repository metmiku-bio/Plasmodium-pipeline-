#!/usr/bin/env Rscript

# plink_pca.R - PCA and tree construction from PLINK distance matrix
# Updated for Nextflow pipeline integration

# Load required libraries
library(showtext)
library(dplyr)
library(ggplot2)
library(ape)
library(optparse)

showtext_auto()

# Parse command-line arguments
parser <- OptionParser(description = 'Run PCA and tree construction from PLINK distance matrix')
parser$add_argument('--workdir', required=TRUE, help='Working directory containing PLINK distance files')
parser$add_argument('--prefix', required=TRUE, help='Prefix for PLINK .dist and .dist.id files')
parser$add_argument('--metadata', required=TRUE, help='Metadata file path (tab-separated)')
parser$add_argument('--color_by', required=FALSE, default="country", help="Column to be colored by")
parser$add_argument('--pc_count', required=FALSE, type="integer", default=10, help="Number of PCs to compute")
parser$add_argument('--output_dir', required=FALSE, default=".", help="Output directory for plots")

args <- parser$parse_args()

workdir <- args$workdir
prefix <- args$prefix
metadata_file <- args$metadata
color_by <- args$color_by
pc_count <- args$pc_count
output_dir <- args$output_dir

# Function to calculate variance explained
calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    return(vars)
}

# Function to create PCA plot
create_pca_plot <- function(df, vars, pc_x, pc_y, color_by, output_file) {
    p <- ggplot(data = df, aes(x = !!sym(pc_x), y = !!sym(pc_y),
                                color = !!sym(color_by))) +
        geom_point(size = 3, alpha = 0.7) +
        labs(x = paste0(pc_x, " (", vars[pc_x], "%)"),
             y = paste0(pc_y, " (", vars[pc_y], "%)"),
             title = paste("PCA:", pc_x, "vs", pc_y),
             color = color_by) +
        theme_classic() +
        theme(
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.title = element_text(face = "bold"),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold")
        ) +
        scale_color_viridis_d()
    
    ggsave(plot = p, filename = output_file, width = 10, height = 8, dpi = 300)
    cat(paste("Saved plot:", output_file, "\n"))
}

# Print startup message
cat("\n========================================\n")
cat("PCA ANALYSIS WITH R\n")
cat("========================================\n")
cat(paste("Working directory:", workdir, "\n"))
cat(paste("PLINK prefix:", prefix, "\n"))
cat(paste("Metadata file:", metadata_file, "\n"))
cat(paste("Color by:", color_by, "\n"))
cat(paste("Number of PCs:", pc_count, "\n"))
cat(paste("Output directory:", output_dir, "\n"))
cat("========================================\n\n")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# 1. Read metadata
cat("Reading metadata...\n")
met <- read.table(metadata_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
cat(paste("Loaded metadata for", nrow(met), "samples\n"))

# 2. Read PLINK distance matrix files
dist_file <- file.path(workdir, paste0(prefix, ".dist"))
dist_id_file <- file.path(workdir, paste0(prefix, ".dist.id"))

# Check if files exist
if (!file.exists(dist_file)) {
    stop(paste("Distance matrix file not found:", dist_file))
}
if (!file.exists(dist_id_file)) {
    stop(paste("Distance matrix ID file not found:", dist_id_file))
}

cat("Reading PLINK distance matrix...\n")
dist <- read.table(dist_file, header = FALSE)
id <- read.table(dist_id_file, stringsAsFactors = FALSE)

cat(paste("Loaded distance matrix for", nrow(dist), "samples\n"))

# 3. Merge with metadata
# Try different column name possibilities
if ("Sample" %in% colnames(met)) {
    desc <- id %>% left_join(met, by = c("V1" = "Sample"))
} else if ("sample_id" %in% colnames(met)) {
    desc <- id %>% left_join(met, by = c("V1" = "sample_id"))
} else if ("ID" %in% colnames(met)) {
    desc <- id %>% left_join(met, by = c("V1" = "ID"))
} else {
    cat("Warning: No matching sample ID column found in metadata\n")
    cat("Available columns:", paste(colnames(met), collapse=", "), "\n")
    desc <- id
    desc$country <- "Unknown"
}

# 4. Create distance matrix
dist_m <- as.matrix(dist)
colnames(dist_m) <- desc$V1
rownames(dist_m) <- desc$V1

cat("Distance matrix dimensions:", dim(dist_m)[1], "x", dim(dist_m)[2], "\n")

# 5. Perform PCA (Multidimensional Scaling)
cat("Performing PCA/MDS...\n")
cmd <- cmdscale(dist_m, k = min(pc_count, nrow(dist_m)-1), eig = TRUE, x.ret = TRUE)

# 6. Calculate variance explained
vars <- calc_variance_explained(cmd)
cat("Variance explained by each PC:\n")
for (i in 1:min(length(vars), 10)) {
    cat(paste("  PC", i, ":", vars[paste0("PC", i)], "%\n", sep=""))
}

# 7. Create PCA data frame
df <- as.data.frame(cmd$points, stringsAsFactors = FALSE)
colnames(df) <- paste0("PC", 1:ncol(df))

# Add metadata to the PCA points
if (nrow(df) == nrow(desc)) {
    df$sample_id <- desc$V1
    
    # Add color_by column if it exists in metadata
    if (color_by %in% colnames(met)) {
        df[[color_by]] <- desc[[color_by]]
    } else {
        cat(paste("Warning: Column '", color_by, "' not found in metadata. Using default.\n", sep=""))
        df[[color_by]] <- "Group 1"
    }
} else {
    cat("Warning: PCA points and metadata dimensions don't match\n")
    df$sample_id <- rownames(dist_m)
    df[[color_by]] <- "Group 1"
}

# 8. Create PCA plots
cat("\nCreating PCA plots...\n")

# PC1 vs PC2
if (ncol(df) >= 2) {
    plot_file <- file.path(output_dir, "pca_plot_PC1_PC2.pdf")
    create_pca_plot(df, vars, "PC1", "PC2", color_by, plot_file)
    
    # Also create PNG version
    plot_file_png <- file.path(output_dir, "pca_plot_PC1_PC2.png")
    create_pca_plot(df, vars, "PC1", "PC2", color_by, plot_file_png)
}

# PC1 vs PC3
if (ncol(df) >= 3) {
    plot_file <- file.path(output_dir, "pca_plot_PC1_PC3.pdf")
    create_pca_plot(df, vars, "PC1", "PC3", color_by, plot_file)
    
    plot_file_png <- file.path(output_dir, "pca_plot_PC1_PC3.png")
    create_pca_plot(df, vars, "PC1", "PC3", color_by, plot_file_png)
}

# PC2 vs PC3
if (ncol(df) >= 3) {
    plot_file <- file.path(output_dir, "pca_plot_PC2_PC3.pdf")
    create_pca_plot(df, vars, "PC2", "PC3", color_by, plot_file)
    
    plot_file_png <- file.path(output_dir, "pca_plot_PC2_PC3.png")
    create_pca_plot(df, vars, "PC2", "PC3", color_by, plot_file_png)
}

# 9. Create variance explained bar plot
cat("Creating variance explained plot...\n")
var_df <- data.frame(
    PC = factor(names(vars)[1:min(10, length(vars))], levels = names(vars)[1:min(10, length(vars))]),
    Variance = vars[1:min(10, length(vars))]
)

var_plot <- ggplot(var_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = paste0(Variance, "%")), vjust = -0.5, size = 3.5) +
    labs(title = "Variance Explained by Principal Components",
         x = "Principal Component",
         y = "Variance Explained (%)") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))

var_plot_file <- file.path(output_dir, "variance_explained.pdf")
ggsave(plot = var_plot, filename = var_plot_file, width = 10, height = 6, dpi = 300)

# 10. Create scree plot
scree_plot <- ggplot(var_df, aes(x = PC, y = Variance, group = 1)) +
    geom_line(color = "steelblue", size = 1.2) +
    geom_point(color = "steelblue", size = 3) +
    labs(title = "Scree Plot",
         x = "Principal Component",
         y = "Variance Explained (%)") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))

scree_file <- file.path(output_dir, "scree_plot.pdf")
ggsave(plot = scree_plot, filename = scree_file, width = 8, height = 6, dpi = 300)

# 11. Export distance matrix to Newick format for Neighbour Joining tree
cat("Building Neighbour Joining tree...\n")
tree <- nj(dist_m)
newick_file <- file.path(output_dir, paste0(prefix, ".newick"))
write.tree(phy = tree, file = newick_file)
cat(paste("Tree saved to:", newick_file, "\n"))

# 12. Save PCA results as RDS for future use
rds_file <- file.path(output_dir, paste0(prefix, "_pca_results.rds"))
saveRDS(list(cmd = cmd, vars = vars, df = df), file = rds_file)
cat(paste("PCA results saved to:", rds_file, "\n"))

# 13. Write summary file
summary_file <- file.path(output_dir, "pca_summary.txt")
sink(summary_file)
cat("PCA Analysis Summary\n")
cat("====================\n\n")
cat(paste("Date:", Sys.time(), "\n"))
cat(paste("Working directory:", workdir, "\n"))
cat(paste("PLINK prefix:", prefix, "\n"))
cat(paste("Metadata file:", metadata_file, "\n"))
cat(paste("Color by:", color_by, "\n"))
cat(paste("Number of PCs computed:", pc_count, "\n"))
cat(paste("Number of samples:", nrow(dist_m), "\n\n"))
cat("Variance explained by each PC:\n")
for (i in 1:min(length(vars), pc_count)) {
    cat(sprintf("  PC%d: %.1f%%\n", i, vars[paste0("PC", i)]))
}
cat("\nOutput files:\n")
cat(paste("  - PCA plots:", output_dir, "/pca_plot_*.pdf\n", sep=""))
cat(paste("  - Tree file:", newick_file, "\n"))
cat(paste("  - RDS file:", rds_file, "\n"))
sink()

cat("\n========================================\n")
cat("PCA ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("========================================\n")
cat(paste("Summary file:", summary_file, "\n"))
cat(paste("Results saved to:", output_dir, "\n"))
cat("========================================\n\n")