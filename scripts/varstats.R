#!/usr/bin/env Rscript

# Load required libraries (auto-install if missing)
suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
  library(ggplot2)
  library(patchwork)
  library(optparse)
})

# 1. Define command line options
option_list <- list(
  make_option(c("-o", "--output"), type = "character", default = "sv_comparative_dist.png",
              help = "Path to output image file [default %default]", metavar = "character")
)

parser <- OptionParser(usage = "Usage: %prog file.VCF [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)

if (length(args$args) < 1) {
  print_help(parser)
  stop("Missing input VCF file.", call. = FALSE)
}

vcf_file <- args$args[1]
output_file <- args$options$output

# 2. Robust VCF Reading
cat("Reading VCF:", vcf_file, "...\n")
all_lines <- readLines(vcf_file)
header_idx <- grep("^#CHROM", all_lines)

if (length(header_idx) == 0) {
  stop("Error: Could not find header line starting with #CHROM")
}

vcf_data <- read.table(text = all_lines[header_idx:length(all_lines)], 
                       header = TRUE, sep = "\t", comment.char = "", 
                       check.names = FALSE, stringsAsFactors = FALSE, quote = "")

# 3. Parse INFO for SVTYPE and SVLEN
extract_field <- function(info_column, field) {
  pattern <- paste0(field, "=([^;]+)")
  matches <- regexec(pattern, info_column)
  extracted <- regmatches(info_column, matches)
  sapply(extracted, function(x) if(length(x) > 1) x[2] else NA)
}

vcf_data$SVTYPE <- extract_field(vcf_data$INFO, "SVTYPE")
vcf_data$SVLEN  <- as.numeric(extract_field(vcf_data$INFO, "SVLEN"))

# 4. Filter and Prep Data
vcf_data$ABS_LEN <- abs(vcf_data$SVLEN)
relevant_types <- c("INS", "DEL", "DUP", "INV")
filtered_data <- vcf_data[vcf_data$SVTYPE %in% relevant_types & !is.na(vcf_data$ABS_LEN), ]

if (nrow(filtered_data) == 0) {
  stop("No valid INS, DEL, DUP, or INV variants found with SVLEN tags.")
}

# Calculate global max length for unified X-axis
global_max_len <- max(filtered_data$ABS_LEN, na.rm = TRUE)

# 5. Create Scatterplots with Shared X and Relative Y
plot_list <- list()

for (type in relevant_types) {
  type_data <- subset(filtered_data, SVTYPE == type)
  
  if (nrow(type_data) == 0) {
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = paste("No", type, "found")) +
      theme_void() + ggtitle(paste("Type:", type))
  } else {
    # Generate relative index for the Y axis
    type_data$idx <- 1:nrow(type_data)
    
    p <- ggplot(type_data, aes(x = ABS_LEN, y = idx)) +
      geom_point(alpha = 0.5, color = "midnightblue", size = 1) +
      # Use unified X scale but allow Y to scale to the number of variants
      coord_cartesian(xlim = c(0, global_max_len)) +
      theme_minimal() +
      labs(title = paste("Type:", type),
           x = "Length (bp)",
           y = "Count (Relative Index)") +
      theme(plot.title = element_text(face = "bold"))
  }
  plot_list[[type]] <- p
}

# 6. Combine and Save
cat("Combining plots...\n")
combined_plot <- (plot_list[["INS"]] | plot_list[["DEL"]]) / 
                 (plot_list[["DUP"]] | plot_list[["INV"]]) +
  plot_annotation(title = paste("SV Distribution:", basename(vcf_file)),
                  subtitle = "X-axis (Length) is fixed across all plots; Y-axis is relative to variant count.",
                  theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))

ggsave(output_file, combined_plot, width = 12, height = 9, dpi = 300)
cat("Plot saved to:", output_file, "\n")
