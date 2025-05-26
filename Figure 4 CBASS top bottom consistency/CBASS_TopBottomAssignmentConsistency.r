# Script by Christian R Voolstra, Luigi Colin @ University of Konstanz, Germany; 20250523
# Simple R script for analysis of CBASS ED50 top/bottom assignment consistency
# Reads data and checks correspondence of top/bottom assignment based on ED50s of consecutive CBASS runs
# Script will read csv file

# --- Setup --- #
# Load required libraries
library(pacman)
p_load(readr, Hmisc, data.table, ggplot2, ggpubr, tidyr)

# --- Clear workspace ---
rm(list = ls())

# --- Parameters --- #
# Set working directory and input/output paths. Input and output folder are relative to the working directory. # Default reads and writes all in the working directory.
work_dir <- rstudioapi::selectDirectory()   # Select working directory. if working outside RStudio, set manually. e.g.; work_dir <- "/path/to/working/directory"
if (!exists("work_dir")) work_dir <- "."    # Set working directory if it fails to set in RStudio.
out_path <- "./"                            # Default reads and writes all in the working directory.
input_folder <- "./"                        # Needs trailing "/". # Default reads and writes all in the working directory.

prefix <- ""                                # Can also be empty ""
item_per_group <- 5                         # Number of items in each group (top and bottom)
# --- End Parameters --- #

# --- Functions --- #
# Read data
# Function that reads all the CSV in a folder (for "multiple dataset" processing)
read_data_cbass_com <- function(path_to_csv_folder) {
    csv_files <- list.files(path_to_csv_folder, pattern = "\\.csv$", full.names = TRUE)
    data_list <- lapply(csv_files, function(f) {
        dat <- read_csv(f)
        colnames(dat)[1:3] <- c("ID", "ED50_R1", "ED50_R2")
        dat$meanED50 <- rowMeans(dat[, 2:3])
        dat <- dat[order(dat$ID), ]
        dat
    })
    names(data_list) <- tools::file_path_sans_ext(basename(csv_files))
    return(data_list)
}

# --- plotting theme  --- #
theme <- theme_minimal(base_size = 10) +
        theme(
            strip.text.y = element_blank(),
            strip.text.x = element_text(size = 10),
            panel.grid.minor = element_blank(),  # Hide minor grid lines
            panel.grid.major = element_blank(),  # Hide minor grid lines
            legend.key.size = unit(1.5, "lines"),  # Increase legend key size
            panel.spacing.y = unit(1.2, "lines"),  # Increase spacing between vertical panels
            panel.border = element_blank(),        # Remove full box border
            axis.ticks.length = unit(0.1, "cm"),
            axis.ticks = element_line(color = "black"),
            axis.line = element_line(color = "black") # Only x and y axis lines
        )
# --- End plotting theme --- #
# --- End Parameters --- #
# --- End Setup --- #

data_sets <- read_data_cbass_com(paste0(work_dir, "/", input_folder))

results_list <- list()
for (dataset in names(data_sets)) {
    # --- Initialize Output Data Frame --- #
    ResOutputall <- data.frame(
        numsamples = numeric(),
        numoverlapping = numeric(),
        Ttestpvalue = numeric(),
        Wilcoxtestpvalue = numeric()
    )

    # --- Main Simulation Loop --- #
    # --- Random Sampling and Analysis --- #
    for (i in ((item_per_group * 2):nrow(data_sets[[dataset]]))) {
        for (j in 1:1000) {
            # Randomly sample i rows
            x <- data_sets[[dataset]][sample(nrow(data_sets[[dataset]]), i), ]
            # gets the top and bottom N genets. Adding up the number of overlapping genets between the two groups.
            R1_top <- x[order(x$ED50_R1, decreasing = TRUE)[1:item_per_group], ]
            R1_bot <- x[order(x$ED50_R1, decreasing = FALSE)[1:item_per_group], ]
            R1_top_genets <- R1_top$ID
            R1_bot_genets <- R1_bot$ID
            R2_top5_genets <- x$ID[order(x$ED50_R2, decreasing = TRUE)[1:item_per_group]]
            R2_bot5_genets <- x$ID[order(x$ED50_R2, decreasing = FALSE)[1:item_per_group]]
            numoverlaptopbottom <- mean(c(
                sum(R1_top_genets %in% R2_bot5_genets),
                sum(R1_bot_genets %in% R2_top5_genets)
            ))

        # Store results
        ResOutputall <- rbind(
            ResOutputall,
            data.frame(
                numsamples = i,
                numoverlapping = numoverlaptopbottom
            )
            )
        }
    }
    results_list[[dataset]] <- ResOutputall
}

Output_Summary_all_list <- list()

for (dataset in names(data_sets)) {
    # --- Summarize Results --- #
    OutputSummaryall <- data.frame(
        numsamples = (item_per_group * 2):nrow(data_sets[[dataset]]),
        meanoverlap = aggregate(numoverlapping ~ numsamples, mean, data = ResOutputall)[[2]],
        minoverlap = aggregate(numoverlapping ~ numsamples, min, data = ResOutputall)[[2]],
        maxoverlap = aggregate(numoverlapping ~ numsamples, max, data = ResOutputall)[[2]],
        medianoverlap = aggregate(numoverlapping ~ numsamples, median, data = ResOutputall)[[2]],
        stddeverroroverlap = aggregate(numoverlapping ~ numsamples, function(x) sd(x) / sqrt(length(x)), data = ResOutputall)[[2]]
    )

# --- Save Results to file list --- #
Output_Summary_all_list[[dataset]] <- OutputSummaryall
}

# --- Plotting  --- #
for (dataset in names(data_sets)) {
  
  # Errorbar plot for mean overlap ± 1 stddev
  p3 <- ggplot(Output_Summary_all_list[[dataset]], aes(x = numsamples, y = meanoverlap)) +
    geom_line(linewidth = 0.3) +
    geom_point(size = 0.5) +
    geom_errorbar(
      aes(
        ymin = meanoverlap - stddeverroroverlap,
        ymax = meanoverlap + stddeverroroverlap
      ),
      width = 0.2
    ) +
    ylim(0,2.5)+
    theme +
    labs(
      x = "Number of colonies",
      y = "Mean overlap ± 1 stderror",
      title = paste0("Mean overlap top ", item_per_group," vs. bottom ", item_per_group," ED50s")
    )
  
  # 4. Save plots
  pdf_name <- if (prefix == "") {
    paste0(work_dir, "/", out_path, dataset, "TopVsBottom.pdf")
#  } else {
#    paste0(work_dir, "/", out_path, prefix, "_", dataset, ".pdf")
  }
  
  # Save each plot separately
  pdf(pdf_name, width = 6, height = 4)
  print(p3)
  dev.off()
}
# --- End of script --- #
