library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(ggrepel)  # for better label positioning
library(plotly)   # for interactive plots
library(RColorBrewer)

# Your ethanol extract GC-MS data
ethanol_extract_data <- data.frame(
  Peak = 1:20,
  R.Time = c(17.711, 17.821, 19.891, 21.111, 23.542, 24.671, 25.74, 25.952,
             26.921, 26.985, 27.333, 28.001, 28.234, 28.545, 28.972, 29.001,
             29.752, 33.793, 34.718, 34.94),
  Area = c(3324567, 4835205, 3974012, 585421, 3792104, 458015, 3324576, 544532,
           748731, 674847, 5327889, 3957392, 6797872, 20958615, 6520475, 55346572,
           720547, 3467297, 2084567, 44856387),
  Area_Percent = c(1.93, 2.81, 2.31, 0.34, 2.20, 0.27, 1.93, 0.32,
                   0.43, 0.39, 3.09, 2.30, 3.95, 12.16, 3.78, 32.12,
                   0.42, 2.01, 1.21, 26.03),
  Height = c(1247345, 2511025, 327563, 234110, 237502, 128075, 1021572, 27185,
             354221, 254626, 78569, 1254314, 2136250, 5432372, 2725371, 25537252,
             34528, 1533473, 482352, 14342175),
  Height_Percent = c(2.08, 4.19, 0.55, 0.39, 0.40, 0.21, 1.71, 0.05,
                     0.59, 0.43, 0.13, 2.09, 3.57, 9.07, 4.55, 42.63,
                     0.06, 2.56, 0.81, 23.94),
  A_H_Ratio = c(2.67, 1.93, 12.13, 2.50, 15.97, 3.58, 3.25, 20.03,
                2.11, 2.65, 67.81, 3.16, 3.18, 3.86, 2.39, 2.17,
                20.87, 2.26, 4.32, 3.13),
  Compound_Name = c("Pentanoic acid, 4-methyl-, pentyl ester",
                    "2-Hydroxy-1-(1'-pyrrolidiyl)-1-buten-3-one",
                    "Neophytadiene",
                    "Octadecyne-3",
                    "13-Methyltetradecanal",
                    "trans-2-Dodecen-1-ol",
                    "1-Hexadecyne",
                    "n-Hexadecanoic acid",
                    "Eicosane",
                    "Lactose",
                    "1-Trimethylsilyl-3-(dimethyl-n-pentylsilyl) but-1-ene",
                    "Phytol",
                    "Cyclodecasiloxane, eicosamethyl",
                    "Octadecanoic acid",
                    "Cyclododecanone",
                    "Oleic Acid",
                    "9-Octadecenamide",
                    "Hexadecanoic acid, 2-hydroxy-1-(hydroxymethyl)ethyl ester",
                    "Glycerol 1-palmitate",
                    "Stigmasterol")
)

# Function to create publication-quality GC-MS chromatogram
create_gcms_chromatogram <- function(data,
                                     plot_type = "area", # "area" or "height"
                                     color_by_intensity = TRUE,
                                     show_labels = "major", # "all", "major", "none"
                                     label_threshold = 5.0, # minimum area% for labeling
                                     title = "GC-MS Chromatogram - Ethanol Extract") {
 
  # Select intensity column based on plot_type
  if (plot_type == "area") {
    data$intensity <- data$Area
    data$intensity_percent <- data$Area_Percent
    y_label <- "Peak Area"
  } else {
    data$intensity <- data$Height
    data$intensity_percent <- data$Height_Percent
    y_label <- "Peak Height"
  }
 
  # Create color palette
  if (color_by_intensity) {
    colors <- colorRampPalette(c("lightblue", "steelblue", "darkblue", "red"))(nrow(data))
    data$color <- colors[rank(data$intensity_percent)]
  } else {
    data$color <- "steelblue"
  }
 
  # Create base plot
  p <- ggplot(data, aes(x = R.Time, y = intensity)) +
    geom_col(aes(fill = intensity_percent), width = 0.1, alpha = 0.8) +
    scale_fill_gradient2(low = "lightblue", mid = "steelblue", high = "red",
                         midpoint = median(data$intensity_percent),
                         name = paste0(tools::toTitleCase(plot_type), "\n(%)")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray60"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.position = "right",
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    ) +
    labs(
      title = title,
      subtitle = paste("Total Compounds Detected:", nrow(data)),
      x = "Retention Time (min)",
      y = y_label,
      caption = paste("Data source: Ethanol extract analysis •",
                      "Major peaks (≥5%):", sum(data$intensity_percent >= 5))
    ) +
    scale_x_continuous(breaks = pretty_breaks(n = 10),
                       limits = c(min(data$R.Time) - 0.5, max(data$R.Time) + 0.5)) +
    scale_y_continuous(labels = scales::scientific_format(digits = 2))
 
  # Add labels based on show_labels parameter
  if (show_labels == "all") {
    label_data <- data
  } else if (show_labels == "major") {
    label_data <- data[data$intensity_percent >= label_threshold, ]
  } else {
    label_data <- data.frame() # empty for no labels
  }
 
  if (nrow(label_data) > 0) {
    # Truncate long compound names for better display
    label_data$short_name <- sapply(label_data$Compound_Name, function(x) {
      if (nchar(x) > 25) {
        paste0(substr(x, 1, 22), "...")
      } else {
        x
      }
    })
   
    p <- p +
      geom_text_repel(data = label_data,
                      aes(x = R.Time, y = intensity,
                          label = paste0(short_name, "\n(",
                                         round(intensity_percent, 1), "%)")),
                      size = 3, color = "black",
                      box.padding = 0.3, point.padding = 0.2,
                      max.overlaps = 15, seed = 42)
  }
 
  return(p)
}

# Function to create compound distribution pie chart
create_compound_pie_chart <- function(data, top_n = 10) {
  # Select top N compounds and group others
  data_sorted <- data[order(-data$Area_Percent), ]
 
  if (nrow(data_sorted) > top_n) {
    top_compounds <- data_sorted[1:top_n, ]
    other_percent <- sum(data_sorted[(top_n + 1):nrow(data_sorted), ]$Area_Percent)
   
    pie_data <- rbind(top_compounds[, c("Compound_Name", "Area_Percent")],
                      data.frame(Compound_Name = paste("Others (", nrow(data_sorted) - top_n, " compounds)"),
                                 Area_Percent = other_percent))
  } else {
    pie_data <- data_sorted[, c("Compound_Name", "Area_Percent")]
  }
 
  # Truncate names
  pie_data$short_name <- sapply(pie_data$Compound_Name, function(x) {
    if (nchar(x) > 20) {
      paste0(substr(x, 1, 17), "...")
    } else {
      x
    }
  })
 
  # Create pie chart
  p <- ggplot(pie_data, aes(x = "", y = Area_Percent, fill = reorder(short_name, -Area_Percent))) +
    geom_col(width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      legend.position = "right"
    ) +
    labs(
      title = paste("Compound Distribution - Top", min(top_n, nrow(data))),
      fill = "Compounds"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set3")
 
  return(p)
}

# Function to create horizontal bar chart (alternative to pie chart)
create_compound_bar_chart <- function(data, top_n = 10) {
  # Select top N compounds
  data_sorted <- data[order(-data$Area_Percent), ][1:min(top_n, nrow(data)), ]
 
  # Truncate names
  data_sorted$short_name <- sapply(data_sorted$Compound_Name, function(x) {
    if (nchar(x) > 30) {
      paste0(substr(x, 1, 27), "...")
    } else {
      x
    }
  })
 
  # Create bar chart
  p <- ggplot(data_sorted, aes(x = reorder(short_name, Area_Percent), y = Area_Percent)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = paste("Top", min(top_n, nrow(data)), "Compounds by Peak Area"),
      x = "Compound Name",
      y = "Peak Area (%)"
    ) +
    geom_text(aes(label = paste0(round(Area_Percent, 1), "%")),
              hjust = -0.1, size = 3)
 
  return(p)
}

# Function to create summary statistics table
create_summary_table <- function(data) {
  cat("\n=== GC-MS ANALYSIS SUMMARY ===\n")
  cat("Dataset: Ethanol Extract Analysis\n")
  cat("Total number of compounds detected:", nrow(data), "\n")
  cat("Retention time range:", round(min(data$R.Time), 2), "-", round(max(data$R.Time), 2), "minutes\n")
  cat("Total peak area:", format(sum(data$Area), big.mark = ","), "\n")
  cat("Average peak area:", format(round(mean(data$Area)), big.mark = ","), "\n\n")
 
  # Major compounds (≥5% area)
  major_compounds <- data[data$Area_Percent >= 5, ]
  cat("=== MAJOR COMPOUNDS (≥5% Peak Area) ===\n")
  if (nrow(major_compounds) > 0) {
    for (i in 1:nrow(major_compounds)) {
      cat(sprintf("Peak %d: %s\n", major_compounds$Peak[i], major_compounds$Compound_Name[i]))
      cat(sprintf("  Retention Time: %.3f min\n", major_compounds$R.Time[i]))
      cat(sprintf("  Peak Area: %.2f%%\n", major_compounds$Area_Percent[i]))
      cat(sprintf("  Peak Height: %.2f%%\n", major_compounds$Height_Percent[i]))
      cat("\n")
    }
  }
 
  # Top 5 by area
  top5 <- data[order(-data$Area_Percent), ][1:5, ]
  cat("=== TOP 5 COMPOUNDS BY PEAK AREA ===\n")
  print(top5[, c("Peak", "R.Time", "Area_Percent", "Compound_Name")], row.names = FALSE)
 
  return(invisible())
}

# Function to create interactive plotly visualization
create_interactive_chromatogram <- function(data, plot_type = "area") {
 
  if (plot_type == "area") {
    data$intensity <- data$Area
    data$intensity_percent <- data$Area_Percent
    y_label <- "Peak Area"
  } else {
    data$intensity <- data$Height
    data$intensity_percent <- data$Height_Percent
    y_label <- "Peak Height"
  }
 
  # Create custom hover text
  data$hover_text <- paste(
    "Peak:", data$Peak,
    "\nCompound:", data$Compound_Name,
    "\nRetention Time:", data$R.Time, "min",
    "\nPeak Area:", round(data$Area_Percent, 2), "%",
    "\nPeak Height:", round(data$Height_Percent, 2), "%",
    "\nA/H Ratio:", round(data$A_H_Ratio, 2)
  )
 
  # Create ggplot with points and segments instead of columns
  p <- ggplot(data, aes(x = R.Time, y = intensity, text = hover_text)) +
    # Add vertical lines from baseline to peak height
    geom_segment(aes(x = R.Time, xend = R.Time, y = 0, yend = intensity,
                     color = intensity_percent),
                 linewidth = 2, alpha = 0.8) +
    # Add points at peak tops
    geom_point(aes(color = intensity_percent),
               size = 3, alpha = 0.9) +
    scale_color_gradient2(low = "lightblue", mid = "steelblue", high = "red",
                          midpoint = median(data$intensity_percent),
                          name = paste0(tools::toTitleCase(plot_type), "\n(%)")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 10)
    ) +
    labs(
      title = "Interactive GC-MS Chromatogram - Ethanol Extract",
      x = "Retention Time (min)",
      y = y_label
    ) +
    scale_y_continuous(labels = scales::scientific_format())
 
  # Convert to plotly
  interactive_plot <- ggplotly(p, tooltip = "text")
 
  # Customize plotly layout
  interactive_plot <- interactive_plot %>%
    layout(
      hovermode = 'closest',
      showlegend = TRUE,
      title = list(
        text = "Interactive GC-MS Chromatogram - Ethanol Extract",
        font = list(size = 16)
      )
    )
 
  return(interactive_plot)
}

# Alternative function using native plotly (no ggplot conversion)
create_native_plotly_chromatogram <- function(data, plot_type = "area") {
 
  if (plot_type == "area") {
    data$intensity <- data$Area
    data$intensity_percent <- data$Area_Percent
    y_label <- "Peak Area"
  } else {
    data$intensity <- data$Height
    data$intensity_percent <- data$Height_Percent
    y_label <- "Peak Height"
  }
 
  # Create color scale
  colors <- colorRampPalette(c("lightblue", "steelblue", "red"))(100)
  color_indices <- round(scales::rescale(data$intensity_percent, to = c(1, 100)))
  data$color <- colors[color_indices]
 
  # Create plotly figure
  fig <- plot_ly(
    data = data,
    x = ~R.Time,
    y = ~intensity,
    type = 'bar',
    marker = list(color = ~color),
    hovertemplate = paste(
      "<b>%{customdata[0]}</b><br>",
      "Peak: %{customdata[1]}<br>",
      "Retention Time: %{x:.3f} min<br>",
      "Peak Area: %{customdata[2]:.2f}%<br>",
      "Peak Height: %{customdata[3]:.2f}%<br>",
      "A/H Ratio: %{customdata[4]:.2f}<br>",
      "<extra></extra>"
    ),
    customdata = ~cbind(Compound_Name, Peak, Area_Percent, Height_Percent, A_H_Ratio)
  ) %>%
    layout(
      title = list(
        text = "Interactive GC-MS Chromatogram - Ethanol Extract",
        font = list(size = 16)
      ),
      xaxis = list(title = "Retention Time (min)"),
      yaxis = list(title = y_label, tickformat = ".2e"),
      showlegend = FALSE,
      hovermode = 'closest'
    )
 
  return(fig)
}

# Function to export analysis results
export_analysis_results <- function(data, filename_prefix = "ethanol_extract_analysis") {
 
  # Export detailed compound data
  export_data <- data %>%
    select(Peak, R.Time, Area, Area_Percent, Height, Height_Percent,
           A_H_Ratio, Compound_Name) %>%
    rename(
      `Peak Number` = Peak,
      `Retention Time (min)` = R.Time,
      `Peak Area` = Area,
      `Peak Area (%)` = Area_Percent,
      `Peak Height` = Height,
      `Peak Height (%)` = Height_Percent,
      `Area/Height Ratio` = A_H_Ratio,
      `Compound Name` = Compound_Name
    )
 
  write_csv(export_data, paste0(filename_prefix, "_detailed_results.csv"))
 
  # Export summary for major compounds
  major_compounds <- data[data$Area_Percent >= 5, ]
  if (nrow(major_compounds) > 0) {
    write_csv(major_compounds, paste0(filename_prefix, "_major_compounds.csv"))
  }
 
  cat("Results exported:\n")
  cat("- Detailed results:", paste0(filename_prefix, "_detailed_results.csv"), "\n")
  if (nrow(major_compounds) > 0) {
    cat("- Major compounds:", paste0(filename_prefix, "_major_compounds.csv"), "\n")
  }
}

# =====================================
# ANALYSIS EXECUTION
# =====================================

cat("=== ETHANOL EXTRACT GC-MS ANALYSIS ===\n")
cat("Loading and analyzing your GC-MS data...\n\n")

# Display summary statistics
create_summary_table(ethanol_extract_data)

# Create main chromatogram (by area)
cat("\nCreating chromatogram visualization...\n")
chromatogram_area <- create_gcms_chromatogram(
  ethanol_extract_data,
  plot_type = "area",
  show_labels = "major",
  label_threshold = 3.0,  # Show compounds ≥3% area
  title = "GC-MS Chromatogram - Ethanol Extract (Peak Area)"
)

print(chromatogram_area)

# Create chromatogram by height
chromatogram_height <- create_gcms_chromatogram(
  ethanol_extract_data,
  plot_type = "height",
  show_labels = "major",
  label_threshold = 3.0,
  title = "GC-MS Chromatogram - Ethanol Extract (Peak Height)"
)

# Create compound distribution pie chart
pie_chart <- create_compound_pie_chart(ethanol_extract_data, top_n = 8)
print(pie_chart)

# Create interactive visualization
cat("\nCreating interactive chromatogram...\n")
interactive_plot <- create_interactive_chromatogram(ethanol_extract_data, plot_type = "area")

# Display plots
cat("\nDisplaying visualizations...\n")
cat("1. Static chromatogram (area) - displayed above\n")
cat("2. Pie chart - displayed above\n")
cat("3. Interactive plot available in 'interactive_plot' object\n")

# Uncomment to show interactive plot
# print(interactive_plot)

# Chemical class analysis
cat("\n=== CHEMICAL CLASS ANALYSIS ===\n")
fatty_acids <- grep("acid", ethanol_extract_data$Compound_Name, ignore.case = TRUE)
esters <- grep("ester", ethanol_extract_data$Compound_Name, ignore.case = TRUE)
alcohols <- grep("ol$", ethanol_extract_data$Compound_Name, ignore.case = TRUE)
aldehydes <- grep("al$", ethanol_extract_data$Compound_Name, ignore.case = TRUE)
sterols <- grep("sterol", ethanol_extract_data$Compound_Name, ignore.case = TRUE)

cat("Fatty acids and derivatives:", length(fatty_acids), "compounds\n")
cat("Esters:", length(esters), "compounds\n")
cat("Alcohols:", length(alcohols), "compounds\n")
cat("Aldehydes:", length(aldehydes), "compounds\n")
cat("Sterols:", length(sterols), "compounds\n")

# Export results (uncomment to export)
# export_analysis_results(ethanol_extract_data)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Available objects:\n")
cat("- ethanol_extract_data: Your original data\n")
cat("- chromatogram_area: Static area chromatogram\n")
cat("- chromatogram_height: Static height chromatogram\n")
cat("- pie_chart: Compound distribution chart\n")
cat("- interactive_plot: Interactive chromatogram\n")
