# =============================================================================
# Daya Bay fish community multivariate analysis
# Final version
#
# Analyses:
#   1. log(x + 1) transformation
#   2. Bray-Curtis dissimilarity
#   3. NMDS
#      - k = 2
#      - trymax = 100
#      - all available samples
#   4. PERMANOVA
#      - Season as explanatory factor
#      - 9,999 permutations
#      - paired stations only
#      - permutations restricted within Station
#   5. PERMDISP
#      - homogeneity of multivariate dispersion
#      - 9,999 permutations
#   6. ANOSIM
#      - complementary rank-based analysis
#      - paired stations
#      - permutations restricted within Station
#
# Figure S2:
#   A Fish community (individuals): NMDS + PERMANOVA + PERMDISP
#   B Fish community (eggs):        NMDS + PERMANOVA + PERMDISP
#   C Fish community (individuals): ANOSIM
#   D Fish community (eggs):        ANOSIM
#
# Sampling design:
#
# Fish individuals
# Spring:
#   S1 S2 S3 S4 S5 S6 S8 S9 S10 S11 S13 S14
#
# Autumn:
#   S1 S2 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14
#
# Paired stations:
#   S1 S2 S4 S5 S6 S8 S9 S10 S11 S13 S14
#   n = 11 pairs
#
# Fish eggs:
#   Spring: S1-S14
#   Autumn: S1-S14
#   n = 14 pairs
#
# =============================================================================


# =============================================================================
# 1. Install and load packages
# =============================================================================

packages <- c(
  "vegan",
  "permute",
  "ggplot2",
  "patchwork"
)

new_packages <- packages[
  !(packages %in% installed.packages()[, "Package"])
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

library(vegan)
library(permute)
library(ggplot2)
library(patchwork)


# =============================================================================
# 2. Working directory
# =============================================================================
#
# Modify this path according to the actual location of your files.
#
# The two input files should be placed in this folder:
#
# 非度量多维标度（NMDS）分析 鱼.txt
# 非度量多维标度（NMDS）分析 鱼卵.txt
#
# =============================================================================

work_dir <- "C:/Users/Administrator/Desktop/NMDS_PERMANOVA"

if (!dir.exists(work_dir)) {
  stop(
    "Working directory does not exist:\n",
    work_dir
  )
}

setwd(work_dir)


# Output directory

output_dir <- file.path(
  work_dir,
  "Multivariate_results"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# =============================================================================
# 3. Input files
# =============================================================================

fish_file <- file.path(
  work_dir,
  "非度量多维标度（NMDS）分析 鱼.txt"
)

egg_file <- file.path(
  work_dir,
  "非度量多维标度（NMDS）分析 鱼卵.txt"
)


# Check files

cat("\n")
cat("============================================================\n")
cat("INPUT FILE CHECK\n")
cat("============================================================\n")

cat(
  "Fish file:\n",
  fish_file,
  "\nExists:",
  file.exists(fish_file),
  "\n\n"
)

cat(
  "Egg file:\n",
  egg_file,
  "\nExists:",
  file.exists(egg_file),
  "\n"
)


if (!file.exists(fish_file)) {
  stop(
    "Cannot find fish-individual file:\n",
    fish_file
  )
}

if (!file.exists(egg_file)) {
  stop(
    "Cannot find fish-egg file:\n",
    egg_file
  )
}


# =============================================================================
# 4. Function for reading community matrix
# =============================================================================

read_community <- function(file) {
  
  dat <- read.delim(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # First column = OTU/species ID
  taxa <- dat[[1]]
  
  # Remaining columns = samples
  mat <- dat[
    ,
    -1,
    drop = FALSE
  ]
  
  # Convert to numeric
  mat[] <- lapply(
    mat,
    function(x) {
      as.numeric(
        as.character(x)
      )
    }
  )
  
  # Transpose:
  # rows    = samples
  # columns = taxa
  mat <- t(
    as.matrix(mat)
  )
  
  colnames(mat) <- taxa
  
  
  # ---------------------------------------------------------------------------
  # Remove taxa with zero abundance across all samples
  # ---------------------------------------------------------------------------
  
  keep_taxa <- colSums(
    mat,
    na.rm = TRUE
  ) > 0
  
  mat <- mat[
    ,
    keep_taxa,
    drop = FALSE
  ]
  
  
  # ---------------------------------------------------------------------------
  # Check empty samples
  # ---------------------------------------------------------------------------
  
  sample_sum <- rowSums(
    mat,
    na.rm = TRUE
  )
  
  if (any(sample_sum == 0)) {
    
    empty_samples <- rownames(mat)[
      sample_sum == 0
    ]
    
    warning(
      "Empty samples detected and removed: ",
      paste(
        empty_samples,
        collapse = ", "
      )
    )
    
    mat <- mat[
      sample_sum > 0,
      ,
      drop = FALSE
    ]
  }
  
  
  return(mat)
}


# =============================================================================
# 5. Read data
# =============================================================================

fish <- read_community(
  fish_file
)

egg <- read_community(
  egg_file
)


cat("\n")
cat("============================================================\n")
cat("RAW DATA DIMENSIONS\n")
cat("============================================================\n")

cat(
  "Fish individuals:",
  nrow(fish),
  "samples x",
  ncol(fish),
  "taxa\n"
)

cat(
  "Fish eggs:",
  nrow(egg),
  "samples x",
  ncol(egg),
  "taxa\n"
)


# =============================================================================
# 6. Construct sample metadata
# =============================================================================


# -----------------------------------------------------------------------------
# 6.1 Fish individuals
# -----------------------------------------------------------------------------

fish_meta <- data.frame(
  
  SampleID = paste0(
    "Sample",
    1:25
  ),
  
  Season = c(
    rep("Spring", 12),
    rep("Autumn", 13)
  ),
  
  Station = c(
    
    # Spring
    # S7 and S12 were not sampled
    
    "S1",
    "S2",
    "S3",
    "S4",
    "S5",
    "S6",
    "S8",
    "S9",
    "S10",
    "S11",
    "S13",
    "S14",
    
    # Autumn
    # S3 was not sampled
    
    "S1",
    "S2",
    "S4",
    "S5",
    "S6",
    "S7",
    "S8",
    "S9",
    "S10",
    "S11",
    "S12",
    "S13",
    "S14"
  ),
  
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 6.2 Fish eggs
# -----------------------------------------------------------------------------

egg_meta <- data.frame(
  
  SampleID = paste0(
    "Sample",
    1:28
  ),
  
  Season = c(
    rep("Spring", 14),
    rep("Autumn", 14)
  ),
  
  Station = c(
    paste0("S", 1:14),
    paste0("S", 1:14)
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 7. Match metadata and community matrices
# =============================================================================

match_metadata <- function(
    comm,
    meta
) {
  
  idx <- match(
    rownames(comm),
    meta$SampleID
  )
  
  
  if (any(is.na(idx))) {
    
    stop(
      "The following samples were not found in metadata: ",
      paste(
        rownames(comm)[is.na(idx)],
        collapse = ", "
      )
    )
  }
  
  
  meta <- meta[
    idx,
    ,
    drop = FALSE
  ]
  
  
  if (!all(
    rownames(comm) ==
    meta$SampleID
  )) {
    
    stop(
      "Sample order mismatch."
    )
  }
  
  
  meta$Season <- factor(
    meta$Season,
    levels = c(
      "Spring",
      "Autumn"
    )
  )
  
  
  meta$Station <- factor(
    meta$Station
  )
  
  
  return(meta)
}


fish_meta <- match_metadata(
  fish,
  fish_meta
)

egg_meta <- match_metadata(
  egg,
  egg_meta
)


# =============================================================================
# 8. Check sampling design
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("FISH INDIVIDUAL SAMPLING DESIGN\n")
cat("============================================================\n")

print(
  table(
    fish_meta$Season
  )
)

print(
  table(
    fish_meta$Station,
    fish_meta$Season
  )
)


cat("\n")
cat("============================================================\n")
cat("FISH EGG SAMPLING DESIGN\n")
cat("============================================================\n")

print(
  table(
    egg_meta$Season
  )
)

print(
  table(
    egg_meta$Station,
    egg_meta$Season
  )
)


# =============================================================================
# 9. Main analysis function
# =============================================================================

run_multivariate_analysis <- function(
    comm,
    meta,
    community_name,
    output_dir
) {
  
  
  cat("\n\n")
  cat("############################################################\n")
  cat("# ", community_name, "\n")
  cat("############################################################\n")
  
  
  # ===========================================================================
  # A. log(x + 1) transformation
  # ===========================================================================
  
  comm_log <- log1p(
    comm
  )
  
  
  # ===========================================================================
  # B. Bray-Curtis dissimilarity
  # ===========================================================================
  
  bray_all <- vegdist(
    comm_log,
    method = "bray"
  )
  
  
  # ===========================================================================
  # C. NMDS
  #
  # NMDS uses all available samples.
  #
  # k = 2
  # trymax = 100
  # distance = Bray-Curtis
  # ===========================================================================
  
  set.seed(
    12345
  )
  
  
  nmds <- metaMDS(
    
    comm_log,
    
    distance = "bray",
    
    k = 2,
    
    trymax = 100,
    
    autotransform = FALSE,
    
    trace = FALSE
  )
  
  
  nmds_stress <- nmds$stress
  
  
  cat("\n")
  cat("---------------- NMDS ----------------\n")
  
  cat(
    "Stress =",
    nmds_stress,
    "\n"
  )
  
  
  # ===========================================================================
  # D. Extract NMDS coordinates
  # ===========================================================================
  
  nmds_scores <- as.data.frame(
    
    scores(
      nmds,
      display = "sites"
    )
  )
  
  
  nmds_scores$SampleID <- rownames(
    nmds_scores
  )
  
  
  # Keep the original sample order
  
  nmds_scores <- nmds_scores[
    match(
      meta$SampleID,
      nmds_scores$SampleID
    ),
    ,
    drop = FALSE
  ]
  
  
  nmds_scores$Season <- meta$Season
  
  nmds_scores$Station <- meta$Station
  
  
  write.csv(
    
    nmds_scores,
    
    file.path(
      output_dir,
      paste0(
        community_name,
        "_NMDS_scores.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  # ===========================================================================
  # E. Identify paired stations
  #
  # Only stations sampled in both Spring and Autumn are retained for
  # PERMANOVA, PERMDISP and ANOSIM.
  # ===========================================================================
  
  station_season <- table(
    meta$Station,
    meta$Season
  )
  
  
  paired_stations <- rownames(
    station_season
  )[
    station_season[, "Spring"] == 1 &
      station_season[, "Autumn"] == 1
  ]
  
  
  cat("\n")
  cat("---------------- PAIRED DESIGN ----------------\n")
  
  cat(
    "Paired stations:\n",
    paste(
      paired_stations,
      collapse = ", "
    ),
    "\n"
  )
  
  cat(
    "Number of paired stations =",
    length(
      paired_stations
    ),
    "\n"
  )
  
  
  # ===========================================================================
  # F. Extract paired samples
  # ===========================================================================
  
  keep <- as.character(
    meta$Station
  ) %in% paired_stations
  
  
  comm_pair <- comm_log[
    keep,
    ,
    drop = FALSE
  ]
  
  
  meta_pair <- droplevels(
    meta[
      keep,
      ,
      drop = FALSE
    ]
  )
  
  
  # ===========================================================================
  # G. Check paired structure
  # ===========================================================================
  
  pair_check <- table(
    meta_pair$Station,
    meta_pair$Season
  )
  
  
  cat("\nPaired sample structure:\n")
  
  print(
    pair_check
  )
  
  
  if (
    any(
      pair_check[, "Spring"] != 1
    ) ||
    any(
      pair_check[, "Autumn"] != 1
    )
  ) {
    
    stop(
      "Paired design error: each paired station must contain ",
      "exactly one Spring and one Autumn observation."
    )
  }
  
  
  cat(
    "\nSamples used for paired analyses =",
    nrow(
      comm_pair
    ),
    "\n"
  )
  
  
  # ===========================================================================
  # H. Bray-Curtis for paired samples
  # ===========================================================================
  
  bray_pair <- vegdist(
    comm_pair,
    method = "bray"
  )
  
  
  # ===========================================================================
  # I. Restricted permutation design
  #
  # Spring and Autumn labels are only permuted within each sampling station.
  # ===========================================================================
  
  perm_design <- how(
    
    nperm = 9999,
    
    blocks = meta_pair$Station
  )
  
  
  # ===========================================================================
  # J. PERMANOVA
  #
  # Outputs:
  #   pseudo-F
  #   R2
  #   permutation-based P
  # ===========================================================================
  
  set.seed(
    12345
  )
  
  
  permanova <- adonis2(
    
    bray_pair ~ Season,
    
    data = meta_pair,
    
    permutations = perm_design,
    
    by = NULL
  )
  
  
  cat("\n")
  cat("---------------- PERMANOVA ----------------\n")
  
  print(
    permanova
  )
  
  
  permanova_F <- as.numeric(
    permanova[
      1,
      "F"
    ]
  )
  
  
  permanova_R2 <- as.numeric(
    permanova[
      1,
      "R2"
    ]
  )
  
  
  permanova_P <- as.numeric(
    permanova[
      1,
      "Pr(>F)"
    ]
  )
  
  
  # ===========================================================================
  # K. PERMDISP
  #
  # Tests homogeneity of multivariate dispersion between seasons.
  #
  # Outputs:
  #   F
  #   permutation-based P
  # ===========================================================================
  
  dispersion <- betadisper(
    
    bray_pair,
    
    group = meta_pair$Season,
    
    type = "centroid",
    
    bias.adjust = TRUE
  )
  
  
  set.seed(
    12345
  )
  
  
  permdisp_perm <- permutest(
    
    dispersion,
    
    permutations = 9999
  )
  
  
  cat("\n")
  cat("---------------- PERMDISP ----------------\n")
  
  print(
    permdisp_perm
  )
  
  
  permdisp_table <- permdisp_perm$tab
  
  
  permdisp_F <- as.numeric(
    permdisp_table[
      1,
      "F"
    ]
  )
  
  
  permdisp_P <- as.numeric(
    permdisp_table[
      1,
      "Pr(>F)"
    ]
  )
  
  
  # ===========================================================================
  # L. ANOSIM
  #
  # Complementary rank-based analysis.
  # Uses the same paired subset.
  # Permutations are restricted by sampling station.
  # ===========================================================================
  
  set.seed(
    12345
  )
  
  
  anosim_result <- anosim(
    
    bray_pair,
    
    grouping = meta_pair$Season,
    
    permutations = 9999,
    
    strata = meta_pair$Station
  )
  
  
  cat("\n")
  cat("---------------- ANOSIM ----------------\n")
  
  print(
    anosim_result
  )
  
  
  anosim_R <- as.numeric(
    anosim_result$statistic
  )
  
  
  anosim_P <- as.numeric(
    anosim_result$signif
  )
  
  
  # ===========================================================================
  # M. Summary table
  # ===========================================================================
  
  summary_table <- data.frame(
    
    Community = community_name,
    
    Total_samples_NMDS = nrow(
      comm_log
    ),
    
    Paired_stations = length(
      paired_stations
    ),
    
    Samples_paired_analysis = nrow(
      comm_pair
    ),
    
    NMDS_dimensions = 2,
    
    NMDS_random_starts = 100,
    
    NMDS_stress = nmds_stress,
    
    PERMANOVA_pseudo_F = permanova_F,
    
    PERMANOVA_R2 = permanova_R2,
    
    PERMANOVA_P = permanova_P,
    
    PERMDISP_F = permdisp_F,
    
    PERMDISP_P = permdisp_P,
    
    ANOSIM_R = anosim_R,
    
    ANOSIM_P = anosim_P
  )
  
  
  # ===========================================================================
  # N. Export results
  # ===========================================================================
  
  write.csv(
    
    summary_table,
    
    file.path(
      output_dir,
      paste0(
        community_name,
        "_statistical_summary.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  capture.output(
    
    permanova,
    
    file = file.path(
      output_dir,
      paste0(
        community_name,
        "_PERMANOVA.txt"
      )
    )
  )
  
  
  capture.output(
    
    permdisp_perm,
    
    file = file.path(
      output_dir,
      paste0(
        community_name,
        "_PERMDISP.txt"
      )
    )
  )
  
  
  capture.output(
    
    anosim_result,
    
    file = file.path(
      output_dir,
      paste0(
        community_name,
        "_ANOSIM.txt"
      )
    )
  )
  
  
  # ===========================================================================
  # O. Return results
  # ===========================================================================
  
  return(
    
    list(
      
      NMDS = nmds,
      
      NMDS_scores = nmds_scores,
      
      NMDS_stress = nmds_stress,
      
      PERMANOVA = permanova,
      
      PERMANOVA_F = permanova_F,
      
      PERMANOVA_R2 = permanova_R2,
      
      PERMANOVA_P = permanova_P,
      
      PERMDISP = permdisp_perm,
      
      PERMDISP_F = permdisp_F,
      
      PERMDISP_P = permdisp_P,
      
      ANOSIM = anosim_result,
      
      ANOSIM_R = anosim_R,
      
      ANOSIM_P = anosim_P,
      
      Bray_pair = bray_pair,
      
      Meta_pair = meta_pair,
      
      summary = summary_table,
      
      paired_stations = paired_stations
    )
  )
}


# =============================================================================
# 10. Run fish-individual analysis
# =============================================================================

fish_result <- run_multivariate_analysis(
  
  comm = fish,
  
  meta = fish_meta,
  
  community_name = "Fish_individuals",
  
  output_dir = output_dir
)


# =============================================================================
# 11. Run fish-egg analysis
# =============================================================================

egg_result <- run_multivariate_analysis(
  
  comm = egg,
  
  meta = egg_meta,
  
  community_name = "Fish_eggs",
  
  output_dir = output_dir
)


# =============================================================================
# 12. Combine final statistical results
# =============================================================================

final_summary <- rbind(
  
  fish_result$summary,
  
  egg_result$summary
)


write.csv(
  
  final_summary,
  
  file.path(
    output_dir,
    "FINAL_multivariate_statistical_summary.csv"
  ),
  
  row.names = FALSE
)


# =============================================================================
# 13. Print manuscript-ready statistics
# =============================================================================

cat("\n\n")
cat("============================================================\n")
cat("MANUSCRIPT-READY STATISTICS\n")
cat("============================================================\n")


for (
  i in seq_len(
    nrow(final_summary)
  )
) {
  
  x <- final_summary[
    i,
    ,
    drop = FALSE
  ]
  
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  
  cat(
    x$Community,
    "\n"
  )
  
  cat("------------------------------------------------------------\n")
  
  
  cat(
    
    sprintf(
      "NMDS: Stress = %.3f\n",
      x$NMDS_stress
    )
  )
  
  
  cat(
    
    sprintf(
      paste0(
        "PERMANOVA: pseudo-F = %.3f, ",
        "R2 = %.3f, ",
        "P = %.4f\n"
      ),
      x$PERMANOVA_pseudo_F,
      x$PERMANOVA_R2,
      x$PERMANOVA_P
    )
  )
  
  
  cat(
    
    sprintf(
      "PERMDISP: F = %.3f, P = %.4f\n",
      x$PERMDISP_F,
      x$PERMDISP_P
    )
  )
  
  
  cat(
    
    sprintf(
      "ANOSIM: R = %.3f, P = %.4f\n",
      x$ANOSIM_R,
      x$ANOSIM_P
    )
  )
}


# =============================================================================
# 14. Function for formatting P values
# =============================================================================

format_p <- function(p) {
  
  if (is.na(p)) {
    
    return("= NA")
  }
  
  
  if (p < 0.001) {
    
    return("< 0.001")
  }
  
  
  return(
    paste0(
      "= ",
      sprintf(
        "%.3f",
        p
      )
    )
  )
}


# =============================================================================
# 15. Prepare ANOSIM plotting data
# =============================================================================

prepare_anosim_plot_data <- function(
    bray,
    meta
) {
  
  dmat <- as.matrix(
    bray
  )
  
  
  n <- nrow(
    dmat
  )
  
  
  result <- data.frame(
    
    Dissimilarity = numeric(0),
    
    Category = character(0),
    
    stringsAsFactors = FALSE
  )
  
  
  for (
    i in 1:(n - 1)
  ) {
    
    for (
      j in (i + 1):n
    ) {
      
      
      season_i <- as.character(
        meta$Season[i]
      )
      
      
      season_j <- as.character(
        meta$Season[j]
      )
      
      
      # -----------------------------------------------------------------------
      # Between seasons
      # -----------------------------------------------------------------------
      
      if (
        season_i != season_j
      ) {
        
        category <- "Between"
        
        
        # -----------------------------------------------------------------------
        # Within Autumn
        # -----------------------------------------------------------------------
        
      } else if (
        season_i == "Autumn" &&
        season_j == "Autumn"
      ) {
        
        category <- "Autumn"
        
        
        # -----------------------------------------------------------------------
        # Within Spring
        # -----------------------------------------------------------------------
        
      } else {
        
        category <- "Spring"
      }
      
      
      result <- rbind(
        
        result,
        
        data.frame(
          
          Dissimilarity = dmat[
            i,
            j
          ],
          
          Category = category,
          
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # Rank all Bray-Curtis dissimilarities
  # ---------------------------------------------------------------------------
  
  result$Rank <- rank(
    
    result$Dissimilarity,
    
    ties.method = "average"
  )
  
  
  result$Category <- factor(
    
    result$Category,
    
    levels = c(
      "Between",
      "Autumn",
      "Spring"
    )
  )
  
  
  return(
    result
  )
}


# =============================================================================
# 16. Generate ANOSIM plotting data
# =============================================================================

fish_anosim_data <- prepare_anosim_plot_data(
  
  fish_result$Bray_pair,
  
  fish_result$Meta_pair
)


egg_anosim_data <- prepare_anosim_plot_data(
  
  egg_result$Bray_pair,
  
  egg_result$Meta_pair
)


# =============================================================================
# 17. NMDS plotting function
# =============================================================================

make_nmds_panel <- function(
    result,
    panel_letter,
    panel_title,
    show_y_title = TRUE,
    show_legend = FALSE
) {
  
  
  df <- result$NMDS_scores
  
  
  # ---------------------------------------------------------------------------
  # Statistical annotation line 1:
  # Stress + PERMANOVA
  # ---------------------------------------------------------------------------
  
  stat_line1 <- paste0(
    
    "Stress = ",
    
    sprintf(
      "%.3f",
      result$NMDS_stress
    ),
    
    "      ",
    
    "PERMANOVA: pseudo-F = ",
    
    sprintf(
      "%.3f",
      result$PERMANOVA_F
    ),
    
    ", R² = ",
    
    sprintf(
      "%.3f",
      result$PERMANOVA_R2
    ),
    
    ", P ",
    
    format_p(
      result$PERMANOVA_P
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Statistical annotation line 2:
  # PERMDISP
  # ---------------------------------------------------------------------------
  
  stat_line2 <- paste0(
    
    "PERMDISP: F = ",
    
    sprintf(
      "%.3f",
      result$PERMDISP_F
    ),
    
    ", P ",
    
    format_p(
      result$PERMDISP_P
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  
  p <- ggplot(
    
    df,
    
    aes(
      x = NMDS1,
      y = NMDS2,
      colour = Season,
      shape = Season
    )
  ) +
    
    
    # 95% confidence ellipses
    
    stat_ellipse(
      
      aes(
        group = Season
      ),
      
      type = "norm",
      
      level = 0.95,
      
      linewidth = 0.75,
      
      show.legend = FALSE
    ) +
    
    
    # Sample points
    
    geom_point(
      size = 2.8
    ) +
    
    
    # Point shapes
    
    scale_shape_manual(
      
      values = c(
        Spring = 16,
        Autumn = 17
      )
    ) +
    
    
    labs(
      
      title = paste0(
        panel_letter,
        "    ",
        panel_title
      ),
      
      subtitle = paste0(
        stat_line1,
        "\n",
        stat_line2
      ),
      
      x = "NMDS1",
      
      y = if (
        show_y_title
      ) {
        "NMDS2"
      } else {
        NULL
      },
      
      colour = NULL,
      
      shape = NULL
    ) +
    
    
    theme_classic(
      base_size = 11
    ) +
    
    
    theme(
      
      plot.title = element_text(
        hjust = 0.5,
        size = 12,
        face = "plain"
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 8.2,
        lineheight = 1.15
      ),
      
      axis.title = element_text(
        size = 11
      ),
      
      axis.text = element_text(
        size = 9
      ),
      
      legend.position = if (
        show_legend
      ) {
        "right"
      } else {
        "none"
      },
      
      legend.title = element_blank(),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
  
  
  return(
    p
  )
}


# =============================================================================
# 18. ANOSIM plotting function
# =============================================================================

make_anosim_panel <- function(
    data,
    result,
    panel_letter,
    panel_title,
    show_y_title = TRUE
) {
  
  
  stat_text <- paste0(
    
    "ANOSIM: R = ",
    
    sprintf(
      "%.3f",
      result$ANOSIM_R
    ),
    
    "        P ",
    
    format_p(
      result$ANOSIM_P
    )
  )
  
  
  p <- ggplot(
    
    data,
    
    aes(
      x = Category,
      y = Rank,
      fill = Category
    )
  ) +
    
    
    geom_boxplot(
      
      width = 0.70,
      
      outlier.shape = NA,
      
      linewidth = 0.55
    ) +
    
    
    labs(
      
      title = paste0(
        panel_letter,
        "    ",
        panel_title
      ),
      
      subtitle = stat_text,
      
      x = NULL,
      
      y = if (
        show_y_title
      ) {
        "Rank of dissimilarity"
      } else {
        NULL
      }
    ) +
    
    
    theme_classic(
      base_size = 11
    ) +
    
    
    theme(
      
      plot.title = element_text(
        hjust = 0.5,
        size = 12,
        face = "plain"
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 9
      ),
      
      axis.title = element_text(
        size = 11
      ),
      
      axis.text = element_text(
        size = 9
      ),
      
      legend.position = "none",
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
  
  
  return(
    p
  )
}


# =============================================================================
# 19. Generate Figure S2 panels
# =============================================================================


# -----------------------------------------------------------------------------
# Panel A
# Fish individuals
# NMDS + PERMANOVA + PERMDISP
# -----------------------------------------------------------------------------

plot_A <- make_nmds_panel(
  
  result = fish_result,
  
  panel_letter = "A",
  
  panel_title = "Fish community (individuals)",
  
  show_y_title = TRUE,
  
  show_legend = FALSE
)


# -----------------------------------------------------------------------------
# Panel B
# Fish eggs
# NMDS + PERMANOVA + PERMDISP
# -----------------------------------------------------------------------------

plot_B <- make_nmds_panel(
  
  result = egg_result,
  
  panel_letter = "B",
  
  panel_title = "Fish community (eggs)",
  
  show_y_title = FALSE,
  
  show_legend = TRUE
)


# -----------------------------------------------------------------------------
# Panel C
# Fish individuals
# ANOSIM
# -----------------------------------------------------------------------------

plot_C <- make_anosim_panel(
  
  data = fish_anosim_data,
  
  result = fish_result,
  
  panel_letter = "C",
  
  panel_title = "Fish community (individuals)",
  
  show_y_title = TRUE
)


# -----------------------------------------------------------------------------
# Panel D
# Fish eggs
# ANOSIM
# -----------------------------------------------------------------------------

plot_D <- make_anosim_panel(
  
  data = egg_anosim_data,
  
  result = egg_result,
  
  panel_letter = "D",
  
  panel_title = "Fish community (eggs)",
  
  show_y_title = FALSE
)


# =============================================================================
# 20. Combine Figure S2
# =============================================================================

figure_S2 <- (
  
  plot_A + plot_B
  
) / (
  
  plot_C + plot_D
  
) +
  
  plot_layout(
    
    widths = c(
      1,
      1
    ),
    
    heights = c(
      1,
      1.10
    ),
    
    guides = "collect"
  )


# =============================================================================
# 21. Display Figure S2
# =============================================================================

print(
  figure_S2
)


# =============================================================================
# 22. Export Figure S2 as PDF
# =============================================================================

ggsave(
  
  filename = file.path(
    output_dir,
    "Figure_S2_NMDS_PERMANOVA_PERMDISP_ANOSIM.pdf"
  ),
  
  plot = figure_S2,
  
  width = 11,
  
  height = 11,
  
  units = "in",
  
  device = cairo_pdf
)


# =============================================================================
# 23. Export Figure S2 as 600-dpi PNG
# =============================================================================

ggsave(
  
  filename = file.path(
    output_dir,
    "Figure_S2_NMDS_PERMANOVA_PERMDISP_ANOSIM.png"
  ),
  
  plot = figure_S2,
  
  width = 11,
  
  height = 11,
  
  units = "in",
  
  dpi = 600
)


# =============================================================================
# 24. Export Figure S2 statistics
# =============================================================================

figure_S2_statistics <- data.frame(
  
  Community = c(
    "Fish community (individuals)",
    "Fish community (eggs)"
  ),
  
  NMDS_total_samples = c(
    nrow(fish),
    nrow(egg)
  ),
  
  Paired_stations = c(
    length(
      fish_result$paired_stations
    ),
    length(
      egg_result$paired_stations
    )
  ),
  
  NMDS_Stress = c(
    fish_result$NMDS_stress,
    egg_result$NMDS_stress
  ),
  
  PERMANOVA_pseudo_F = c(
    fish_result$PERMANOVA_F,
    egg_result$PERMANOVA_F
  ),
  
  PERMANOVA_R2 = c(
    fish_result$PERMANOVA_R2,
    egg_result$PERMANOVA_R2
  ),
  
  PERMANOVA_P = c(
    fish_result$PERMANOVA_P,
    egg_result$PERMANOVA_P
  ),
  
  PERMDISP_F = c(
    fish_result$PERMDISP_F,
    egg_result$PERMDISP_F
  ),
  
  PERMDISP_P = c(
    fish_result$PERMDISP_P,
    egg_result$PERMDISP_P
  ),
  
  ANOSIM_R = c(
    fish_result$ANOSIM_R,
    egg_result$ANOSIM_R
  ),
  
  ANOSIM_P = c(
    fish_result$ANOSIM_P,
    egg_result$ANOSIM_P
  )
)


write.csv(
  
  figure_S2_statistics,
  
  file.path(
    output_dir,
    "Figure_S2_statistical_results.csv"
  ),
  
  row.names = FALSE
)


# =============================================================================
# 25. Print Figure S2 statistical results
# =============================================================================

cat("\n\n")
cat("============================================================\n")
cat("FIGURE S2 STATISTICAL RESULTS\n")
cat("============================================================\n")

print(
  figure_S2_statistics
)


# =============================================================================
# 26. Interpretation of PERMDISP
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("PERMDISP INTERPRETATION\n")
cat("============================================================\n")


# Fish individuals

cat("\nFish individuals:\n")

if (
  fish_result$PERMDISP_P < 0.05
) {
  
  cat(
    "PERMDISP is significant (P < 0.05).\n",
    "Seasonal differences detected by PERMANOVA should be interpreted ",
    "with caution because multivariate dispersion also differs between seasons.\n"
  )
  
} else {
  
  cat(
    "PERMDISP is not significant (P >= 0.05).\n",
    "There is no evidence of a significant difference in multivariate ",
    "dispersion between seasons.\n"
  )
}


# Fish eggs

cat("\nFish eggs:\n")

if (
  egg_result$PERMDISP_P < 0.05
) {
  
  cat(
    "PERMDISP is significant (P < 0.05).\n",
    "Seasonal differences detected by PERMANOVA should be interpreted ",
    "with caution because multivariate dispersion also differs between seasons.\n"
  )
  
} else {
  
  cat(
    "PERMDISP is not significant (P >= 0.05).\n",
    "There is no evidence of a significant difference in multivariate ",
    "dispersion between seasons.\n"
  )
}


# =============================================================================
# 27. Save R session information
# =============================================================================

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    output_dir,
    "R_session_information.txt"
  )
)


# =============================================================================
# 28. Finish
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat(
  "All results were saved to:\n",
  output_dir,
  "\n"
)

cat("\nMain output files:\n")

cat(
  "1. FINAL_multivariate_statistical_summary.csv\n"
)

cat(
  "2. Figure_S2_statistical_results.csv\n"
)

cat(
  "3. Figure_S2_NMDS_PERMANOVA_PERMDISP_ANOSIM.pdf\n"
)

cat(
  "4. Figure_S2_NMDS_PERMANOVA_PERMDISP_ANOSIM.png\n"
)

cat("============================================================\n")