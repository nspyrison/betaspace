# GUIDING WORKFLOW -----

## {Phyloseq} --> {vegan} dist/ord ----

#' @examples
#' library(phyloseq)
#' library(vegan)
#' library(dplyr)
#' data("GlobalPatterns")
#' 
#' # 1. Aggregate counts at the Class level
#' gp_class <- tax_glom(GlobalPatterns, taxrank = "Genus")
#' # 2. Extract the OTU table (samples x taxa) and transpose for vegan
#' otu_class <- as.data.frame(t(otu_table(gp_class)))
#' # 3. Extract the Class labels for each sample (from the taxonomic table)
#' #    We'll use the most abundant Class in each sample as a simple proxy.
#' tax_class <- tax_table(gp_class)[, "Class", drop = FALSE]
#' # For this example, let's assign each sample to its most abundant Class
#' sample_classes <- apply(otu_class, 1, function(x) {
#'   idx <- which.max(x)
#'   tax_class[idx, "Class"]
#' })
#' sample_classes <- as.factor(sample_classes)
#' # 4. Remove Classes with very low total abundance (e.g., < 10 counts total)
#' class_sums <- colSums(otu_class)
#' otu_class_filt <- otu_class[, class_sums >= 10]
#' 
#' 
#' # Compute Bray-Curtis distance and perform CAP with Class constraint
#' bray_dist <- vegdist(otu_class_filt, method = "bray")
#' 
#' # CAP: constrained by the sample's dominant Class
#' cap_result <- capscale(bray_dist ~ sample_classes)
#' 
#' # 6. Extract the constrained axes (CA1, CA2) scores
#' cap_scores <- scores(cap_result, display = "sites", choices = 1:2)
#' #cap_scores <- scores(cap_result, display = "all", choices = 1:20)
#' cap_scores_df <- as.data.frame(cap_scores)
#' #colnames(cap_scores_df) <- c("CA1", "CA2")
#' cap_scores_df$Class <- sample_classes
#' 
#' # # 7. Examine the proportion of variance explained by the constraint
#' # summary(cap_result)$concont$importance
#' 
#' library(ggplot2)
#' 
#' # Plot constrained ordination with ellipses
#' p <- ggplot(cap_scores_df, aes(x = CAP1, y = CAP2, color = Class)) +
#'   geom_point(size = 3, alpha = 0.8) +
#'   stat_ellipse(level = 0.95, linewidth = 0.8, alpha = 0.5) +
#'   labs(
#'     title = "CAP of GlobalPatterns (Bray-Curtis, constrained by dominant Class)",
#'     subtitle = paste("Constrained variance:", 
#'                      round(summary(cap_result)$concont$importance[2, 1] * 100, 1), "%"),
#'     x = paste0("CA1 (", 
#'                round(summary(cap_result)$concont$importance[2, 1] * 100, 1), "% constrained)"),
#'     y = paste0("CA2 (", 
#'                round(summary(cap_result)$concont$importance[2, 2] * 100, 1), "% constrained)")
#'   ) +
#'   theme_minimal() +
#'   theme(legend.position = "right")
#' 
#' print(p)


## {Phyloseq} thru CAP ----
#' @examples
#' library(phyloseq)
#' library(vegan)
#' library(dplyr)
#' library(ggplot2)
#' 
#' 
#' # 1. Aggregate counts at the Class level
#' gp_class <- tax_glom(GlobalPatterns, taxrank = "Genus")
#' 
#' # Manually compute dbrda (which is what ordinate() should do)
#' otu_table_mat <- as(otu_table(gp_class), "matrix")
#' bray_dist <- vegdist(t(otu_table_mat), method = "bray")
#' class_var <- sample_data(gp_class)$Class
#' 
#' # This gives you what ordinate(..., method="dbrda", ...) would return
#' gp_bray_cap0 <- dbrda(bray_dist ~ class_var)
#' gp_bray_cap <- ordinate(gp_class, method = "CCA", distance = "bray", formula = ~ Class)



### ORIG -----
# You now have a proper vegan dbrda object that your utilities expect

# # 2. Derive the sample class proxy (most abundant Class)
# # We extract matrices to perform the calculation
# otu_mat <- as(otu_table(gp_class), "matrix")
# tax_mat <- as(tax_table(gp_class), "matrix")
# 
# sample_classes <- apply(otu_mat, 1, function(x) {
#   idx <- which.max(x)
#   tax_mat[idx, "Class"]
# })
# 
# # 3. Add the proxy to sample_data so ordinate() can find it via formula
# # This is the crucial step for the phyloseq workflow
# df_metadata <- data.frame(
#   sample_id = sample_names(gp_class),
#   Class = sample_data(gp_class)$SampleType,
#   row.names = sample_names(gp_class)
# )
# sample_data(gp_class) <- sample_data(df_metadata)
# 
# # 4. Perform constrained ordination using phyloseq::ordinate()
# # method = "dbrda" is the distance-based RDA (equivalent to capscale)
# gp_bray_cap <- ordinate(gp_class, method = "dbrda", distance = "bray", formula = ~ Class)


### ALT ----
# library(vegan)
# 
# # Get OTU table and convert to numeric matrix
# otu_table_mat <- as(otu_table(gp_class), "matrix")
# 
# # Calculate Bray-Curtis distance on samples (transpose so samples are rows)
# bray_dist <- vegdist(t(otu_table_mat), method = "bray")
# 
# # Extract Class variable
# class_var <- sample_data(gp_class)$Class
# 
# # Perform dbrda :(
# gp_bray_cap <- vegan::dbrda(bray_dist ~ class_var)
# 
# # 5. Extract scores and visualize
# # The scores() function works on the object returned by ordinate()
# scores_df <- as.data.frame(scores(gp_bray_cap, display = "sites"))
# scores_df$Class <- sample_data(gp_class)$Class
# 
# ggplot(scores_df, aes(x = dbRDA1, y = dbRDA2, color = Class)) +
#   geom_point(size = 3, alpha = 0.8) +
#   stat_ellipse(level = 0.95, linewidth = 0.8, alpha = 0.5) +
#   theme_minimal() +
#   labs(
#     title = "Phyloseq ordinate (dbrda) of GlobalPatterns",
#     subtitle = "Constrained by most abundant Class"
#   )


## {Vegan} binned continuous classes, no taxa ----
#' @examples 
#' library(vegan)
#' data(mite, mite.env)
#' 
#' mite.env$SubsDens %>% hist
#' # Create a categorical constraint from substrate density
#' mite.env$SubsDens_cat <- cut(
#'   mite.env$SubsDens, breaks = 3, labels = c("Low", "Medium", "High"))
#' mite.env$SubsDens_cat %>% table()
#' 
#' # Bray-Curtis distance
#' mite_bray <- vegdist(mite, method = "bray")
#' 
#' # CAP constrained by substrate density category
#' mite_cap <- capscale(mite_bray ~ SubsDens_cat, data = mite.env)
#' 
#' # Permutation test (more reliable with n=70)
#' anova(mite_cap, permutations = 999)
#' 
#' # Extract site scores
#' cap_scores <- scores(mite_cap, display = "sites", choices = c(1, 2))
#' cap_df <- as.data.frame(cap_scores)
#' colnames(cap_df) <- c("CA1", "CA2")
#' cap_df$SubsDens <- mite.env$SubsDens_cat
#' 
#' # Variance explained by constrained axes
#' var_exp <- summary(mite_cap)$concont$importance
#' cat("Constrained variance:", round(var_exp[2, 1] * 100, 1), "%\n")
#' 
#' # Simple scatter plot
#' ggplot(cap_df, aes(x = CA1, y = CA2, color = SubsDens)) +
#'   geom_point(size = 3) +
#'   labs(title = "CAP of mite communities (constrained by substrate density)",
#'        x = paste0("CA1 (", round(var_exp[2, 1] * 100, 1), "% constrained)"),
#'        y = paste0("CA2 (", round(var_exp[2, 2] * 100, 1), "% constrained)")) +
#'   theme_minimal()


## EMP data -- not working  ------
#' @examples
#' # Install if needed
#' if (!require("BiocManager", quietly = TRUE))
#'   install.packages("BiocManager")
#' if (!require("curatedMetagenomicData", quietly = TRUE))
#'   BiocManager::install("curatedMetagenomicData")
#' 
#' library(curatedMetagenomicData)
#' library(phyloseq)
#' 
#' # Load a specific study (e.g., "NielsenHB_2014" - Danish gut microbiome)
#' tictoc::tic("start")
#' emp_data <- curatedMetagenomicData("NielsenHB_2014.metaphlan_bugs_list.stool", dryrun = FALSE)
#' ps_emp <- emp_data[[1]]  # Extract the phyloseq object
#' 
#' # Check dimensions
#' ps_emp
#' # otu_table()   OTU Table:         [ 552 taxa and 396 samples ]
#' # tax_table()   Taxonomy Table:    [ 552 taxa by 7 taxonomic ranks ]


# Bigger example scripts -----

## {vegan} Bray, PCoA, GGally -----
#' @examples
#' library(vegan)
#' data(dune) ## Numeric X
#' data(dune.env) ## Site metadata, like class ($Management)
#' dist_type <- "bray" # Distance method
#' clas <- dune.env$Management
#' 
#' ## CAP, Constrained Analysis of Principal coordinates, (Constrained)
#' dune_bray_cap <- capscale(dune ~ clas, distance = dist_type)
#' # dune_bray_cap0 <- capscale(dune ~ Management, data = dune.env, distance = dist_type)
#' # (dune_bray_cap$Ybar == dune_bray_cap0$Ybar) %>% table ## TRUE 280
#'
#' 
#' dune_bray_cap_scree <- scree_df(.ordination = dune_bray_cap)
#' debugonce(scree_plot)
#' scree_plot(dune_bray_cap_scree)
#' 
#' # Estimate number of dimensions to embed -- ON RAW DATA
#' (dune_ide <- ide_math(dune, inc_slow = FALSE))
#' (accepted_d_raw <- dune_ide %>% median(na.rm = TRUE) %>% floor())
#' (diff_d_embed_dist_less_raw <-
#'   dune_bray_cap_scree$meta_df$d_embedded - accepted_d_raw)
#' .d <- dune_bray_cap_scree$meta_df$d_embedded
#' 
#' # Visualized accepted embedding space
#' # Extract site coordinates on the constrained (supervised) axes
#' site_scores <- scores(dune_bray_cap, display = "sites", choices = 1:.d)
#' 
#' # Extract taxon/variable scores
#' species_scores <- scores(dune_bray_cap, display = "species", choices = 1:.d)
#' 

#' 
#' pairs(site_scores)
#' ## Color on class
#' pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
#' ## RColorBrewer::brewer.pal(n = 4, name = "Dark2") %>% dput()
#' col_vec <- pal[as.numeric(clas)]
#' pairs(site_scores, main = "Site Scores -- Dune Bray CAP",
#'       pch = 21, bg = col_vec)
#' GGally::ggpairs(site_scores, ggplot2::aes(color = clas))


## Manual GGally pairs with ellipse -----
#' @examples 
#' library(ggplot2)
#' library(patchwork)
#' 
#' {
#'   tictoc::tic("MANUAL GG PAIRS")
#'   # Define the variables you want to compare
#'   
#'   vars <- colnames(site_scores)
#'   n <- length(vars)
#'   
#'   # Create a list to store plots
#'   plot_list <- list()
#'   
#'   for (i in 1:n) {
#'     for (j in 1:n) {
#'       # Calculate index for a 1D list to use with patchwork
#'       idx <- (i - 1) * n + j
#'       
#'       if (i == j) {
#'         # Diagonal: Density plots or similar
#'         plot_list[[idx]] <- ggplot(df_plot, aes(x = .data[[vars[i]]], fill = clas)) +
#'           geom_density(alpha = 0.5) +
#'           theme_minimal() +
#'           theme(legend.position = "none")
#'         
#'       } else if (i > j) {
#'         # Lower triangle: Scatter + Ellipse
#'         plot_list[[idx]] <- ggplot(df_plot, aes(x = .data[[vars[j]]], y = .data[[vars[i]]], color = clas)) +
#'           geom_point(size = 1) +
#'           stat_ellipse() +
#'           #theme_minimal() +
#'           theme(legend.position = "none")
#'         
#'       } else {
#'         # Upper triangle: Correlation or just empty
#'         plot_list[[idx]] <- ggplot() + theme_void()
#'       }
#'     }
#'   }
#'   
#'   # Wrap the list into a grid using patchwork
#'   wrap_plots(plot_list, ncol = n) %>% print()
#'   tictoc::toc()
#' }
#' {
#'   tictoc::tic("GGALLY GG PAIRS")
#'   GGally::ggpairs(site_scores, ggplot2::aes(color = clas)) %>% print()
#'   tictoc::toc()
#' }







## Different coordinate embedding and fortification -----
#' @examples
#' library(vegan)
#' data(dune) ## Numeric X
#' data(dune.env) ## Site metadata, like class ($Management)
#' dist_type <- "bray" # Distance method
#' 
#' # Embedding distance matrix
#' ## PCoA, Principal Coordinates Analysis, (Unconstrained) on precomputed dist
#' dune_bray <- vegan::vegdist(dune, method = "dist_type)
#' dune_bray_pcoa <- stats::cmdscale(dune_bray, k = nrow(dune) - 1, eig = TRUE)
#' ## PCO, Principal COordinates analysis, aka PCoA, the {vegan} format
#' dune_bray_pco <- vegan::pco(dune, dist = dist_type, sqrt.dist = TRUE)
#' ## CAP, Constrained Analysis of Principal coordinates, (Constrained)
#' ### on data described by external class, internal distance matrix calculation
#' dune_bray_cap <- capscale(dune ~ Management, data = dune.env, distance = dist_type)
#' ## dbRDA, distance-based ReDundancy Analysis, (Constrained)
#' ### on data described by external class, internal distance matrix calculation
#' dune_bray_dbrda <- dbrda(dune ~ Management, data = dune.env, distance = dist_type)
#' 
#' # Compare the structures:
#' str(dune_bray_pcoa)    # Unconstrained, {stats}, len 5: $points, $eig, $GOF
#' str(dune_bray_pco)     # Unconstrained: len 12
#' str(dune_bray_cap)     # Constrained: len 14 $CA$eig (unconstrained), $CA$u? (constrained)
#' str(dune_bray_dbrda)   # Constrained: len 13$CA (unconstrained), $CCA (constrained)
#' names(dune_bray_pco)   # Unconstrained: len 12
#' names(dune_bray_cap)   # Constrained: len 14 $CA$eig (unconstrained), $CA$u? (constrained)
#' names(dune_bray_dbrda) # Constrained: len 13$CA (unconstrained), $CCA (constrained)
#' class(dune_bray_pco)   # Unconstrained: len 12
#' class(dune_bray_cap)   # Constrained: len 14 $CA$eig (unconstrained), $CA$u? (constrained)
#' class(dune_bray_dbrda) # Constrained: len 13$CA (unconstrained), $CCA (constrained)
#' dune_bray_pco$Ybar %>% dim
#' 
#' # All eigenvalues (constrained + unconstrained)
#' dune_bray_cap_eig <- c(dune_bray_cap$CCA$eig, dune_bray_cap$CA$eig)
#' dune_bray_cap_scree <- dune_bray_cap_eig %>% scree_df()
#' scree_plot(dune_bray_cap_scree)
#' 
#' 
#' # Screeplot
#' dune_bray_pcoa_scree <- dune_bray_pcoa$eig %>% scree_df()
#' scree_plot(dune_bray_pcoa_scree)
#' 
#' # Estimate number of dimensions to embed
#' (dune_ide <- ide_math(dune, inc_slow = FALSE))
#' (accepted_d <- dune_ide %>% median(na.rm = TRUE) %>% round(0))
#' 
#' # Visualized accepted embedding space
#' dune_beta_embed <- dune_bray_pcoa$points[, 1:accepted_d] %>% as_tibble()
#' 
#' 
#' pairs(dune_beta_embed)
#' GGally::ggpairs(dune_beta_embed)
