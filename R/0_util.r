# GUIDING WORKFLOW -----
## SMALL n, LARGE P... :( ----
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


## binned classed, no taxa ----
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

## Bray, PCoA, GGally -----
#' @examples
#' library(vegan)
#' library(phyloseq)
#' 
#' # Load example data
#' data(dune)
#' dune_bray <- vegdist(dune, method = "bray")
#' 
#' 
#' # Perform PCoA (Principal Coordinates Analysis)
#' # dune_bray_pcoa <- cmdscale(dune_bray, k = nrow(dune) - 1, eig = TRUE)
#' # dune_bray_pcoa$points %>% dim() ## note p by truncated-p dimensions
#' # Constrained PCoA (CAP) with Management as the grouping variable
#' dune_bray_cap <- capscale(dune ~ Management, data = dune.env, distance = "bray")
#' 
#' # Compare the structures:
#' str(dune_bray_pcoa)   # Unconstrained: $points, $eig, $GOF
#' str(dune_cap)         # Constrained: $CA (unconstrained), $CCA (constrained)
#' 
#' 
#' # All eigenvalues (constrained + unconstrained)
#' dune_bray_cap_eig <- c(dune_cap$CCA$eig, dune_cap$CA$eig)
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


# Functions -----
ensure_vegan_ordination <- function(.ordination){
  .clas <- any(c("vegan_pco", "capscale", "dbrda") %in% class(.ordination))
  #.len <- length(.ordination) > 10
  .nms <- all(c("Ybar", "method", "call", "CA") %in% names(.ordination))
  if(!all(.clas, .nms))
    stop("Supplied .ordination wasn't made with {vegan}, please use a return from pco(), dbrda(), or capscale().")
  invisible()
}

fortify_positive_eigenvalues <- function(
    .eig, .if_neg_eigen #= c("rmv", "add")
){
  # .if_neg_eigen <- match.arg(.if_neg_eigen, several.ok = FALSE)
  .eig_neg_idx <- .eig < 0
  
  ## If negative eigen values, deal with them
  if(any(.eig_neg_idx)){
    if(.if_neg_eigen == "rmv"){
      .eig <- .eig[!.eig_neg_idx]
    }
    if(.if_neg_eigen == "add"){ ## Add most negative about, so all positive
      .eig <- .eig - min(.eigen)
    }
    
    ## Report
    .interp <- dplyr::case_when(
      .if_neg_eigen == "rmv" ~
        ".if_neg_egin is 'rmv'; removed negative eigenvalues.",
      .if_neg_eigen == "add" ~
        ".if_neg_egin is 'add'; add to all eigenvalues, now the lowest is 0.",
      TRUE ~ "<!!Issue with mapping!!>")
    .n <- sum(.eig_neg_idx)
    .N <- length(.eig_neg_idx)
    .msg <- paste0(
      "Note: contained negative Eigenvalues (dissimilarity not fully Eclidian, fairly common); ", .n, " of ", .N, 
      " eigen values were negative (or ", round(100 * .n / .N, 1), "%), ", .interp)
    warning(.msg)
  }
  
  .eig
}




#' Common distance/dissimilarity metrics for biodiversity analysis
#'
#' A curated list of methods available in R, primarily via `vegan::vegdist()`
#' and `stats::dist()`, useful for beta diversity and ordination.
#'
#' @format A data frame with 3 columns:
#' \describe{
#'   \item{method}{Character, the method name as used in `vegdist` or `dist`.}
#'   \item{type}{Character, ecological or general statistical.}
#'   \item{description}{Character, brief context on typical use case.}
#' }
#' @export
#' @exmaples
#' distances_expl()
distances_expl <- function(){
  tibble::tribble(
    ~method,          ~type,           ~description,
    "bray",           "ecological",    "Bray-Curtis; abundance-based, robust to zeros.",
    "jaccard",        "ecological",    "Jaccard; presence/absence, binary dissimilarity.",
    "euclidean",      "general",       "Euclidean; geometric distance in variable space.",
    "manhattan",      "general",       "Manhattan; sum of absolute differences.",
    "gower",          "general",       "Gower; mixed variable types, range-scaled.",
    "kulczynski",     "ecological",    "Kulczynski; abundance-weighted, similar to Bray.",
    "horn",           "ecological",    "Morisita-Horn; emphasizes dominant species.",
    "mountford",      "ecological",    "Mountford; for presence/absence, rare species sensitive.",
    "raup",           "ecological",    "Raup-Crick; probabilistic, null-model adjusted.",
    "chao",           "ecological",    "Chao; adjusts for unseen shared species (estimation).",
    "aitchison",      "compositional", "Aitchison; for compositional data (log-ratio).",
    "mahalanobis",    "general",       "Mahalanobis; accounts for covariance structure.",
    "unifrac",        "phylogenetic",  "Unweighted UniFrac; presence/absence, requires phylogenetic tree (phyloseq/GUniFrac).",
    "wunifrac",       "phylogenetic",  "Weighted UniFrac; abundance-weighted, requires phylogenetic tree (phyloseq/GUniFrac)."
  ) %>% dplyr::arrange(type, method)
}


#' #' @examples
#' #' library(vegan)
#' #' 
#' #' # Load example data
#' #' data(dune)
#' #' dune_bray <- vegdist(dune, method = "bray")
#' #' is_dist(dune_bray)
#' #' dune_pcoa <- cmdscale(dune_bray, k = 10, eig = TRUE)
#' is_dist <- function(x) {
#'   inherits(x, "dist")
#' }
#' #' @examples
#' #' is.matrix(dune_bray)
#' #' dim(dune_bray) ## ! NULL
#' #' length(dune_bray) ## 190
#' #' NROW(dune_bray) ## ! 190
#' #' (p <- attributes(dune_bray)$Size) ## 20; GOOD
#' #' p * (p - 1) / 2 == length(dune_bray) ## TRUE; FITS EXPECTED LENGTH, GOOD.




#' @examples
#' library(vegan)
#' data(dune) ## Numeric X
#' data(dune.env) ## Site metadata, like class ($Management)
#' dist_type <- "bray" # Distance method
#' 
#' ## CAP, Constrained Analysis of Principal coordinates, (Constrained)
#' dune_bray_cap <- capscale(dune ~ Management, data = dune.env, distance = dist_type)
#' 
#' scree_df(.ordination = dune_bray_cap)
scree_df <- function(
  .ordination, .if_neg_eigen = c("rmv", "add")
){
  .if_neg_eigen <- match.arg(.if_neg_eigen, several.ok = FALSE)
  ensure_vegan_dist(.ordination)
  .eig_orig <- c(.ordination$CCA$eig, .ordination$CA$eig)
  .eig <- fortify_positive_eigenvalues(.eig_orig, .if_neg_eigen) %>%
    sort(decreasing = T)
  .dist_type <- .ordination$inertia
  .ord_type <- .ordination$method
  
  
  ## IDE
  .ordination %>% names
  ### Get the original data in the parent environment
  .form <- .ordination$call$formula
  .form_lhs <- all.vars(.form[[2]])[1]
  .original_data <- get(.form_lhs, envir = parent.frame())
  .ide <- ide_math(.original_data, inc_slow = F)
  .accepted_d <- median(.ide, na.rm = TRUE) %>% floor()
  
  ## Embedding space
  ## TODO NEED TO REVIEW MADISON's EXAMPLES TO SEE WHAT THEY DO.,. they use a completely different pkg...
  ## phyloseq::ordinate(full_meta, method='PCoA',distance='bray')
  #' @examples
  #' if(F){
  #'   ## Species/sample and site scores
  #'   .ordination$CA$u
  #'   scores(.ordination)
  #'   
  #'   # For the constrained axes (the supervised embedding)
  #'   embedding <- scores(.ordination, display = "lc", choices = 1:.accepted_d)
  #'   head(embedding)
  #'   
  #'   # If you want the full set of scores (constrained + unconstrained)
  #'   all_scores <- scores(.ordination, display = c("lc", "wa"), choices = 1:.accepted_d)
  #'   head(all_scores)
  #' }
  
  
  ## On to the scree tbl
  .scree_df <- tibble(
    dist_type = .dist_type,
    ord_type = .ord_type,
    if_neg_eigen = .if_neg_eigen,
    comp_num = 1:NROW(.eig),
    comp_nm = names(.eig),
    comp_type = names(.eig) %>% gsub("[0-9]", "", .),
    var = .eig,
    pct_var = .eig / sum(.eig) * 100,
    cumsum_pct_var = cumsum(pct_var)
  )
  #' @examples
  #' ## Check tolerance
  #' tol = .01; max(.scree_df$cumsum_pct_var) %>% between(100-tol, 100+tol)
  #' ## Print version, char annoyingly
  #' scree_df %>%
  #'   mutate(
  #'     var = format(var, digits = 4, nsmall = 4),
  #'     pct_var = format(pct_var, digits = 2, nsmall = 2),
  #'     cumsum_pct_var = format(cumsum_pct_var, digits = 2, nsmall = 2)
  #'   )
  
  
  #' @examples
  #' ## DEV annotate 80 pct (var) line
  #' ### For each method, find the first component where cumulative >= 80%
  #' .max_pct <- max(.scree_df$pct_eigen, na.rm = TRUE)
  #' .scale_cum <- .max_pct / 100
  #' .annotate80_df <- .scree_df %>%
  #'   group_by(dist_type, ord_type, if_neg_eigen) %>%
  #'   summarize(
  #'     # Component index where cumulative first reaches or exceeds 80%
  #'     n_components = which.max(cumsum_pct_var >= 80),
  #'     # Position for text: rightmost component in the plot
  #'     x_pos = max(comp_num),
  #'     # Position for text: just above the 80% line
  #'     y_pos = 80 * .scale_cum + (0.02 * .max_pct)
  #'   ) %>%
  #'   ungroup()
  
  
  ## Metadata
  .meta_df <- tibble(
    dist_type = .dist_type,
    ord_type = .ord_type,
    if_neg_eigen = .if_neg_eigen,
    d_embeded = .accepted_d,
    n_raw_eigen = length(.eig_orig),
    n_fortified_eigen = length(.eig),
    n_constrained = length(.ordination$CCA$eig),
    n_unconstrained = length(.ordination$CA$eig),
    ord_call = .ordination$call %>% deparse1(),
    scree_call = match.call() %>% deparse1()
  )
  
  
  ## Return
  return(
    list(
      scree_df = .scree_df
      meta_df = .meta_df,
      annotate80pct_df = .annotate80_df
    )
  )
}







#' @examples
#' library(vegan)
#' data(dune)
#' 
#' # Compute Bray-Curtis dissimilarity matrix
#' dune_bray <- vegdist(dune, method = "bray")
#' # Perform PCoA (Principal Coordinates Analysis)
#' dune_bray_pcoa <- cmdscale(dune_bray, k = nrow(dune) - 1, eig = TRUE)
#' # Extract eigenvalues for scree plot
#' dune_eigen <- dune_bray_pcoa$eig
#' dune_scree <- scree_df(.eigen = dune_eigen, .method = "Bray-Curtis")
#' scree_plot(dune_scree)
scree_plot <- function(
    .scree_df = NULL, .components = c(1, 2)
){
  requireNamespace("patchwork")
  
  ## Init ---
  # Determine the maximum of the individual eigenvalue percentages
  max_pct <- max(.scree_df$pct_egien, na.rm = TRUE)
  # Scaling factor: secondary axis will show 0–100
  scale_cum <- max_pct / 100
  
  ## Plot ---
  # g1 <- spinifex::ggtour(bas, data) + spinifex::proto_default(...) + 
  #   ggplot2::labs(x = paste0("Comp", .components[1], ", ", var_df$pct_exp[.components[1]], "% eigen"),
  #                 y = paste0("Comp", .components[2], ", ", var_df$pct_exp[.components[2]], "% eigen"))
  g2 <- 
    ggplot(.scree_df) +
    geom_col(aes(component, pct_egien),
             color = "black", fill = "steelblue3") +
    geom_line(aes(component, cumsum_pct_egien * scale_cum),
              color = "darkorange3", linewidth = 1, linetype = 2) +
    geom_point(aes(component, cumsum_pct_egien * scale_cum),
               color = "darkorange3", size = 3) +
    ## Hline @ 80%
    geom_hline(yintercept = 80 * scale_cum,
               linetype = 3, color = "grey40", linewidth = .5) +
    # Add annotation text
    geom_text(data = annotate80_df,
              aes(x = x_pos, y = y_pos,
                  label = paste0("80% eigen: keep ", n_components, " components")),
              hjust = 1, vjust = 0, size = 4, color = "gray40") +
    scale_y_continuous(
      labels = scales::percent_format(scale = 1),
      sec.axis = sec_axis(~ . / scale_cum,
                          name = "Cummulative Eigenvalue",
                          breaks = seq(0, 100, by = 20),
                          labels = scales::percent_format(scale = 1))
    ) +
    ggplot2::labs(y = "Component Eigenvalue",
                  x = "Component") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.y.left  = element_text(color = "steelblue3", face = "bold"),
          axis.title.y.right = element_text(color = "darkorange3", face = "bold"),
          axis.text.y.left   = element_text(color = "steelblue3"),
          axis.text.y.right  = element_text(color = "darkorange3")
    )
  # Add faceting only if 'method' column exists and has >1 unique value
  if ("method" %in% names(.scree_df)
      #&& length(unique(.scree_df$method)) > 1
  ) {
    g2 <- g2 + facet_grid(cols = vars(method))
  }
  
  #g1 + 
  return(g2)
}




#' Intrisic Dimensionality Estimation
#'
#' A curated list of methods available in R, primarily via `Rdimtools`,
#' mostly differnet mathematically volume pack filling methods. Used to suggest
#' a number of components to keep.
#'
#' @format A data frame with 3 columns:
#' \describe{
#'   \item{method}{Character, the method name as used in `Rdimtools` or `dist`.}
#'   \item{type}{Character, ecological or general statistical.}
#'   \item{description}{Character, brief context on typical use case.}
#' }
#' @exmaples
#' distances_expl()
#' Interactive tour of beta diversity under different distance metrics
#'
#' @param data a numeric data frame.
#' @param inc_slow T/F, whether or not to include some slower `Rdimtool` methods.
#' @return A named vector of suggested intrisic data dimensionality.
#' @export
#' @examples
#' library(vegan)
#' 
#' # Load example data
#' data(dune)
#' dune_bray <- vegdist(dune, method = "bray")
#' 
#' 
#' # Perform PCoA (Principal Coordinates Analysis)
#' # dune_bray_pcoa <- cmdscale(dune_bray, k = nrow(dune) - 1, eig = TRUE)
#' # dune_bray_pcoa$points %>% dim() ## note p by truncated-p dimensions
#' # Constrained PCoA (CAP) with Management as the grouping variable
#' dune_bray_cap <- capscale(dune ~ Management, data = dune.env, distance = "bray")
#' 
#' # Compare the structures:
#' str(dune_bray_pcoa)   # Unconstrained: $points, $eig, $GOF
#' str(dune_cap)         # Constrained: $CA (unconstrained), $CCA (constrained)
#' 
#' 
#' # All eigenvalues (constrained + unconstrained)
#' dune_bray_cap_eig <- c(dune_cap$CCA$eig, dune_cap$CA$eig)
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
ide_math <- function(
    data, inc_slow = FALSE
){
  ## Init which functions
  ls_funcs <- list(
    Rdimtools::est.boxcount, Rdimtools::est.correlation, 
    Rdimtools::est.made, Rdimtools::est.mle2, Rdimtools::est.twonn)
  nms <- c(
    "est.boxcount", "est.correlation", "est.made", "est.mle2", "est.twonn")
  if (inc_slow == TRUE) {
    ls_funcs <- c(ls_funcs, list(
      Rdimtools::est.clustering, Rdimtools::est.danco,
      Rdimtools::est.gdistnn, Rdimtools::est.incisingball,
      Rdimtools::est.mindkl, Rdimtools::est.Ustat))
    nms <- c(nms, "est.clustering", "est.danco", "est.gdistnn",
             "est.incisingball", "est.mindkl", "est.Ustat")
  }
  
  ## Apply functions
  ret <- numeric(length(ls_funcs)) ## Better allocation
  for(i in seq_along(ls_funcs)){
    ret[i] <- tryCatch(
      ls_funcs[[i]](data)$estdim,
      error = function(cond) {
      message(paste0(
        "Error in an Rdimtools::est.* function (reported as NA), function: `",
        nms[i], "`, original error: ", cond$message))
      return(NA)
    })
  }
  ## Fix and return
  names(ret) <- nms
  
  ret
}




# EXPERIMENTAL BELOW ----
is_distance_matrix(dune_bray)
## RETURNS TRUE FOR: inherits dist (TRIANGLE) OR MATRIX that is square, symmetric, and diag near 0.
is_distance_matrix <- function(x) {
  if (inherits(x, "dist")) return(TRUE)
  
  # Check if it's a matrix
  if (!is.matrix(x)) return(FALSE)
  
  # Check dimensions
  if (nrow(x) != ncol(x)) return(FALSE)
  
  # Check symmetry (within tolerance)
  is_symmetric <- all(abs(x - t(x)) < .Machine$double.eps^0.5, na.rm = TRUE)
  
  # Check diagonal is zero (or NA for some ecological distances)
  diag_zero <- all(abs(diag(x)) < .Machine$double.eps^0.5, na.rm = TRUE)
  
  is_symmetric && diag_zero
}



#' Common distance/dissimilarity metrics for biodiversity analysis
#'
#' A curated list of methods available in R, primarily via `vegan::vegdist()`
#' and `stats::dist()`, useful for beta diversity and ordination.
#'
#' @format A data frame with 3 columns:
#' \describe{
#'   \item{method}{Character, the method name as used in `vegdist` or `dist`.}
#'   \item{type}{Character, ecological or general statistical.}
#'   \item{description}{Character, brief context on typical use case.}
#' }
#' @exmaples
#' distances_expl()
#' Interactive tour of beta diversity under different distance metrics
#'
#' @param comm Community matrix (sites x species).
#' @param distances Character vector of distance methods (must be valid for `vegdist`).
#' @param ord_method Ordination method: "PCoA" or "NMDS".
#' @param color Optional factor for coloring sites.
#' @param height,width Plot dimensions for the widget.
#' @return An interactive htmlwidget (plotly + crosstalk).
#' @export
#' @examples
#' library(vegan)
#' library(plotly)
#' library(crosstalk)
#' 
#' # Load example data
#' data(dune)
#' # Run the tour
#' tour_plot <- DEV_tour_beta(
#'   comm = dune,
#'   distances = c("bray", "jaccard", "euclidean"),
#'   ord_method = "PCoA",
#'   color = dune.env$Management
#' )
#' # Display the interactive plot
#' tour_plot
DEV_tour_beta <- function(comm, distances = c("bray", "jaccard", "euclidean"),
                      ord_method = "PCoA", color = NULL, height = 500, width = 800) {
  
  # Pre-compute distance matrices and ordinations for each method
  ord_list <- list()
  for (d in distances) {
    dist_mat <- vegan::vegdist(comm, method = d)
    if (ord_method == "PCoA") {
      ord <- cmdscale(dist_mat, k = 2, eig = FALSE)
    } else if (ord_method == "NMDS") {
      ord <- metaMDS(dist_mat, k = 2, trace = 0)$points
    }
    ord_list[[d]] <- as.data.frame(ord)
    colnames(ord_list[[d]]) <- c("Axis1", "Axis2")
    if (!is.null(color)) ord_list[[d]]$color <- color
  }
  
  # Create a shared data object for crosstalk linking
  shared_data <- crosstalk::SharedData$new(data.frame(site_id = seq_len(nrow(comm))))
  
  # Build a plotly plot for each distance, linked by site_id
  plots <- lapply(distances, function(d) {
    df <- ord_list[[d]]
    df$site_id <- seq_len(nrow(comm))
    
    p <- plotly::plot_ly(df, x = ~Axis1, y = ~Axis2, 
                         color = if (!is.null(color)) ~color else NULL,
                         key = ~site_id, 
                         type = 'scatter', mode = 'markers',
                         text = ~paste("Site:", site_id, "<br>Dist:", d),
                         hoverinfo = 'text',
                         showlegend = FALSE)
    p <- p %>% plotly::layout(title = paste(ord_method, "with", d),
                              xaxis = list(title = "Axis 1"),
                              yaxis = list(title = "Axis 2"))
    p
  })
  
  # Arrange plots in a grid
  subplot <- plotly::subplot(plots, nrows = 1, shareX = TRUE, shareY = FALSE, titleX = TRUE)
  
  # Add a dropdown to highlight a specific site across all plots
  subplot <- subplot %>% plotly::highlight(on = "plotly_click", off = "plotly_doubleclick")
  
  return(subplot)
}


library(vegan)
library(plotly)
library(crosstalk)

# Load example data
data(dune)
data(dune.env)

# Run the tour
tour_plot <- tour_beta(
  comm = dune,
  distances = c("bray", "jaccard", "euclidean"),
  ord_method = "PCoA",
  color = dune.env$Management
)

# Display the interactive plot
tour_plot


#' Suggest elbow point based on scree differences
#'
#' @param eig Numeric vector of eigenvalues (positive only).
#' @param method Character, distance method name (for labeling).
#' @return A tibble with the suggested elbow component and the difference curve.
#' @examples
#' library(vegan)
#' data(dune)
#' 
#' # Compute Bray-Curtis dissimilarity matrix
#' dune_bray <- vegdist(dune, method = "bray")
#' # Perform PCoA (Principal Coordinates Analysis)
#' dune_bray_pcoa <- cmdscale(dune_bray, k = nrow(dune) - 1, eig = TRUE)
#' # Extract eigenvalues for scree plot
#' dune_eigen <- dune_bray_pcoa$eig
#' ide_scree(dune_eigen)
ide_scree_elbow <- function(eigen, method = "unknown") {
  # Ensure eigenvalues are sorted descending
  eig <- sort(eigen, decreasing = TRUE)
  # First differences (drop in eigenvalue)
  diff1 <- -diff(eigen)
  # Second differences (change in drop)
  diff2 <- diff(diff1)
  # Elbow: component where the second difference is maximized (sharpest bend)
  elbow_idx <- which.max(diff2) + 1  # +1 because diff2 is one shorter
  
  tibble::tibble(
    method = method,
    n_components = length(eig),
    elbow_component = elbow_idx,
    elbow_eigenvalue = eig[elbow_idx],
    diff1 = list(diff1),
    diff2 = list(diff2)
  )
}
