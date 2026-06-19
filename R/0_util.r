# Functions -----

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
#' data(dune) ## Numeric X
#' data(dune.env) ## Site metadata, like class ($Management)
#' dist_type <- "bray" # Distance method
#' 
#' ## CAP, Constrained Analysis of Principal coordinates, (Constrained)
#' dune_bray_cap <- capscale(dune ~ Management, data = dune.env, distance = dist_type)
#' 
#' "Typically ide_math is an internal function, but we can do it manually too:"
#' # Estimate number of dimensions to embed
#' # Constrained site scores (same as $CCA$u)
#' site_scores <- scores(dune_bray_cap, display = "lc", choices = 1:999)
#' # Constrained species scores (same as $CCA$v)
#' species_scores <- scores(dune_bray_cap, display = "species", choices = 1:999)
#' # Unconstrained site scores (residual ordination)
#' uncon_scores <- scores(dune_bray_cap, display = "wa", choices = 1:999)
#' 
#' (dune_ide_site <- ide_math(species_scores, inc_slow = FALSE))
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
          "Rdimtools::", nms[i], " failed (reported as NA),",
          " original error:\n", cond$message))
        return(NA)
      })
  }
  ## Fix and return
  names(ret) <- nms
  
  ret
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




fortify_positive_eigenvalues <- function(
    .eig, .if_neg_eigen #= c("rmv", "add"), will from from arg in scree_df, no need for default here
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
















# EXPERIMENTAL BELOW ----

#' @examples
#' library(vegan)
#'
#' # Load example data
#' data(dune)
#' dune_bray <- vegdist(dune, method = "bray")
#' is_dist(dune_bray)
#' dune_pcoa <- cmdscale(dune_bray, k = 10, eig = TRUE)
#' is_vegan_dist(dune_pcoa)
#' 
#' is_vegan_dist(gp_bray_cap)
is_vegan_dist <- function(x) {
  inherits(x, "dist")
}


## TODO finish this thought if needed
is_diag_mat <- function(mat){
  if(is_vegan_dist(mat)){
    ret <- all(names(attributes(mat)) ==
      c("maxdist", "Size", "Labels", "Diag", "Upper", "method", "call", "class"))
  } else {
    if(is_philoseq_dist(mat)){
      ret <- "## TODO"
    }
  } else {
    error("Passed matrix not identifiable to {vegan} or {phyloseq}.")
  }
  
  ret
}

diag_mat_len <- function(daimat){
  if()
  is.matrix(dune_bray)
  dim(dune_bray) ## ! NULL
  length(dune_bray) ## 190
  NROW(dune_bray) ## ! 190
  (p <- attributes(dune_bray)$Size) ## 20; GOOD
  p * (p - 1) / 2 == length(dune_bray) ## TRUE; FITS EXPECTED LENGTH, GOOD.
}
#' @examples
#' is.matrix(dune_bray)
#' dim(dune_bray) ## ! NULL
#' length(dune_bray) ## 190
#' NROW(dune_bray) ## ! 190
#' (p <- attributes(dune_bray)$Size) ## 20; GOOD
#' p * (p - 1) / 2 == length(dune_bray) ## TRUE; FITS EXPECTED LENGTH, GOOD.



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


