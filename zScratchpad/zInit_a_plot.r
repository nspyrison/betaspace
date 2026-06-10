## Setup -----
#install.packages("vegan")
library(vegan)
library(spinifex)
library(dplyr)
library(ggplot2)
library(fcuk)

data(dune)
dat <- dune %>%
  as_tibble()

## Fucntions -----
#' @examples
#' library(vegan)
#' data(dune)
#' 
#' # Compute Bray-Curtis dissimilarity matrix
#' dist_matrix <- vegdist(dune, method = "bray")
#' # Perform PCoA (Principal Coordinates Analysis)
#' pcoa_result <- cmdscale(dist_matrix, k = nrow(dune) - 1, eig = TRUE)
#' # Extract eigenvalues for scree plot
#' eigenvalues <- pcoa_result$eig
#' scree_df(.eigen = eigenvalues, .method = "Bray-Curtis")
scree_df <- function(
  .eigen = NULL, .method = NULL, class = NULL,
  .if_neg_eigen = c("rmv", "less_min")
){
  .if_neg_eigen <- match.arg(.if_neg_eigen, several.ok = FALSE)
  .eignen_neg_idx <- .eigen < 0
  if(any(.eignen_neg_idx)){
    ## Deal with it
    if(.if_neg_eigen == "rmv"){
      .eigen <- .eigen[!.eignen_neg_idx]
    }
    if(.if_neg_eigen == "less_min"){
      .eigen <- .eigen - min(.eigen)
    }
  
    ## Report
    .interp <- case_when(
      .if_neg_eigen == "rmv" ~
        ".if_neg_egin is 'rmv'; removed negative eigenvalues.",
      .if_neg_eigen == "less_min" ~
        ".if_neg_egin is 'less_min'; subtracted the most negative from all eigenvalues.",
      TRUE ~ "<!!Issue with mapping!!>")
    .n <- sum(.eignen_neg_idx)
    .N <- length(.eignen_neg_idx)
    .msg <- paste0(
      "!Negative Eigenvalues; There were ", .n, " of ", .N, 
      " (or ", round(100 * .n / .N, 1), "% of columns), ", .interp)
    warning(.msg)
  }
  
  ## On to the scree tbl
  scree_df <- tibble(
    method = .method,
    component = 1:NROW(.eigen),
    class = class,
    var = .eigen,
    pct_var = .eigen / sum(.eigen) * 100,
    cumsum_var = cumsum(pct_var) 
  )
  #' @examples 
  #' tol = .01; max(scree_df$cumsum_var) %>% between(100-tol, 100+tol)
  #' ## Print version, char annoyingly
  #' scree_df %>%
  #'   mutate(
  #'     var = format(var, digits = 4, nsmall = 4),
  #'     pct_var = format(pct_var, digits = 2, nsmall = 2),
  #'     cumsum_var = format(cumsum_var, digits = 2, nsmall = 2)
  #'   )
  
  return(scree_df)
}

#' @examples
#' library(vegan)
#' data(dune)
#' 
#' # Compute Bray-Curtis dissimilarity matrix
#' dist_matrix <- vegdist(dune, method = "bray")
#' # Perform PCoA (Principal Coordinates Analysis)
#' pcoa_result <- cmdscale(dist_matrix, k = nrow(dune) - 1, eig = TRUE)
#' # Extract eigenvalues for scree plot
#' eigenvalues <- pcoa_result$eig
#' (idd <- ide(mtcars))
#' median(idd)
ide <- function(
    data, inc_slow = FALSE
){
  ls_funcs <- list(Rdimtools::est.boxcount, Rdimtools::est.correlation, 
                   Rdimtools::est.made, Rdimtools::est.mle2, Rdimtools::est.twonn)
  nms <- c("est.boxcount", "est.correlation", "est.made", "est.mle2", 
           "est.twonn")
  if (inc_slow == TRUE) {
    ls_funcs <- c(ls_funcs, list(Rdimtools::est.clustering, 
                                 Rdimtools::est.danco, Rdimtools::est.gdistnn, Rdimtools::est.incisingball, 
                                 Rdimtools::est.mindkl, Rdimtools::est.Ustat))
    nms <- c(nms, "est.clustering", "est.danco", "est.gdistnn", 
             "est.incisingball", "est.mindkl", "est.Ustat")
  }
  ret <- sapply(1:length(ls_funcs), function(i) {
    tryCatch(ls_funcs[[i]](data)$estdim, error = function(cond) {
      message("Error in an Rdimtools::est.* function (reported as NA)")
      message(cond)
      return(NA)
    })
  })
  names(ret) <- c(nms)
  return(ret)
}

# The loop ----
.dists <- c("bray", "jaccard", "euclidean")
.scree_list <- list()
for(i in seq_along(.dists)){
  .meth <- .dists[i]
  .dist <- vegdist(dat, method = .meth) %>% tibble()
  .pcoa <- cmdscale(.dist, k = nrow(dat) - 1, eig = TRUE)
  
  # Create scree dataframe
  .scree <- scree_df(.pcoa$eig, .method = .meth)
  
  # Store with a meaningful name
  .scree_list[[.meth]] <- .scree
}
# Bind all rows, preserving the method as an ID column
.all_scree <- bind_rows(.scree_list, .id = "method")

spinifex:::ide()

## Visual
ggplot(.all_scree) +
  geom_col(aes(x = component, y = pct_var, fill = method)) +
  geom_line(aes(x = component, y = cumsum_var)) +
  geom_point(aes(component, cumsum_var)) +
  facet_wrap(~ method, scales = "free_y") +
  labs(
    x = paste0("PCo1, ", round(.all_scree$pct_var[1], 1), "% var"),
    y = paste0("PCo2, ", round(.all_scree$pct_var[2], 1), "% var")
  ) +
  theme_minimal() +
  theme(legend.position = "none")




## Boilerplate code ----
#' @examples
#' library(vegan)
#' data(dune)
#' 
#' # Compute Bray-Curtis dissimilarity matrix
#' dist_matrix <- vegdist(dune, method = "bray")
#' # Perform PCoA (Principal Coordinates Analysis)
#' pcoa_result <- cmdscale(dist_matrix, k = nrow(dune) - 1, eig = TRUE)
#' # Extract eigenvalues for scree plot
#' eigenvalues <- pcoa_result$eig
#' 
#' # Calculate proportion of variance explained
#' variance_explained <- eigenvalues / sum(eigenvalues)
#' 
#' # Cumulative variance explained
#' cumulative_variance <- cumsum(variance_explained)
#' print(cumulative_variance)  # Show first five axes
#' 
#' # Scree plot
#' plot(variance_explained,
#'      type = "b",
#'      main = "PCoA Scree Plot",
#'      xlab = "Principal Coordinate",
#'      ylab = "Proportion of Variance Explained",
#'      pch = 19)
#' abline(h = 0.05, lty = 2, col = "gray")  # Optional: a 5% variance threshold line
#' 
#' 
#' # Calculate alpha diversity (Shannon index)
#' alpha_shannon <- diversity(dune, index = "shannon")
#' 
#' # Simple boxplot of alpha diversity across sites
#' boxplot(alpha_shannon,
#'         main = "Alpha Diversity (Shannon Index) across Sites",
#'         ylab = "Shannon Index")
#' 
#' # For beta diversity, we already have the distance matrix
#' # Let's visualize it with a heatmap
#' heatmap(as.matrix(dist_matrix),
#'         main = "Beta Diversity (Bray-Curtis) Heatmap",
#'         xlab = "Sites", ylab = "Sites",
#'         col = heat.colors(100))
