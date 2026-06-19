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
  ensure_vegan_ordination(.ordination)
  .eig_orig <- c(.ordination$CCA$eig, .ordination$CA$eig)
  .eig <- fortify_positive_eigenvalues(.eig_orig, .if_neg_eigen) %>%
    sort(decreasing = T)
  .dist_type <- .ordination$inertia
  .ord_type <- .ordination$method
  
  
  ## IDE
  ### IDE of original data, cheeky grab from parent environment
  .form <- .ordination$call$formula
  .form_lhs <- all.vars(.form[[2]])[1]
  .original_data <- get(.form_lhs, envir = parent.frame())
  .ide_data <- ide_math(.original_data, inc_slow = F)
  .accepted_d_data <- median(.ide_data, na.rm = TRUE) %>% floor()
  ### IDE of row embedding
  #### for dune; 12 rather than raw data's 6. we are usign CAPS doe.
  .row_embed <- scores(.ordination, display = "sites", choices = 1:length(.eig))
  .ide_row_embed <- ide_math(.row_embed, inc_slow = F)
  .accepted_d_row_embed <- median(.ide_row_embed, na.rm = TRUE) %>% floor()
  # ### IDE of col embedding
  # #### !!embedding returns NULL
  # "DO WE LOSE THIS IN CAPS?"
  # .col_embed <- scores(.ordination, display = "species", choices = 1:3)
  # .ide_col_embed <- ide_math(.col_embed, inc_slow = F)
  # .accepted_d_col_embed <- median(.ide_col_embed, na.rm = TRUE) %>% floor()
  ### IDE of BOTH embedding
  #### for dune; 11 rather than raw data's 6. we are usign CAPS doe.
  .all_embed <- scores(.ordination, display = "all", choices = 1:length(.eig))
  .ide_all_embed <- ide_math(
    bind_cols(.all_embed$sites, .all_embed$constraints), inc_slow = F)
  .accepted_d_all_embed <- median(.ide_col_embed, na.rm = TRUE) %>% floor()
  
  
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
    d_embedded = .accepted_d,
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
      scree_df = .scree_df,
      meta_df = .meta_df,
      annotate80pct_df = .annotate80_df
    )
  )
}







#' @examples
#' library(vegan)
#' data(dune) ## Numeric X
#' data(dune.env) ## Site metadata, like class ($Management)
#' dist_type <- "bray" # Distance method
#' 
#' ## CAP, Constrained Analysis of Principal coordinates, (Constrained)
#' dune_bray_cap <- capscale(dune ~ Management, data = dune.env, distance = dist_type)
#' 
#' dune_bray_cap_scree <- scree_df(.ordination = dune_bray_cap)
#' scree_plot(dune_bray_cap_scree)
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
