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


#' @examples
#' peng <- datasets::penguins %>%
#'   as_tibble() %>%
#'   filter(!is.na(bill_len)) %>% 
#'   select(species, bill_len, bill_dep, flipper_len, body_mass)
#' peng %>% complete.cases %>% all
#' 
#' splom_base(peng)
#' splom_base(peng, color = "species", shape = "species",
#'       title = "Penguins Splom")
#' 
#' ## TODO Make examples for CAP case. 
splom_base <- function(
  df, color = NULL, shape = NULL, size = 2, title = ""
){
  # Setup data and variables
  df <- as_tibble(df)
  # Get numeric columns (exclude color_col if it's in the dataframe)
  if (!is.null(color)) {
    if(is.character(size)){.size <- size}else{.size <- NULL}
    vars <- setdiff(names(df), c(color, shape, .size))
  } else {
    vars <- names(df)
  }
  
  # Filter to numeric columns only
  numeric_vars <- vars[sapply(df[vars], is.numeric)]
  d <- length(numeric_vars)
  if(d < 2) stop("Need at least two numeric variables to plot pairs.")
  
  # Save old par and restore on exit (only save modifiable parameters)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  # Set up the grid layout with room for labels
  par(mfrow = c(d, d), mar = c(1.5, 1.5, 1.5, 1.5), oma = c(2, 2, 2, 2))
  
  # palette setup
  pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
           "#66A61E", "#E6AB02", "#A6761D", "#666666",
           ## "Set2" pal
           "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
           "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  
  if (!is.null(color) && color %in% names(df)) {
    color_factor <- as.factor(df[[color]])
    levels_col <- levels(color_factor)
    # Use Dark2 palette, recycling if needed
    pal <- rep(pal, length.out = length(levels_col))
    point_colors <- pal[as.numeric(color_factor)]
  } else {
    point_colors <- "black"
    color_factor <- NULL
  }
  
  # Handle shape and size parameters
  if (is.character(shape) && shape %in% names(df)) {
    shape_factor <- as.factor(df[[shape]])
    shape_levels <- levels(shape_factor)
    # Map shape levels to pch values: 16, 17, 18, 15, 19, etc.
    shape_pch <- c(16, 17, 18, 15, 19, 21, 22, 23, 24, 25)
    shape_pch <- shape_pch[seq_along(shape_levels)]
    point_shapes <- shape_pch[as.numeric(shape_factor)]
  } else {
    point_shapes <- 19
    shape_factor <- NULL
  }
  
  if (is.character(size) && size %in% names(df)) {
    point_sizes <- df[[size]]
    # Scale to reasonable point size range (0.5 to 3)
    point_sizes <- 0.5 + 2.5 * (point_sizes - min(point_sizes)) / (max(point_sizes) - min(point_sizes))
  } else if (is.numeric(size)) {
    point_sizes <- rep(size, nrow(df))
  } else {
    point_sizes <- rep(2, nrow(df))
  }
  
  # Create the pairwise grid
  for (i in 1:d) {
    for (j in 1:d) {
      
      if (i == j) {
        # Diagonal: density plots
        x_var <- numeric_vars[i]
        # Calculate density for overall range
        dens_overall <- density(df[[x_var]])
        
        # Plot with fixed y-axis at 1
        plot(dens_overall,
             main = "",
             xlab = "",
             ylab = "",
             axes = FALSE,
             col = NA,
             ylim = c(0, 1))
        # Add grid in background like theme_bw()
        grid(nx = NULL, ny = NULL, col = "gray90", lty = 1, lwd = 0.5)
        box(col = "gray30", lwd = 0.5)
        # Add label at top
        mtext(x_var, side = 3, line = 0.2, cex = 0.85, font = 1)
        
        # Add density curves for each group if color_col is provided
        if (!is.null(color_factor)) {
          for (k in seq_along(levels_col)) {
            lvl <- levels_col[k]
            idx <- color_factor == lvl
            sub_data <- df[[x_var]][idx]
            
            if (length(sub_data) >= 2) {
              dens <- density(sub_data)
              lines(dens, col = pal[k], lwd = 2)
              # Soft fill with alpha proportional to number of classes
              alpha <- 0.2 + (0.3 / length(levels_col))
              polygon(dens, col = adjustcolor(pal[k], alpha.f = alpha), border = NA)
            }
          }
        } else {
          dens <- density(df[[x_var]])
          lines(dens, col = "black", lwd = 2)
        }
        
      } else if (i > j) {
        # Lower triangle: scatter + ellipse
        x_var <- numeric_vars[j]
        y_var <- numeric_vars[i]
        
        plot(df[[x_var]], df[[y_var]],
             col = point_colors,
             pch = point_shapes,
             cex = point_sizes,
             xlab = "", ylab = "",
             axes = FALSE,
             xlim = range(df[[x_var]], na.rm = TRUE),
             ylim = range(df[[y_var]], na.rm = TRUE))
        # Add grid in background like theme_bw()
        grid(nx = NULL, ny = NULL, col = "gray90", lty = 1, lwd = 0.5)
        box(col = "gray30", lwd = 0.5)
        # Add y-axis label on left
        mtext(y_var, side = 2, line = 0.2, cex = 0.85, font = 1)
        # Add x-axis label on bottom
        mtext(x_var, side = 1, line = 0.2, cex = 0.85, font = 1)
        
        # Add ellipses for each group if color_col is provided
        if (!is.null(color_factor)) {
          for (lvl in levels_col) {
            idx <- color_factor == lvl
            sub_df <- df[idx, ]
            
            if (nrow(sub_df) >= 3) {
              points_sub <- as.matrix(sub_df[, c(x_var, y_var)])
              
              # Calculate ellipse using covariance
              center <- colMeans(points_sub)
              cov_mat <- cov(points_sub)
              ev <- eigen(cov_mat)
              
              # Create ellipse boundary
              theta <- seq(0, 2*pi, length.out = 100)
              a <- sqrt(ev$values[1])
              b <- sqrt(ev$values[2])
              ellipse_pts <- cbind(a * cos(theta), b * sin(theta))
              
              # Rotate by eigenvector
              ellipse_pts <- ellipse_pts %*% t(ev$vectors)
              
              # Translate to center
              ellipse_pts <- sweep(ellipse_pts, 2, center, "+")
              
              # Draw ellipse with Dark2 color
              lines(ellipse_pts, col = pal[which(levels_col == lvl)], lwd = 1.5)
            }
          }
        }
        
      } else {
        # Empty plots for upper triangle
        plot.new()
      }
    }
  }
  
  # Reset layout
  par(mfrow = c(1, 1))
  
  # Add title if provided
  if (title != "") {
    mtext(title, outer = TRUE, side = 3, line = 0.5, cex = 1.2, font = 2)
  }
  
  invisible(NULL)
}

#' Scatterplot Matrix using ggplot2 (GGAlly-style but faster)
#'
#' Creates an interoperable ggplot2-based scatterplot matrix that can be
#' modified with standard ggplot2 operations. Optimized for speed using
#' vectorized operations instead of looping through individual plots.
#'
#' @param df A data frame or tibble
#' @param color Optional column name for color mapping (categorical or numeric)
#' @param shape Optional column name for shape mapping (categorical)
#' @param size Optional column name or numeric value for point size
#' @param title Optional title for the plot grid
#' @param alpha Transparency of points (0-1)
#' @param lower_fn Function to apply to lower triangle ("scatter", "density", "blank")
#' @param diag_fn Function to apply to diagonal ("density", "histogram", "blank")
#' @param upper_fn Function to apply to upper triangle ("scatter", "density", "blank")
#'
#' @return A patchwork object combining all plots
#'
#' @examples
#' peng <- datasets::penguins %>%
#'   as_tibble() %>%
#'   filter(!is.na(bill_len)) %>% 
#'   select(species, bill_len, bill_dep, flipper_len, body_mass)
#' 
#' splom_gg(peng)
#' splom_gg(peng, color = "species", shape = "species",
#'   title = "betaspace::splom_gg penguins")
#' 
#' # Alternative to:
#' GGally::ggpairs(peng[, -1], aes(color = peng$species, shape = peng$species)) +
#'   ggtitle("GGally::ggpairs penguins")
#' graphics::pairs(peng[2:5], main = "graphics::pairs penguins",
#'   pch = 21, bg = c("red", "green3", "blue")[unclass(peng$species)])
#' splom_base(peng, color = "species", shape = "species",
#'   title = "betaspace::splom_base penguins")
#' @export
splom_gg <- function(
  df, color = NULL, shape = NULL, size = 1.5, title = NULL,
  alpha = 1
) {
  library(ggplot2)
  library(patchwork)
  
  # Setup data and variables
  df <- as_tibble(df)
  
  # Get numeric columns (exclude color_col, shape_col, size_col if they're non-numeric)
  if (!is.null(color)) {
    if (is.character(size)) {
      .size <- size
    } else {
      .size <- NULL
    }
    vars <- setdiff(names(df), c(color, shape, .size))
  } else {
    vars <- names(df)
  }
  
  # Filter to numeric columns only
  numeric_vars <- vars[sapply(df[vars], is.numeric)]
  d <- length(numeric_vars)
  if (d < 2) stop("Need at least two numeric variables to plot pairs.")
  
  # Prepare data
  plot_data <- df |>
    select(all_of(c(numeric_vars, color, shape)),
           everything())
  
  # Determine if color and shape map to the same variable (for unified legend)
  unified_legend <- !is.null(color) && !is.null(shape) && color == shape
  
  # Create list to hold all plots
  plot_list <- vector("list", d * d)
  
  # Create the pairwise grid
  idx <- 1
  for (i in 1:d) {
    for (j in 1:d) {
      
      if (i == j) {
        # Diagonal: density or histogram
        x_var <- numeric_vars[i]
        
        ## Diag density
        {
          p <- ggplot(plot_data, aes(x = .data[[x_var]])) +
            geom_density(fill = "gray70", alpha = 1, linewidth = 0.8) +
            theme_bw() +
            labs(x = "", y = "") +
            theme(
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              #axis.title = element_blank(),
              plot.margin = margin(2, 2, 2, 2, "pt"),
              panel.grid.minor = element_blank()
            )
          
          # Add colored densities if color is specified
          if (!is.null(color) && color %in% names(plot_data)) {
            p <- p +
              geom_density(aes(fill = .data[[color]]), alpha = 1, linewidth = 0.6) +
              scale_fill_brewer(palette = "Dark2", guide = "none")
          }
        }
        
      } else if (i > j) {
        # Lower triangle: scatter plot
        x_var <- numeric_vars[j]
        y_var <- numeric_vars[i]
        
        p <- ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]]))
        
        # Build aesthetic mappings
        aes_params <- list(alpha = alpha)
        
        if (!is.null(color) && color %in% names(plot_data)) {
          aes_params[["color"]] <- as.symbol(color)
        }
        if (!is.null(shape) && shape %in% names(plot_data)) {
          aes_params[["shape"]] <- as.symbol(shape)
        }
        if (is.character(size) && size %in% names(plot_data)) {
          aes_params[["size"]] <- as.symbol(size)
        } else if (is.numeric(size)) {
          aes_params[["size"]] <- size
        } else {
          aes_params[["size"]] <- 1.5
        }
        
        # Build aes call dynamically
        aes_call <- do.call(aes, aes_params[sapply(
          aes_params, function(x) !is.null(x) && x != "")])
        
        # Determine if this is the first (bottom-left) plot for legend placement
        is_legend_plot <- (i == d && j == 1)
        
        p <- p + geom_point(aes_call, na.rm = TRUE)
        
        # Add color scale if needed
        if (!is.null(color) && color %in% names(plot_data)) {
          if (is.numeric(plot_data[[color]])) {
            p <- p + scale_color_viridis_c(guide = if(
              is_legend_plot && unified_legend) "legend" else "none")
          } else {
            p <- p + scale_color_brewer(palette = "Dark2", guide = if(
              is_legend_plot && unified_legend) "legend" else "none")
          }
        }
        
        # Add shape scale if needed
        if (!is.null(shape) && shape %in% names(plot_data)) {
          p <- p + scale_shape_manual(
            values = c(16, 17, 18, 15, 19, 21, 22, 23, 24, 25),
            guide = if(is_legend_plot && unified_legend) "legend" else "none"
          )
        }
        
        # Add size scale if needed
        if (is.character(size) && size %in% names(plot_data)) {
          p <- p + scale_size_continuous(range = c(1, 5), guide = if(
            is_legend_plot) "legend" else "none")
        }
        
        p <- p +
          theme_bw() +
          labs(x = "", y = "") +
          theme(
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            #axis.title = element_blank(),
            plot.margin = margin(2, 2, 2, 2, "pt"),
            panel.grid.minor = element_blank(),
            legend.position = if(is_legend_plot) "right" else "none",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.3, "cm")
          )
        
      } else {
        # Upper triangle: blank
        p <- ggplot() + theme_void()
      }
      
      plot_list[[idx]] <- p
      idx <- idx + 1
    }
  }
  
  # Combine plots using patchwork
  final_plot <- Reduce(function(x, y) x + y, plot_list) +
    plot_layout(ncol = d, nrow = d, guides = "collect")
  
  # Add title if provided
  if (!is.null(title)) {
    final_plot <- final_plot + plot_annotation(title = title)
  }
  
  return(final_plot)
}
