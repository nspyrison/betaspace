
#' @examples
#' library(phyloseq)
#' library(vegan)
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
#' is_vegan_ordination(gp_bray_cap0)
#' gp_bray_cap <- ordinate(gp_class, method = "CCA", distance = "bray", formula = ~ Class)
#' is_vegan_ordination(gp_bray_cap)
is_vegan_ordination <- function(.ordination){
  .clas <- any(c("vegan_pco", "capscale", "dbrda") %in% class(.ordination))
  .nms <- all(c("inertia", "method", "call", "CA") %in% names(.ordination))
  #.len <- length(.ordination) > 10
  
  return(.clas & .nms)
}



#' @examples
#' library(phyloseq)
#' library(vegan)
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
#' is_phyloseq_ordination(gp_bray_cap0)
#' gp_bray_cap <- ordinate(gp_class, method = "CCA", distance = "bray", formula = ~ Class)
#' is_phyloseq_ordination(gp_bray_cap)
is_phyloseq_ordination <- function(.ordination){
  .clas <- any(c("cca") %in% class(.ordination))
  .nms <- all(c("inertia", "method", "call", "CA") %in% names(.ordination))
  #.len <- length(.ordination) > 10
  
  return(.clas & .nms)
}
