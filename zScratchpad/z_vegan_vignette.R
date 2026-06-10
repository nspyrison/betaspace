library(vegan)


## 1 Load data ----
# Load community matrix (species abundance)
data(dune)
# Load environmental metadata
data(dune.env)

head(dune[, 1:5], n = 3)  # View first 5 species for context
head(dune.env, n = 3)    # View metadata (Management, Moisture, Use, Manure)


## 2 Calc Distance ----
# Compute Bray-Curtis dissimilarity matrix
dune_beta <- vegdist(dune, method = "bray")

# This outputs a symmetric distance matrix comparing every site to every other site
# 0 = completely identical species profiles; 1 = share absolutely zero species


## 3 Unconstrained Ordination: PCoA (calssical MDS) -----
# Run PCoA (Classical Multidimensional Scaling)
dune_pcoa <- cmdscale(dune_beta, k = 10, eig = TRUE)

# Extract coordinates for plotting
pcoa_scores <- as.data.frame(dune_pcoa$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")

# Calculate explained variance per axis
ev <- dune_pcoa$eig / sum(dune_pcoa$eig) * 100
cat(sprintf("Axis 1 explains %.1f%%, Axis 2 explains %.1f%% variance.\n", ev[1], ev[2]))


## 4 Constrained Ordination: db-RDA / CAP -----

"We use Distance-Based Redundancy Analysis (dbrda), which is the modern mathematical equivalent to CAP."

# Run db-RDA using the formula interface (Beta Diversity ~ Environment)
dune_cap <- dbrda(dune ~ Management, data = dune.env, distance = "bray")

# Review the constraint breakdown
summary(dune_cap)


## 5 Statistical Hypothesis Testing -----

# Is the separation we see between management styles statistically meaningful, or just random noise? We run a permutation test (PERMANOVA framework) on our CAP model.

# Global significance test of the constraints
set.seed(42) # For reproducibility
anova(dune_cap, permutations = 999)

# Marginal significance test to see which specific management types matter
anova(dune_cap, by = "margin", permutations = 999)


## 6 Base-R Visualization vs. The "Tidy" Gap -----
# Quick Base-R CAP Plot
plot(dune_cap, type = "n")
points(dune_cap, display = "sites", pch = 21, bg = as.numeric(dune.env$Management))
ordihull(dune_cap, groups = dune.env$Management, draw = "polygon", col = 1:4, alpha = 30)



## Biplot (or Triplot). -----

### For Step 3 (PCoA Biplot) ----
plot(dune_pcoa, type = "t") # "t" mixes text for species and points for sites
## error, so;

# 1. Extract the raw site coordinates from your original cmdscale output
site_scores <- dune_pcoa$points
# 2. Set up the empty plot window using the coordinate limits
plot(site_scores, type = "n", xlab = "PCoA 1", ylab = "PCoA 2")
# 3. Add the site points manually
points(site_scores, pch = 21, bg = "steelblue")
# 4. Add the site labels
text(site_scores, labels = rownames(site_scores), pos = 3, cex = 0.8)


### For Step 4 (db-RDA Triplot) -----
plot(dune_cap, scaling = "species")

### Getting to a ggplot triplot --
library(vegan)
library(ggplot2)

# 1. Run the model (from the previous steps)
data(dune)
data(dune.env)
dune_cap <- dbrda(dune ~ Management, data = dune.env, distance = "bray")

# 2. Extract Site Scores and bind with Metadata
site_scores <- as.data.frame(scores(dune_cap, display = "sites"))
site_scores$Management <- dune.env$Management
site_scores$Use <- dune.env$Use

# 3. Extract Environmental Centroids (for Categorical Management)
centroid_scores <- as.data.frame(scores(dune_cap, display = "bp")) 
# Note: For categorical variables, vegan stores centroids here:
centroid_scores <- as.data.frame(scores(dune_cap, display = "centroids"))
centroid_scores$Variable <- rownames(centroid_scores)

# 4. Extract Species Scores (just taking the top 5 so the plot isn't a mess)
species_scores <- as.data.frame(scores(dune_cap, display = "species"))
species_scores$Species <- rownames(species_scores)
species_scores <- head(species_scores, 5)

# 5. Build the Plot Layer by Layer with Colors and Shapes
library("ggpubr") ## For stat_chull
ggplot() +
  # Layer 1: Site Points (Colored by Management, Shaped by Use)
  geom_point(data = site_scores, aes(x = dbRDA1, y = dbRDA2, color = Management, shape = Use), size = 3) +
  
  # Layer 2: Convex Hulls around the Management groups
  stat_chull(data = site_scores, aes(x = dbRDA1, y = dbRDA2, fill = Management, color = Management), 
             alpha = 0.1, geom = "polygon") +
  
  # Layer 3: Species Labels (The Gray Text)
  geom_text(data = species_scores, aes(x = dbRDA1, y = dbRDA2, label = Species), color = "darkgray", fontface = "italic") +
  
  # Layer 4: Environmental Centroids (The Big Targets)
  geom_point(data = centroid_scores, aes(x = dbRDA1, y = dbRDA2), color = "black", size = 4, shape = 13) +
  geom_text(data = centroid_scores, aes(x = dbRDA1, y = dbRDA2, label = Variable), vjust = -1, fontface = "bold") +
  
  # Styling
  theme_minimal() +
  labs(title = "A True, Functional db-RDA Triplot", x = "dbRDA 1", y = "dbRDA 2")


# Possible package workflow:
#' @examples 
#' ## NOT YET IMPLEMENTED
#' ggtriplot(dune_cap, metadata = dune.env) +
#'   geom_site_points(aes(color = Management, shape = Use)) +
#'   geom_ordihull(aes(fill = Management), alpha = 0.2) +
#'   geom_env_vectors(color = "darkred", arrow = arrow()) +
#'   geom_species_labels(color = "gray40", check_overlap = TRUE)
