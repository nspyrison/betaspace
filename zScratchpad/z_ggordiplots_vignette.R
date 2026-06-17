# Install and load packages
# install.packages(c("vegan", "ggplot2", "ggordiplots"))
library(vegan)
library(ggplot2)
library(ggordiplots)

# Load example data
data("dune")
data("dune.env")

# Standardize the data using Hellinger transformation
dune.hel <- decostand(dune, method = "hellinger")

# Run Principal Components Analysis (PCA)
ord <- rda(dune.hel)

# Create PCA plot with hulls and spider lines for treatment groups
pca_plot <- gg_ordiplot(ord, 
                        groups = dune.env$Management, 
                        hull = TRUE, 
                        spiders = FALSE, 
                        ellipse = TRUE)

# Customize the plot using standard ggplot2 functions
pca_plot$plot + 
  theme_bw() + 
  labs(color = "Management Group", 
       # x = "PCA Axis 1", 
       # y = "PCA Axis 2", 
       title = "Dune Vegetation Management") +
  theme(plot.title = element_text(hjust = 0.5))
