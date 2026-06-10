# `betaspace`

An R package exploring Ecology Beta Diversity, especially across dissimilarity metrics and viewing more than 2 components.

# The Scope

The package has a clear opportunity to shine by filling a distinct niche:

- Expects known data structures (like a `vegan` ordination objects paired with a `physeq` metadata data frame)
- "Tidy" bridge for macro- and micro-ecologists alike. Instead of forcing users into a massive ecosystem like `phyloseq` or dealing with the clunky object outputs of `ggordiplots`, your package can focus on providing a seamless, modular grammar for community ordination.
- Interoperable `ggplot` objects, easy to extend and with familar `ggplot2` functions
- Lower the barrier to entry for doing interactive visuals like tooltips and linked brsuhing across visuals with `plotly` and `crosstalk`
- Urge Ecologists to shift toward CAP (Constrained Analysis of Principal coordinates, related to LDA, components by descending cluster separation) and away from a PCoA (related to PCA, Components ordered by descending % Var Explained)
- Better facilitate the Intrinsic Data Dimensionality estimation and encourage using more components in visuals instead stoping at PC1 by PC2 with some small fraction of the separation


## Design process

### Step 1: The Best Architecture to Start With

Do not start by writing `geom_` or `stat_` layers from scratch. Creating true custom `ggplot2` extension `geoms` (using the `ggproto` system) is notoriously complex and doesn't translate perfectly to `plotly` if you are generating synthetic data (like convex hulls or ellipses) on the fly.

Instead, utilize a **Tidy Data-to-Data Pipeline** layout.


#### The Vector-to-Dataframe Paradigm

Design a low-level, tidy computational engine first. Write pure functions that accept a vegan object and return a highly clean, augmented data.frame (or tibble) containing both the original metadata and the calculated coordinates (axes, hulls, ellipses).

Why this is the best starting place:

Tidyverse Compatibility: Once you have a clean data frame with explicit group columns, users can natively pipe it into `ggplot()`, `facet_wrap()`, or `facet_grid()`.

`Plotly`/`Crosstalk` Ready: `crosstalk::SharedData$new()` requires a data frame. If your functions spit out structured tables instead of sealed plot objects, implementing interactive linking becomes trivial for the end-user.


### Step 2: Vetting the Interactivity Intersection

When you begin introducing `plotly`, `crosstalk`, and `facet_wrap()`, you are entering a known "danger zone" in R's graphic ecosystem. You need to vet your data structure against three specific technical constraints:

#### 1. The highlight_key Constraint

For `crosstalk` to link an NMDS point to a row in a data table or a point in another facet, it needs a unique identifier key.

- The Test: Ensure your veg_augment() function automatically retains a row_id or sample_id column.

- When vetting, check if `plotly::highlight_key()` plays nicely with your data frame before passing it to `ggplot2`.

```
R
library(crosstalk)
library(plotly)

# Vet this exact pipeline early in your prototyping
shared_eco <- SharedData$new(ordination_data, key = ~sample_id)

p <- ggplot(shared_eco, aes(x = NMDS1, y = NMDS2, color = Treatment)) +
  geom_point() +
  facet_wrap(~Region)

ggplotly(p) %>% highlight(on = "plotly_selected")
```


#### 2. The Synthetic Data Hull Bug (Crucial for vegan concepts)

This is where `ggordiplots` falls short, and where you can win. If you generate convex hulls (`ordihull`) or ellipses, those geometries are synthetic groups calculated from multiple rows of data.

- The Problem: If a user selects a single point via `crosstalk`, a standard client-side JavaScript filter will instantly break the polygon drawn for the hull because the rest of the points making up that polygon are now filtered out.

- The Design Solution: Your package needs a workflow or documentation explaining how to instantiate a `SharedData` object for the points, but keep the hulls/ellipses attached to a static, unlinked layer so the geometry doesn't collapse into a glitchy mess when a user clicks a point.

#### 3. Plotly + `facet_wrap` Scaling Errors

`plotly::ggplotly()` handles standard `facet_wrap()` reasonably well, but it can struggle with `scales = "free"` or complex polygon layers mapped across facets.

- The Vet: Test how `ggplotly()` serializes your custom shapes. If it fails or looks distorted, your package should offer a secondary native `plotly` function layout (using `plotly::plot_ly()` and `plotly::subplot()`) rather than trying to force everything through `ggplotly()`.


### Step 3: Your Development Checklist & Sandbox

To start vetting this immediately, set up a local R playground script using this exact blueprint:

1. The Sandbox Dataset: Use `data(dune)` and `data(dune.env)` from `vegan`. It’s the gold standard community dataset everyone uses to benchmark.

2. The "Tidy" Test: Try to write a function that takes a `metaMDS(dune)` object and yields a single data frame where `ordihull` coordinates are calculated per group, ready for `geom_polygon()`.

3. The Interactivity Test: Wrap that data frame in `SharedData$new()`, build a faceted ggplot, convert it using `ggplotly()`, and verify that clicking a point in Facet A highlights the corresponding point in Facet B.

If your data pipeline survives that three-step pipeline with clean syntax, you have successfully designed a package architecture that completely outclasses `ggordiplots`.


# Related packages

## 1. The Calculators (The Heavy Lifters)

These packages focus on the hard math—calculating distance matrices, running permutations, and extracting eigenvectors—rather than making pretty graphics.

### `vegan`

- Context & Scope: The absolute titan of community ecology in R. If an ecologist is running an NMDS, PCA, or DCA, they are almost certainly using `vegan::metaMDS` or `vegan::rda`. It handles alpha/beta diversity, null models, and environmental data fitting (`envfit`).
- Visual Limits: It relies heavily on R’s base graphics. While functions like `ordiplot()`, `ordihull()`, and `ordisurf()` are incredibly functional, styling them to modern, publication-quality standards requires a massive amount of tedious, low-level coding.

### `adespatial` & `ade4`

- Context & Scope: Focused on multivariate data analysis with a massive emphasis on *spatial* ecology (e.g., dbMEM, principal coordinates of neighbor matrices). It provides rigorous tools for partitioning beta diversity across space and time.
- Visual Limits: Like `vegan`, its primary output relies on specialized base-R structures, making it difficult for users to casually customize plots or add `ggplot2` layers.

## 2. The All-in-One Ecosystems (Data Holders + Wrappers)

These packages are massive pipelines designed around a custom data structure (often an object that binds an abundance table, metadata, and a phylogenetic tree together). They include their own visualization functions.

### `phyloseq` (Bioconductor)

- Context & Scope: Built specifically for microbial ecology (high-throughput sequencing data). It bundles data into an S4 `phyloseq` object and provides convenient wrapper functions like `plot_ordination()`.
- Visual Limits: It uses ggplot2 under the hood, which is great. However, its visualization functions are rigid wrappers. If a user wants to customize a plot beyond the built-in arguments, they often have to extract the underlying data frame manually or fight the wrapper function. It is also strictly geared toward microbiome data, alienating macro-ecologists (e.g., botanists, marine biologists).

### `microeco`

- Context & Scope: A modern, highly flexible pipeline built on the R6 class system. It is essentially an all-in-one suite for preprocessing, alpha/beta diversity testing, network analysis, and machine learning on microbiome data.  Visual Limits: It generates gorgeous ggplot2 graphics natively, but because it relies on its own proprietary R6 classes (like trans_beta), users are forced to learn the microeco object workflow to use any of its visualization capabilities.

## 3. The Visual Bridges (Where This Package Lives)

These packages explicitly sit between the calculators (`vegan`) and the visualization engine (`ggplot2`).

### `ggord`

- Context & Scope: A lightweight package specifically designed to take ordination outputs from `vegan`, `ade4`, or `stats` and plot them using `ggplot2`. It automatically extracts the scores (axes) and maps them to aesthetics.
- The Gap: It’s great for a quick, standard biplot, but it can be rigid when trying to handle complex community metrics, overlaying custom hulls, or mapping specific beta diversity properties neatly.

### `ggordiplots`

- Context & Scope: This package exists for one highly specific reason: to bring `vegan`'s base-graphics features—like `ordihull()`, `ordiellipse()`, and `ordispider()`—into the ggplot2 world.The Gap: It works by returning a list containing both the ggplot object and the underlying data frames used to calculate the hulls/ellipses. While powerful, the syntax can feel clunky and a bit dated for users who want a seamless tidyverse experience.

