# BLADE

**BLADE** (Bayesian Linear Admixture Decomposition and Estimation) is an R package for inferring language admixture proportions and quantifying the genetic, phylogenetic, and geographic contributions in linguistic data. It integrates Bayesian MCMC, F3 statistics, optimal transport (Wasserstein distance), and linear mixed models in a unified framework for comprehensive language evolution studies.

------

## Table of contents

- Key features
- Installation
- Data formats
- Core functions
  - \1. Language admixture decomposition
  - \2. F3 admixture test
  - \3. Spatial Wasserstein distances
  - \4. Linear mixed model for distance matrices
  - \5. Coalescent probability matrix
- Full workflow example
- Output visualisation
- Dependencies
- Contributing
- Citation
- License

------

## Key features

- **Language admixture proportion estimation** – Bayesian linear model with Dirichlet prior on source language contributions, implemented in Stan.
- **Significance testing for linguistic mixing** – F3 statistics with bootstrap (feature resampling) or jackknife (spectral clustering of correlated features) to detect admixture between language varieties.
- **Admixture graph inference for languages** – Builds directed graphs of admixture events from F3 results.
- **Spatial distribution comparison of dialects** – 1‑Wasserstein (Earth Mover’s) distance between point patterns of linguistic subgroups, with projection to planar coordinates.
- **Variance decomposition in linguistic distances** – Linear mixed models partition language distance variation into phylogenetic (genealogical) and spatial (areal) random effects.
- **Coalescent‑based admixture prior** – Compute co‑ancestry probabilities from a language phylogeny to inform admixture models.
- **Publication‑ready plots** – Pie charts, bar plots with credible intervals, boxplots of posterior draws, and admixture graphs.

------

## Installation

You can install the development version of BLADE from GitHub:

r

```
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("yourusername/BLADE")   # replace with actual URL
```



**Important**: The package depends on **rstan**, which requires a C++ compiler.
Please follow the [rstan installation guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for your operating system.
On Windows, you may need [Rtools](https://cran.r-project.org/bin/windows/Rtools/). On macOS, ensure Xcode command line tools are installed.

------

## Data formats

Most functions operate on data frames or matrices of linguistic feature frequencies (e.g., presence/absence of typological features, cognate proportions, or lexical percentages):

- **`freq.dat`** – A data frame where rows are **language varieties** (languages, dialects, or lects) and columns are **linguistic features** (e.g., typological binary features, cognate sets, or phoneme inventories). Values are typically frequencies (0–1).
- **`Y.name`** – Character string naming the target language variety (must be a row name of `freq.dat`).
- **`X.name`** – Character vector naming the source language varieties (must be row names of `freq.dat`).

For spatial functions:

- **`df`** – Data frame containing at least three columns: longitude, latitude, and a subgroup label (e.g., dialect group or language family).

For LMM and coalescent functions:

- **Distance matrices** – Square, symmetric matrices with row/column names matching language variety labels. These could be lexical distances, typological distances, or phonological distances.

------

## Core functions

### 1. Language admixture decomposition (`decomp.admixture.func`)

This is the main function for estimating admixture proportions among language varieties. It fits the model:

text

```
Y = X * β + ε,   β ~ beta(α1, α2),   ε ~ N(0, σ²)
```



where `Y` is the target language’s feature frequency vector, `X` is a matrix of source language frequencies, and `β` is the vector of admixture proportions (summing to 1). In historical linguistics, this models a language as a mixture of contributions from several source languages (e.g., due to language shift, borrowing, or creolisation).

#### Usage

r

```
decomp.admixture.func(
  freq.dat,
  Y.name,
  X.name,
  mcmc.n = 50000,
  col = c("#BF3EFF", "#FF523F", ...),   # colour palette for plots
  test = "wilcox.test",                 # not used internally (legacy)
  prob.not.co.mat = NULL,               # optional correction matrix
  burn.in = 0.5
)
```



#### Arguments

| Argument          | Description                                                  |
| :---------------- | :----------------------------------------------------------- |
| `freq.dat`        | Data frame with language varieties as rows, linguistic features as columns. |
| `Y.name`          | Name of the target language variety (must be a row name).    |
| `X.name`          | Character vector of source language variety names.           |
| `mcmc.n`          | Total number of MCMC iterations (including warmup).          |
| `col`             | Vector of colours for plotting (one per source).             |
| `prob.not.co.mat` | Optional correction matrix for non‑coancestry (see `Coal.prob.func`). If provided, the function multiplies `X` by its inverse. |
| `burn.in`         | Proportion of iterations to discard as warmup (default 0.95). |

#### Value

A list with the following components:

| Component         | Description                                                  |
| :---------------- | :----------------------------------------------------------- |
| `admixture.prop`  | Posterior mean proportions (named vector).                   |
| `admixture.CI`    | 95% credible interval for each proportion (2.5% and 97.5%).  |
| `admixture.SD`    | Posterior standard deviations.                               |
| `burn.left.scale` | Matrix of posterior draws after burn‑in (rows = iterations, cols = sources). |
| `pie.plt`         | `ggplot2` pie chart object.                                  |
| `bar.plt`         | `ggplot2` bar plot with error bars (credible interval).      |
| `box.plt`         | `ggplot2` boxplot of posterior draws.                        |
| `pie.dat`         | Data frame used for pie/bar plots (proportions and CI).      |

#### Example

r

```
# Assume freq_mat is a data frame with row names = language varieties
# and columns = typological features (binary or frequency)
result <- decomp.admixture.func(
  freq.dat = freq_mat,
  Y.name = "ModernLanguage",
  X.name = c("SourceLangA", "SourceLangB", "SourceLangC"),
  mcmc.n = 50000,
  burn.in = 0.5
)

# Print estimates
result$admixture.prop
#   SourceLangA   SourceLangB   SourceLangC 
#          0.65          0.25          0.10

# Show plots
print(result$pie.plt)
print(result$box.plt)
```



------

### 2. F3 admixture test (`F3.admixture.func`)

The F3 statistic measures whether a target language variety `C` is admixed from two source varieties `A` and `B`:

text

```
F3(C; A, B) = mean((C - A) * (C - B))
```



A significantly negative value indicates admixture (i.e., the target is a mixture of the two sources). This function computes all possible triples of language varieties, tests significance (optionally), and builds an admixture graph.

#### Usage

r

```
F3.admixture.func(
  data,
  coord = NULL,
  group.name = "Language family",
  sig = FALSE,
  method = "bootstrap",
  resample.n = 1000
)
```



#### Arguments

| Argument     | Description                                                  |
| :----------- | :----------------------------------------------------------- |
| `data`       | Data frame where rows are language varieties and columns are linguistic features. The column `group.name` identifies language groups (e.g., language family or subgroup). |
| `coord`      | Optional data frame with columns `"X"` and `"Y"` (or equivalent) for geographic coordinates. If provided, only triples whose convex hulls overlap are considered (spatial constraint – languages that are geographically separated may not plausibly admix). |
| `group.name` | Name of the column in `data` that contains language group labels. |
| `sig`        | Logical; if `TRUE`, perform significance testing (bootstrap or jackknife). |
| `method`     | `"bootstrap"` (resample linguistic features with replacement) or `"jackknife"` (cluster features via spectral clustering and drop whole clusters). |
| `resample.n` | Number of resamples (if `sig = TRUE`).                       |

#### Value

If `sig = FALSE`:

- `freq` – Mean feature frequencies per language group.
- `F3.dat` – Data frame with columns `Var1` (target), `Var2`, `Var3` (sources), and `f3`.
- `admixture.direc` – Adjacency matrix (0/1) of admixture directions (sources → target).

If `sig = TRUE`:

- Additional components:
  - `freq.summary.dat` – Includes F3 estimate, standard error, z‑score, and p‑value (Wilcoxon test against zero).
  - `freq.resample.dat` – All resampled F3 values.
  - `freq.summary.geo.constraint.dat` – Subset of triples that satisfy spatial overlap (if `coord` provided).

The function also plots the admixture graph using `igraph::plot.graph`.

#### Example

r

```
# Simple test without significance
f3_out <- F3.admixture.func(
  data = my_linguistic_data,
  group.name = "LanguageFamily",
  sig = FALSE
)

# Significant test using bootstrap
f3_sig <- F3.admixture.func(
  data = my_linguistic_data,
  group.name = "LanguageFamily",
  sig = TRUE,
  method = "bootstrap",
  resample.n = 500
)

# View significant admixture triples (potential contact-induced mixing)
subset(f3_sig$freq.summary.dat, p.value < 0.05 & f3 < 0)
```



------

### 3. Spatial Wasserstein distances (`subgroup.wasserstein.func`)

Computes the 1‑Wasserstein (Earth Mover’s) distance between the spatial distributions of linguistic subgroups (e.g., dialects, language families). This metric quantifies how much mass must be moved to transform one point pattern (geographic locations of speakers) into another, taking into account Euclidean distances in projected coordinates. It is useful for comparing the geographic dispersion of different language groups.

#### Usage

r

```
subgroup.wasserstein.func(
  df,
  lon_col = "Longitude",
  lat_col = "Latitude",
  subgroup_col = "Subgroup",
  crs_proj = 3857
)
```



#### Arguments

| Argument       | Description                                                  |
| :------------- | :----------------------------------------------------------- |
| `df`           | Data frame containing coordinates and subgroup labels (e.g., dialect names or language families). |
| `lon_col`      | Name of longitude column (WGS84 decimal degrees).            |
| `lat_col`      | Name of latitude column.                                     |
| `subgroup_col` | Name of column defining subgroups.                           |
| `crs_proj`     | EPSG code for projected coordinate system (default 3857 = Web Mercator, meters). |

#### Value

A distance matrix `D` of class `matrix` with row/column names equal to the unique subgroups. The unit is meters (depending on the projection).

#### Details

- Input coordinates are transformed from WGS84 (EPSG:4326) to the specified projected CRS.
- Within each subgroup, points are given equal weight (1 / number of points).
- The ground distance is Euclidean in the projected space.
- The function uses `transport::wasserstein()` with `p = 1`.

#### Example

r

```
# Sample data: 100 speaker communities with coordinates and dialect group labels
set.seed(123)
geo_df <- data.frame(
  Longitude = runif(100, -10, 10),
  Latitude  = runif(100, 35, 60),
  Subgroup  = sample(c("DialectA", "DialectB", "DialectC"), 100, replace = TRUE)
)

D <- subgroup.wasserstein.func(geo_df,
  lon_col = "Longitude",
  lat_col = "Latitude",
  subgroup_col = "Subgroup"
)

# View geographic distances between dialect distributions (in meters)
print(D)
```



------

### 4. Linear mixed model for distance matrices (`LMM.func`)

This function partitions the variation in a distance matrix (e.g., lexical distances between languages) into fixed effects (e.g., geographic distance), phylogenetic random effects (genealogical relatedness), and spatial random effects (areal contact). It performs Principal Coordinates Analysis (PCoA) on the distance matrix and fits a separate LMM to each principal coordinate. This allows estimation of how much of the linguistic variation is due to inheritance vs. contact.

#### Usage

r

```
LMM.func(
  Y.dist,
  PC.num = 2,
  X.dist = NULL,
  Phy.id,
  Spatial.id,
  phy.cov.matrix,
  spatial.cov.matrix,
  iter = 5000,
  warmup = 0.9,
  chains = 1,
  seed = 0
)
```



#### Arguments

| Argument             | Description                                                  |
| :------------------- | :----------------------------------------------------------- |
| `Y.dist`             | A distance matrix (response), e.g., lexical distances between language varieties. |
| `PC.num`             | Number of principal coordinates to retain (default 2).       |
| `X.dist`             | Optional distance matrix for fixed effects (e.g., geographic distance). If `NULL`, only an intercept is included. |
| `Phy.id`             | Integer vector of phylogenetic group IDs for each language variety (same order as rows of `Y.dist`). This could represent language family subgroups. |
| `Spatial.id`         | Integer vector of spatial group IDs (e.g., geographical regions). |
| `phy.cov.matrix`     | Covariance matrix for phylogenetic random effects (dimension = number of unique `Phy.id`). Typically derived from a language tree. |
| `spatial.cov.matrix` | Covariance matrix for spatial random effects (e.g., based on geographic distances). |
| `iter`               | Total MCMC iterations.                                       |
| `warmup`             | Proportion of iterations to discard.                         |
| `chains`             | Number of MCMC chains.                                       |
| `seed`               | Random seed.                                                 |

#### Value

A list containing:

- `stan.fit.rlt.list` – List of Stan fit objects (one per principal coordinate).
- `coef.rlt.list` – List of extracted coefficients and R² estimates per PC.
- `R2_full` – Weighted average R² across the selected principal coordinates (weights = eigenvalues / sum of eigenvalues).

#### Model specification (Stan)

text

```
y = Xβ + u_phy + u_spatial + ε
u_phy ~ MVN(0, σ × Σ_phy)
u_spatial ~ MVN(0, σ × Σ_spatial)
ε ~ N(0, σ)
```



The function returns the proportion of variance explained by fixed effects (e.g., geography), phylogeny (inheritance), and spatial random effects (areal contact).

#### Example

r

```
# Suppose we have a lexical distance matrix `lex_dist`
# and a phylogenetic covariance matrix `phy_cov` derived from a language tree
# and a spatial covariance matrix `spat_cov` based on geographic distances
# with group IDs phy_id (language family subgroups) and spat_id (geographic regions)

lmm_out <- LMM.func(
  Y.dist = lex_dist,
  PC.num = 3,
  X.dist = geo_dist,           # geographic distance matrix (fixed effect)
  Phy.id = phy_id,
  Spatial.id = spat_id,
  phy.cov.matrix = phy_cov,
  spatial.cov.matrix = spat_cov,
  iter = 5000
)

# Variance partitioning (R²)
print(lmm_out$R2_full)
#  X.effect.R2 phy.effect.R2 spatial.effect.R2 random.effect.R2
#        0.15          0.55              0.20             0.10
```



------

### 5. Coalescent probability matrix (`Coal.prob.func`)

This function takes a phylogenetic tree of languages (in `ape` or `phylo` format with posterior probabilities) and computes a matrix of co‑ancestry probabilities (i.e., the posterior probability that two language varieties share a common ancestor at a given internal node). It can optionally integrate the results of an LMM to scale the probabilities. The resulting matrix can be used as a prior or correction for admixture analysis.

#### Usage

r

```
Coal.prob.func(
  Tree,
  Y.dist.mat = NULL,
  X.dist.mat = NULL,
  phylo.cor.mat = NULL,
  geo.cor.mat = NULL,
  iter = 1000
)
```



#### Arguments

| Argument                       | Description                                                  |
| :----------------------------- | :----------------------------------------------------------- |
| `Tree`                         | A phylogenetic tree with a `node` and `posterior` column (typically from BEAST or MrBayes output, e.g., from language phylogenetics). |
| `Y.dist.mat`                   | Optional distance matrix (response) for LMM. If provided, the function runs `LMM.func` to estimate the proportion of variance explained by phylogeny (`R2`). |
| `X.dist.mat`                   | Optional fixed‑effect distance matrix for LMM.               |
| `phylo.cor.mat`, `geo.cor.mat` | Covariance matrices for LMM random effects.                  |
| `iter`                         | Number of iterations for LMM (if `Y.dist.mat` provided).     |

#### Value

A list:

- `Adm.prob.mat` – Admixture probability matrix: `1 - Anc.prob.mat * coal.prob.coef`.
  Here `Anc.prob.mat[i,j]` is the posterior probability that language varieties `i` and `j` coalesce before any other variety, and `coal.prob.coef` is either 1 (if no LMM) or the R² of the phylogenetic effect from the LMM.
- `LMM.rlt` – Output from `LMM.func` (if `Y.dist.mat` was supplied), otherwise `NULL`.

This matrix can be passed as `prob.not.co.mat` to `decomp.admixture.func` to adjust the source frequency matrix for non‑coancestry (i.e., down-weighting sources that are too closely related).

#### Example

r

```
# Tree with posterior probabilities (e.g., from a Bayesian language phylogeny)
# Assume tree_df has columns "node", "posterior"
tree <- read.nexus("language_tree.nex")
# Convert to data frame with node posteriors if needed

coal_out <- Coal.prob.func(
  Tree = tree,
  Y.dist.mat = lex_dist,
  X.dist.mat = geo_dist,
  phylo.cor.mat = phy_cov,
  geo.cor.mat = spat_cov,
  iter = 2000
)

# Use the resulting matrix to correct admixture decomposition (avoiding pseudo-admixture from related sources)
admix_result <- decomp.admixture.func(
  freq.dat = feature_freq,
  Y.name = "TargetLang",
  X.name = c("SourceA","SourceB","SourceC"),
  prob.not.co.mat = coal_out$Adm.prob.mat
)
```



------

## Full workflow example

Below is a complete (simulated) example that ties together several functions for a language contact scenario.

r

```
library(BLADE)

# ---------- 1. Simulate linguistic feature frequency data ----------
set.seed(42)
langs <- c("TargetLang", "Source1", "Source2", "Source3")
n_features <- 200
freq_mat <- matrix(runif(length(langs) * n_features, 0, 1),
                   nrow = length(langs), ncol = n_features,
                   dimnames = list(langs, paste0("F", 1:n_features)))

# ---------- 2. Admixture decomposition ----------
admix <- decomp.admixture.func(
  freq.dat = freq_mat,
  Y.name = "TargetLang",
  X.name = c("Source1", "Source2", "Source3"),
  mcmc.n = 50000,
  burn.in = 0.5
)

print(admix$admixture.prop)
# Source1 Source2 Source3 
#    0.52    0.30    0.18

# ---------- 3. F3 admixture test ----------
# Need a data frame with language group labels
freq_long <- as.data.frame(freq_mat)
freq_long$LanguageGroup <- rownames(freq_long)

f3_out <- F3.admixture.func(
  data = freq_long,
  group.name = "LanguageGroup",
  sig = TRUE,
  method = "bootstrap",
  resample.n = 100
)

# Show significant triples (likely contact-induced mixing)
f3_out$freq.summary.dat[f3_out$freq.summary.dat$p.value < 0.05, ]

# ---------- 4. Spatial Wasserstein distances for dialect groups ----------
# Create some dummy geographic data for speaker communities
geo_dat <- data.frame(
  Longitude = runif(100, -5, 5),
  Latitude  = runif(100, 40, 50),
  Subgroup  = sample(c("Source1", "Source2", "Source3"), 100, replace = TRUE)
)
D_geo <- subgroup.wasserstein.func(geo_dat,
  lon_col = "Longitude",
  lat_col = "Latitude",
  subgroup_col = "Subgroup"
)
print(D_geo)

# ---------- 5. LMM for variance decomposition in lexical distances ----------
# Simulate a lexical distance matrix (Euclidean on feature frequencies)
lex_dist <- dist(freq_mat) %>% as.matrix()
# Simulate a geographic distance matrix between language areas
geo_dist <- dist(geo_dat[,1:2]) %>% as.matrix()
# Create phylogenetic and spatial group IDs (all same for demo)
phy_id <- rep(1, nrow(lex_dist))
spat_id <- rep(1, nrow(lex_dist))
# Covariance matrices: identity for demo
phy_cov <- diag(1, nrow(lex_dist))
spat_cov <- diag(1, nrow(lex_dist))

lmm_res <- LMM.func(
  Y.dist = lex_dist,
  PC.num = 2,
  X.dist = geo_dist,
  Phy.id = phy_id,
  Spatial.id = spat_id,
  phy.cov.matrix = phy_cov,
  spatial.cov.matrix = spat_cov,
  iter = 1000
)

print(lmm_res$R2_full)
```



------

## Output visualisation

All plotting functions return `ggplot2` objects that can be further customised:

r

```
# Pie chart with custom theme
admix$pie.plt + theme(legend.position = "right")

# Boxplot with adjusted axis text
admix$box.plt + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Admixture graph (from F3)
plot(graph_from_adjacency_matrix(f3_out$admixture.direc))
```



------

## Dependencies

BLADE imports the following packages (installed automatically):

- **rstan** – Bayesian MCMC
- **dplyr, tibble, reshape2** – Data manipulation
- **ggplot2** – Graphics
- **ape, igraph** – Phylogenetics and graph handling
- **sf, transport** – Spatial operations and Wasserstein distance
- **MASS, kernlab** – Matrix algebra and spectral clustering

**System requirements**: A working C++ compiler (for rstan). On Windows, Rtools; on macOS, Xcode command line tools; on Linux, g++.

------

## Contributing

We welcome contributions! Please:

1. Fork the repository.
2. Create a feature branch.
3. Make your changes, add tests if possible.
4. Run `R CMD check` to ensure no errors.
5. Open a pull request.

For bug reports or feature requests, please use the [GitHub Issues](https://../issues) page.

------

## Citation

If you use BLADE in a publication, please cite:

> 

------

## License

BLADE is released under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).
