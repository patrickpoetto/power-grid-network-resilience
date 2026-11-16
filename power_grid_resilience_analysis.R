## ============================================================================
## Structural and Resilience Analysis of the German Transmission Network
##
## Course:  Statistical Analysis for Network Data
## Authors: Patrick Poetto, Goar Shaboian
##
## Dataset: SciGRID (SciGRID project, funded by the German Federal Ministry of
##          Education and Research) - an automated power-grid network model
##          derived from OpenStreetMap "power" tagging. Snapshot: 2016-07-18.
##          Nodes = substations / plants / generators (geographic centroids).
##          Edges = transmission lines, abstracted as simple links (not the
##          real geographic line path).
##
## Objectives:
##   1. Identify critical nodes/links whose failure disrupts the grid.
##   2. Evaluate the impact of targeted attacks on grid stability.
##   3. Propose resilience-enhancement strategies.
## ============================================================================

library(igraph)
library(sf)
library(dplyr)
library(ggplot2)

set.seed(42)

## ----------------------------------------------------------------------------
## 1. DATA LOADING
##    - SciGRID vertices/edges tables (nodes = stations, edges = lines)
##    - GeoJSON layer with German administrative regions (for spatial context
##      and for validating community-detection results against geography)
## ----------------------------------------------------------------------------

vertices <- read.csv("data/scigrid_vertices.csv", stringsAsFactors = FALSE)
edges    <- read.csv("data/scigrid_edges.csv",    stringsAsFactors = FALSE)

# Administrative regions ("Bundesländer" / voltage-level boundary layer),
# used later to overlay community structure and vulnerability hot-spots.
germany_regions <- st_read("data/germany_admin.geojson", layer = "1_sehr_hoch.geo")
# Simple feature collection: 40 features, 11 fields, MULTIPOLYGON, CRS WGS84

## ----------------------------------------------------------------------------
## 2. GRAPH CONSTRUCTION AND CLEANING
## ----------------------------------------------------------------------------

# Remove duplicated vertices and vertices geographically outside Germany
# (small number of boundary/export artifacts from the OSM extraction).
vertices <- vertices %>% distinct(v_id, .keep_all = TRUE)
vertices <- vertices[!duplicated(vertices$v_id), ]
vertices <- vertices[vertices$lon >= 5.5  & vertices$lon <= 15.5 &
                        vertices$lat >= 47.0 & vertices$lat <= 55.5, ]
# -> 3 duplicated vertices and a handful of out-of-bounds vertices removed

power_network <- graph_from_data_frame(edges, vertices = vertices, directed = FALSE)

# Collapse multi-edges / remove self-loops introduced by the OSM abstraction
power_network <- igraph::simplify(power_network)

cat("Nodes:", vcount(power_network), " Edges:", ecount(power_network), "\n")

## ----------------------------------------------------------------------------
## 3. NETWORK COHESION / CONNECTED COMPONENTS
## ----------------------------------------------------------------------------

is_connected(power_network)   # FALSE -> the grid is NOT a single component

components_info <- components(power_network)
table(components_info$csize)

# Largest connected component ("giant component")
power_subgraph <- decompose(power_network)[[1]]

giant_fraction <- vcount(power_subgraph) / vcount(power_network)
giant_fraction   # ~ 0.9893162

# Result: one giant component holding ~99% of all nodes, plus 5 isolated
# nodes (likely disconnected local generation sites or data artifacts).
# All subsequent structural analysis is carried out on the giant component.

## ----------------------------------------------------------------------------
## 4. CENTRALITY MEASURES
##    degree, closeness, betweenness, eigenvector, local clustering coeff.
## ----------------------------------------------------------------------------

deg        <- degree(power_subgraph)
closen     <- closeness(power_subgraph, normalized = TRUE)
betw       <- betweenness(power_subgraph, normalized = TRUE)
eigen_cent <- eigen_centrality(power_subgraph)$vector
clust_coef <- transitivity(power_subgraph, type = "local", isolates = "zero")

centrality_df <- data.frame(
  node        = V(power_subgraph)$name,
  degree      = deg,
  closeness   = closen,
  betweenness = betw,
  eigenvector = eigen_cent,
  clustering  = clust_coef
)

summary(centrality_df[, -1])

## --- Degree distribution -----------------------------------------------
# Highly right-skewed: most nodes have degree 1-2 (radial / tree-like
# local structure typical of transmission grids). A log-log plot of the
# degree distribution looks roughly linear, suggestive of a power law.
deg_dist <- table(deg)
plot(log(as.numeric(names(deg_dist))), log(as.numeric(deg_dist)),
     xlab = "log(degree)", ylab = "log(frequency)",
     main = "Degree distribution (log-log)")

# Degree vs. average neighbor degree -> disassortative mixing: high-degree
# hub nodes tend to connect to low-degree nodes. This is the same pattern
# seen in the Internet AS-graph and other engineered infrastructure
# networks, and it is a direct consequence of a "hub-and-spoke" topology.
knn_result <- knn(power_subgraph)
plot(knn_result$knn, log = "xy",
     xlab = "degree", ylab = "average neighbor degree",
     main = "Disassortative mixing")

# NOTE: despite the log-log linearity above, a formal power-law fit
# (fit_power_law / poisson.test comparison) does NOT convincingly beat a
# Poisson/exponential null for *degree*. Because this is a spatially
# embedded, (near-)planar network, node degree is physically constrained
# (you cannot route an unbounded number of physical lines into one
# substation) -> the grid is "scale-free-like" in appearance but not a
# true scale-free network in the degree sequence.

## --- Betweenness centrality distribution --------------------------------
# Unlike degree, betweenness DOES show a convincing power-law fit: a small
# number of nodes act as critical, global bottleneck intermediaries on
# shortest paths between otherwise distant parts of the grid.
betw_fit <- fit_power_law(betw[betw > 0])
betw_fit

plot(sort(betw[betw > 0], decreasing = TRUE), log = "xy",
     xlab = "rank", ylab = "betweenness (log)",
     main = "Betweenness centrality (log-log) - linear vs Poisson fit")

## --- Closeness centrality ------------------------------------------------
# Roughly bell-shaped / normal distribution -> fairly uniform reachability
# efficiency across the network (no strongly bimodal "core vs periphery"
# split in terms of how fast a node can reach the rest of the grid).
hist(closen, breaks = 30, main = "Closeness centrality distribution")

## --- Eigenvector centrality and clustering coefficient -------------------
# Both are low/concentrated near zero for most nodes, confirming the
# overall sparseness and low cohesion of the network (few triangles,
# little redundant local connectivity).
hist(eigen_cent, breaks = 30, main = "Eigenvector centrality distribution")
mean(clust_coef)

## ----------------------------------------------------------------------------
## 5. COMMUNITY DETECTION (greedy modularity optimization)
## ----------------------------------------------------------------------------

communities <- cluster_fast_greedy(power_subgraph)
modularity(communities)
length(communities)        # -> 15 clusters identified

V(power_subgraph)$community <- membership(communities)

# Spatial validation: overlaying detected communities on the vertex
# coordinates (and the admin-region GeoJSON) shows denser, more
# heterogeneous clusters in the industrial, densely populated western
# part of Germany, and larger but sparser clusters in less populated
# (mostly eastern/northern) regions.
plot(communities, power_subgraph,
     vertex.label = NA, vertex.size = 3,
     main = "Greedy-modularity communities (15 clusters)")

## ----------------------------------------------------------------------------
## 6. VULNERABILITY / RESILIENCE SIMULATION
##    Node-removal attacks: random vs. targeted-by-centrality.
##    Resilience proxy: size of the giant component after each removal step.
## ----------------------------------------------------------------------------

# Generic attack simulator: removes nodes one at a time (in a given order,
# or re-ranked after every removal for "adaptive" targeted strategies) and
# records the giant-component size at each step.
simulate_attack <- function(graph, order_vec) {
  g <- graph
  n <- vcount(g)
  giant_size <- numeric(n)
  giant_size[1] <- vcount(decompose(g)[[1]])

  for (i in seq_along(order_vec)[-length(order_vec)]) {
    g <- delete_vertices(g, order_vec[i])
    comps <- components(g)
    giant_size[i + 1] <- max(comps$csize)
  }
  giant_size / n
}

n_nodes <- vcount(power_subgraph)
node_names <- V(power_subgraph)$name

## --- Strategy 1: random attack (Monte Carlo average) ---------------------
n_sim <- 100
random_curves <- matrix(NA, nrow = n_sim, ncol = n_nodes)
for (s in seq_len(n_sim)) {
  rand_order <- sample(node_names)
  random_curves[s, ] <- simulate_attack(power_subgraph, rand_order)
}
random_curve <- colMeans(random_curves)

## --- Strategy 2-5: targeted attacks by static centrality ranking ---------
order_degree      <- node_names[order(deg, decreasing = TRUE)]
order_betweenness  <- node_names[order(betw, decreasing = TRUE)]
order_closeness    <- node_names[order(closen, decreasing = TRUE)]
order_eigenvector  <- node_names[order(eigen_cent, decreasing = TRUE)]

curve_degree      <- simulate_attack(power_subgraph, order_degree)
curve_betweenness  <- simulate_attack(power_subgraph, order_betweenness)
curve_closeness    <- simulate_attack(power_subgraph, order_closeness)
curve_eigenvector  <- simulate_attack(power_subgraph, order_eigenvector)

## --- Compare resilience curves --------------------------------------------
resilience_df <- data.frame(
  removed     = 0:(n_nodes - 1),
  random      = random_curve,
  degree      = curve_degree,
  betweenness = curve_betweenness,
  closeness   = curve_closeness,
  eigenvector = curve_eigenvector
)

ggplot(resilience_df, aes(x = removed)) +
  geom_line(aes(y = random,      color = "Random")) +
  geom_line(aes(y = degree,      color = "Degree")) +
  geom_line(aes(y = betweenness, color = "Betweenness")) +
  geom_line(aes(y = closeness,   color = "Closeness")) +
  geom_line(aes(y = eigenvector, color = "Eigenvector")) +
  labs(x = "Nodes removed", y = "Relative giant-component size",
       color = "Attack strategy",
       title = "Resilience of the German transmission network under attack") +
  theme_minimal()

# Number of removals needed to halve the giant component, by strategy
sapply(resilience_df[, -1], function(col) which(col <= 0.5)[1])

## ----------------------------------------------------------------------------
## RESULTS SUMMARY
## ----------------------------------------------------------------------------
# - Random attacks: minimal disruption; the network is robust to random
#   node failure (consistent with a highly heterogeneous degree sequence).
# - Closeness-based attack: steep initial drop - halves the giant component
#   after removing only ~11 nodes - then diminishing returns similar to
#   the random baseline.
# - Degree- and eigenvector-based attacks: comparable to each other,
#   moderately effective, slower than closeness/betweenness initially.
# - Betweenness-based attack: fastest and most efficient collapse of the
#   giant component - the single best attack strategy in terms of fewest
#   removals needed, because betweenness captures a node's global role as
#   a path-intermediary rather than purely local connectivity.
# - Spatial pattern: in densely populated / industrial western regions,
#   disruption requires removing several nodes (redundant paths exist);
#   in sparsely populated regions, removing just 1-2 nodes causes major
#   disruption because there is little to no redundancy.
#
## CONCLUSION
## ----------------------------------------------------------------------------
# The German transmission network displays a small set of highly
# connected / high-betweenness hub nodes among a majority of low-degree
# nodes. This structure makes it robust against random failures but
# clearly vulnerable to targeted attacks on high-betweenness nodes.
# Resilience planning should prioritize protecting these critical
# intermediary nodes, with particular attention to sparsely connected
# (rural) regions, where redundancy - and hence resilience - is lowest.
#
# Future work: extend the analysis to weighted graphs incorporating line
# capacity / power-flow attributes, and replicate the study for other
# European transmission networks as SciGRID-equivalent data becomes
# available.
#
## Key references
## ----------------------------------------------------------------------------
# - Cetinay, H., Devriendt, K., & Van Mieghem, P. (2018). "Nodal
#   Vulnerability to Targeted Attacks in Power Grids."
# - Matke, C., Medjroubi, W., & Kleinhans, D. (2016). SciGRID: An Open
#   Source Reference Model for the European Transmission Network.
# - Medjroubi, W., et al. (2017). "Open Data in Power Grid Modelling: New
#   Approaches Towards Transparent Grid Models."
