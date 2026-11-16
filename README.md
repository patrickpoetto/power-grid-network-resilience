# Structural and Resilience Analysis of the German Transmission Network

**Group project for the course "Statistical Analysis for Network Data"**
Authors: **Patrick Poetto** and **Goar Shaboian**

> This repository contains an applied network-science study of the German
> high-voltage electricity transmission grid, analyzing its structural
> properties and simulating its resilience to random failures and targeted
> attacks.

## Problem

Modern electricity grids are critical infrastructure: the failure of a small
number of key substations or transmission lines can cascade into large-scale
blackouts. This project asks three questions about the German transmission
network:

1. Which nodes (substations, plants, generators) and links are structurally
   *critical* — i.e., their removal disproportionately disrupts the grid?
2. How does the grid's connectivity degrade under random failures versus
   *targeted* attacks on the most "important" nodes, and which notion of
   "importance" makes for the most effective attack?
3. What does the answer imply for where resilience-enhancement investment
   should be prioritized?

## Dataset

We use the **[SciGRID](https://www.scigrid.de/)** power-grid dataset, a
project funded by the German Federal Ministry of Education and Research.
SciGRID automatically derives a network model of the European transmission
grid from OpenStreetMap `power=*` tagging — tagging that, at the time of the
snapshot used here (2016-07-18), was most complete and reliable for Germany.

- **Nodes**: substations, power plants, and generators, represented by their
  geographic centroid.
- **Edges**: transmission lines, abstracted as simple links between nodes
  (not the real-world geographic path of the line).
- The model has been cross-validated by the SciGRID authors against ENTSO-E
  grid maps and the GridKit project.

A supplementary GeoJSON layer of German administrative regions is used to
give geographic context to the network and to visually validate the
community-detection results against real population/industrial density
patterns.

**Preprocessing**: 3 duplicated vertices and a small number of vertices
geographically outside Germany (artifacts of the OSM-based extraction) were
removed before graph construction.

## Methodology

The analysis is implemented in R using the **`igraph`** package, and
proceeds in the following stages (see
[`power_grid_resilience_analysis.R`](power_grid_resilience_analysis.R)):

1. **Graph construction and cleaning** — build an undirected graph from the
   SciGRID vertex/edge tables and simplify it (remove duplicate edges and
   self-loops).
2. **Connectivity / cohesion check** — the raw graph is *not* fully
   connected. It decomposes into one giant component holding ~98.9% of all
   nodes, plus 5 isolated nodes. All subsequent analysis is performed on the
   giant component.
3. **Centrality analysis** — degree, closeness, betweenness, eigenvector
   centrality, and local clustering coefficient are computed for every node
   to characterize each node's structural role.
4. **Community detection** — greedy modularity optimization
   (`cluster_fast_greedy`) partitions the network into 15 clusters, which are
   validated spatially against the administrative-region layer.
5. **Resilience / vulnerability simulation** — nodes are removed one at a
   time under five strategies: random removal (Monte Carlo averaged over
   100 permutations) and targeted removal ranked by degree, closeness,
   betweenness, or eigenvector centrality. After each removal, the size of
   the resulting giant component is recorded as the resilience metric.

## Key findings

- **The degree distribution is highly right-skewed** (most nodes have
  degree 1–2) and looks power-law-like on a log-log plot. However, a formal
  power-law fit does **not** convincingly beat a Poisson/exponential null.
  Because the grid is a spatially embedded, near-planar network, node degree
  is physically constrained — you cannot route unlimited physical lines into
  one substation. The grid is therefore *scale-free in appearance* but not
  truly scale-free in its degree sequence, unlike purely virtual networks
  (e.g. the web graph).
- **Betweenness centrality, by contrast, does show a convincing power-law
  fit.** A small number of nodes act as global, critical intermediaries on
  shortest paths across the grid — these are the structurally irreplaceable
  "bridge" substations.
- **Degree–neighbor-degree correlation is disassortative**: high-degree hub
  nodes tend to connect to low-degree nodes, a pattern also seen in other
  engineered infrastructure networks (e.g., the Internet AS-graph).
- **Community structure correlates with population/industrial density**:
  greedy modularity optimization finds 15 clusters, with denser, more
  heterogeneous communities in industrial western Germany and larger,
  sparser communities in less populated regions.
- **Attack simulations show betweenness-based targeting is the most
  efficient way to disrupt the network** — it collapses the giant component
  with far fewer node removals than degree- or eigenvector-based targeting,
  because it captures a node's *global* path-intermediation role rather than
  purely local connectivity. Closeness-based attacks also cause a very
  steep *initial* drop (halving the giant component after only ~11 node
  removals) before flattening out to behave like the random baseline.
  Random failures, in contrast, cause minimal disruption — the network is
  robust to random node loss.
- **Rural, low-density areas are more vulnerable to single-node failures
  than dense, urban/industrial ones.** In densely populated western
  Germany, redundant paths mean disrupting the grid requires removing
  several nodes; in sparsely populated regions, removing just one or two
  nodes can cause major disruption due to the near-total lack of redundant
  routing.

### Conclusion

The German transmission network exhibits the classic signature of
engineered infrastructure: a few highly connected, high-betweenness hub
substations embedded in a majority of low-degree nodes. This makes the grid
**robust to random failure but clearly vulnerable to targeted attacks**,
especially attacks (or natural failures) concentrated on high-betweenness
nodes. Resilience planning should prioritize protecting these critical
intermediary nodes — with particular attention to sparsely connected rural
regions, where structural redundancy, and therefore resilience, is lowest.

**Future work**: extend the analysis to weighted graphs that incorporate
line capacity / power-flow attributes, and replicate the study for other
European transmission networks as equivalent open data becomes available.

## Repository contents

| File | Description |
|---|---|
| `power_grid_resilience_analysis.R` | Full analysis script: data loading, graph construction, centrality measures, community detection, and the attack-simulation experiments. |
| `README.md` | This document. |
| `LICENSE` | MIT License. |

## Tools

R, with the `igraph`, `sf`, `dplyr`, and `ggplot2` packages.

## References

- Cetinay, H., Devriendt, K., & Van Mieghem, P. (2018). *Nodal Vulnerability
  to Targeted Attacks in Power Grids.*
- Matke, C., Medjroubi, W., & Kleinhans, D. (2016). *SciGRID: An Open Source
  Reference Model for the European Transmission Network.*
- Medjroubi, W., et al. (2017). *Open Data in Power Grid Modelling: New
  Approaches Towards Transparent Grid Models.*
