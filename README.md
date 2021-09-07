# mito-network-sharing
Sharing of contents on mitochondrial encounter networks

Required: R with `igraph`, `brainGraph`, `ggplot2`, and `XML` libraries; `igraph` libraries for C

Wrapper scripts
------

`runcode.sh` calls R code to extract encounter networks from XML trajectory information in `Data/`, then C code to generate other networks for comparison and simulate the "bingo" game on these.

`plots.R` calls R code to produce summary plots of the results.

Code 
----

`bingo.c` is the workhorse C code for network generation and "bingo" simulation. This needs the `igraph` library.

`trajectory-analysis.R` extracts encounter networks from XML files and imposes any required restrictions (for example, truncting trajectory lengths)

The various `bingo-...-script.R` scripts produce visualisations of the different aspects of the simulations. 
