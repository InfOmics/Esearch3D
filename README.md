# Esearch3D: predictor of enhancer activity <img src="https://raw.githubusercontent.com/LucaGiudice/supplementary-files/main/D3SearchE_images/logo.png" width=250 align="right" style="border:4px solid black;" />

[![R-CMD-check](https://github.com/jokergoo/cola/workflows/R-CMD-check/badge.svg)](https://github.com/jokergoo/cola/actions)

## Features

1. It uses the expression levels of genes and a chromatin interaction network to impute the activity score of enhancers represented as nodes
2. It reverse engineers the flow of information modeled by enhancers that act as regulatory sources to increase the rate of transcription of their target genes
3. It leverages the relationship between the 3D organisation of chromatin and global gene expression
4. It represents a novel enhancer associated feature that can be used to predict enhancers
5. It provides a graphical user interface to interact with the network composed by genes, enhancers and intergenic loci associated to thier imputed activity score
6. It includes a parallelized version of the network-based algorithm called random walk with restart that works on sparce matrices to reduce the system usage

## Citation

Manindeer *, Giudice * et al. Esearch3D: Leveraging gene expression and chromatin architecture to predict intergenic enhancers

*equally contributed

## Install

*Esearch3D* is available on R, you can install it by:

```r
install.packages("devtools")
install.packages("htmltools", version="0.5.2", type="source")
devtools::install_github(
  repo="LucaGiudice/Esearch3D",
  ref = "main",
  dependencies = "Depends",
  upgrade = "always",
  quiet = F
  # type = "binary" # usable only in Windows OS
)
```

## Vignettes

There are the following vignettes:

1. [Quick Start with Dummy data](https://github.com/LucaGiudice/Esearch3D/blob/main/vignettes/Vignette.Rmd)
2. [Vignette with TSS data without enhancer annotation](https://github.com/LucaGiudice/D3SearchE/blob/main/vignettes/TSS_vignette.Rmd)
3. [Vignette with WG data without enhancer annotation](https://github.com/LucaGiudice/Esearch3D/blob/main/vignettes/WG_vignette.Rmd)
4. [Vignette with WG data with enhancer annotation and Machine Learning](https://github.com/LucaGiudice/Esearch3D/blob/main/vignettes/WG_ML_vignette.Rmd)


## Representation of chromatin data
<p align="center">
<img src="https://raw.githubusercontent.com/LucaGiudice/supplementary-files/main/D3SearchE_images/dati.png" width=750/>
</p>
Schematic of converting chromatin dynamics into networks:

1. Enhancers are localised and can influence the regulation of promoters directly and indirectly. 
2. 3C can capture chromatin interactions but only in a pairwise manner. 
3. Representation of chromatin fragments as nodes in a network and their interactions as edges preserves indirect associations between interacting chromatin.

## Representation of two step propagation

<p align="center">
<img src="https://raw.githubusercontent.com/LucaGiudice/supplementary-files/main/D3SearchE_images/two%20step.png" width=650 />
</p>

Schematic diagram of the network propagation used to impute activity values at intergenic nodes:

- A. Genes are mapped to nodes representing genic chromatin fragments. Each gene has an associated gene activity value determined by RNA-seq data. 
- B. Gene activity is propagated from gene nodes to genic chromatin nodes in propagation step one. Activity scores are then imputed in intergenic chromatin nodes by propagating the scores from genic chromatin nodes. 
- C. Ranking of non-genic nodes by the imputed activity score to identify high confidence enhancer nodes.

## Difference from two step propagation (aka multi gene propagation) and the single gene propagation

<p align="center">
<img src="https://raw.githubusercontent.com/LucaGiudice/supplementary-files/main/D3SearchE_images/single.png" width=650 />
</p>

A schematic of how multi-gene and single-gene propagation differ in the relative imputed activity scores. 
Multi-gene propagation highlights I6, an enhancer labelled node, with a higher IAS than single-gene propagation. 

### Usage

Few lines of code to run *D3SearchE* prediction:

```r
library(Esearch3D)

#Load and set up the example data ----
data("dummy_data_l")
#gene - fragment interaction network
gf_net=dummy_data_l$gf_net
#gene-fragment-fragment interaction network
ff_net=dummy_data_l$ff_net
#sample profile with starting values for genes and fragments
input_m=dummy_data_l$input_m
#length of chromosomes
chr_len=dummy_data_l$chr_len
#gene annotation
ann_net_b=dummy_data_l$ann_net_b

#Two step propagation -----
#Propagation over the gene-fragment network
gf_prop=rwr_OVprop(g=gf_net,input_m = input_m, no_cores=2, r=0.1)
#Propagation over the gene-fragment-fragment network
ff_prop=rwr_OVprop(g=ff_net,input_m = gf_prop, no_cores=2, r=0.8)

#Create igraph object with all the information included
net=create_net2plot(gf_net,input_m,gf_prop,ann_net_b,frag_pattern="F",ff_net,ff_prop)

#Start GUI
start_GUI(net, ann_net_b, chr_len, example=T)

#Single gene propagation -----
degree = 3
frag_pattern = "F"
gene_in=c("G1")

contrXgene_l=rwr_SGprop(gf_net, ff_net, gene_in, frag_pattern,
                        degree = degree, r1 = 0.1, r2 = 0.8, no_cores = 2)

#Create igraph object with all the information included
sff_prop=as.matrix(contrXgene_l$G1$contr_lxDest$ff_prop[,gene_in])
colnames(sff_prop)=gene_in
#Create igraph object with all the information included
net=create_net2plot(gf_net,input_m,gf_prop,ann_net_b,frag_pattern="F",ff_net,ff_prop)

#Start GUI
start_GUI(net, ann_net_b, chr_len, example=T)
```

### GUI

1. It allows to explore a sample's profile after a network-based propagation
2. It allows to investigate the imputed activity scores obtained by specific genes and their neighbourhood
3. It allows to download the propagated network and to import it in cytoscape

### Legend

<img src="https://github.com/LucaGiudice/supplementary-Simpati/blob/main/images/gui_overall.png" />

  1. Select nodes by chromosomes: allows to filter the network and to keep only those nodes (e.g. genes) that belong to a specific chromosome
  2. Select nodes by genome region: allows to filter the network and to keep only those nodes that belong to a specific genome region
  3. Select or type node + Distance: allows to visualize the neighbourhood of one specific node of interest
  4. Select by propagation ranges: allows to filter the subnetwork generated with a "search". It allows to keep only the nodes with a value that falls inside a specific range
  5. Scale colours: allows to scale the colors of the visualized nodes as if the propagation would have been applied only on them
  6. More/Less: allow to increase or decrease the size of the neighbourhood around a node of interest based on the distance (e.g. first degree neighbourhood, second degree ...)
  7. Open in cytoscape: allows to open the network in cytoscape
  8. Download: allows to download the image of the network visualized and created with the GUI

### Plots

Following plot shows you an example of how to interact with the GUI and its functionalities

<img src="https://github.com/LucaGiudice/supplementary-Simpati/blob/main/images/gui_scene.gif" />

## License

MIT @ Giudice Luca
