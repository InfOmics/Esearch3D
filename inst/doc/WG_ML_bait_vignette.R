## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------

#Clean workspace and memory ----
rm(list=ls())
gc()

#Set working directory ----
gps0=getwd()
gps0=paste(gps0,"/%s",sep="")
rootDir=gps0
setwd(gsub("%s","",rootDir))

#Load libraries ----
suppressWarnings(suppressMessages(
  library("Esearch3D", quietly = T)
  )
)

#Set variables ----
#Set seed to get always the same results out of this vignette
set.seed(8)
#Set number of cores to parallelize the tasks
n_cores=16

## ----loading of wg------------------------------------------------------------

#Load and set up the example data ----
data("wg_ann_data_l")
#gene - fragment interaction network generated from DNase_Prop1_mESC_TSS interactions data
gf_net=wg_ann_data_l$gf_net
#gene-fragment-fragment interaction network generated from mESC_DNase_Net interactions data
ff_net=wg_ann_data_l$ff_net
#sample profile with starting values for genes and fragments generated from mESC_bin_matrix_Prop1
input_m=wg_ann_data_l$input_m

## ----loading related enhancer annotation data---------------------------------

#info dataframe containg for each node and fragment the number of enhancer annotations associated to it
info=wg_ann_data_l$info
info=info[info$Type=="Bait",]
info=info[,-3]
head(info)
#dataframe containg mmu nomenclature about genes
mouse_db=wg_ann_data_l$mouse_db;head(mouse_db);

## ----calculation of centrality measures---------------------------------------
#Process info matrix to get fragments extra information ----
info=get_centr_info(gf_net,ff_net,info);head(info)


## ----nomenclature change------------------------------------------------------

#Convert ENSG to SYMB ----
gf_net[,1]=mapvalues(gf_net[,1],from = mouse_db$ENS, to = mouse_db$Symbol, warn_missing = F)
gf_net[,2]=mapvalues(gf_net[,2],from = mouse_db$ENS, to = mouse_db$Symbol, warn_missing = F)
ens_row=rownames(input_m);sym_row=mapvalues(ens_row,from = mouse_db$ENS, to = mouse_db$Symbol, warn_missing = F);rownames(input_m)=sym_row


## ----two step propagation with automatic tuning of the isolation parameter----

#Tuning of propagation setting -----
res_tuning=tuning_prop_vars(gf_net,ff_net,input_m,n_cores=n_cores);head(res_tuning)

#Two step propagation -----
#Propagated for the network gene-fragment
gf_prop=rwr_OVprop(g=gf_net,input_m = input_m, no_cores=n_cores, 
                   r=res_tuning$best_comb$r1,
                   stop_step = res_tuning$best_comb$stop_iters)
#Propagated for the network fragment-fragment
ff_prop=rwr_OVprop(g=ff_net,input_m = gf_prop, no_cores=n_cores, 
                   r=res_tuning$best_comb$r2,
                   stop_step = res_tuning$best_comb$stop_iters)


## ----enhancer classification, results = "hide"--------------------------------

#Merge propagation results with meta information about enhancer annotations
info=merge_prop_info(ff_prop,info)

#Build enhancer classifier
res_ml=enhancer_classifier(info, n_cores=n_cores)

## ----model explainer----------------------------------------------------------
#Build explainer
res_dalex=explain_classifier(res_ml, n_cores=n_cores)

#Plot results
res_dalex$fi_ranger_df

plot(res_dalex$fi_ranger)

plot(res_dalex$bd_ranger_enh)

plot(res_dalex$bd_ranger_no)

plot(res_dalex$pr_ranger)

