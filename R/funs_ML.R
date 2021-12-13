#' Zero - one standardization
#'
#' This function performs 0-1 standardization on a number vector
#'
#' @param x numeric vector
#' @return 0-1 standardized vector
#'
range01 <- function(x){(
  x-min(x))/(max(x)-min(x))
}

#' Maximum-based ranking normalization
#'
#' This function performs Maximum-based ranking normalization on a number vector
#'
#' @param x numeric vector
#' @return Normalized vector based on max rank
#'
rank_scale=function(v){
  y=rank(v,ties.method = "max")/length(v)
  return(y)
}

#' Get centrality measures about nodes
#'
#' This function determines the centrality measures of the genes and fragments included in the
#' CIN and that have a prior information about being or not being enhancer
#'
#' @param gf_net Edge list of the chromatin interaction network such that first column are genes and second column are "FragX" fragments
#' @param ff_net Edge list of the chromatin interaction network such that first and second column are "FragX" fragments
#' @param info 2-column dataframe, first column is character vector of node names (genes and fragments), second colum is numeric vector of number of enhancer annotations associated to a node
#' @return The info matrix with centrality measures associated to the nodes as new meta information
#' @import igraph
#' @import data.table
#' @import Matrix
#' @import scales
#' @export
get_centr_info=function(gf_net,ff_net,info){

  #Prepare info df
  info$Bin_Enh_Anno=as.numeric(info[,2]>=1)
  colnames(info)=c("ID","Num_Enh_Anno","Bin_Enh_Anno")
  info=as.data.frame(info)

  #Create the joint network
  data_net=ff_net

  #Format and filter from duplicated rows the network edgelist
  colnames(data_net)<-c("V1","V2")
  df=as.data.frame(data_net)
  df=df[!duplicated(t(apply(df,1,sort))),]
  net_m=as.matrix(df);rm(df);

  #Create igraph obj
  net1=igraph::graph_from_edgelist(net_m,directed = F)
  nodes_names_v=V(net1)$name
  frags_nodes_v=intersect(nodes_names_v,unique(info$ID))

  #Conpute centrality measures
  cat("Computing centrality measures of fragments in info matrix \n")
  deg_nodes_v=rank_scale(degree(net1,frags_nodes_v)[frags_nodes_v])
  eig_nodes_v=rank_scale(eigen_centrality(net1)$vector)
  eig_nodes_v=eig_nodes_v[match(frags_nodes_v,names(eig_nodes_v))]
  btw_nodes_v=rank_scale(estimate_betweenness(net1, v = frags_nodes_v,
                                              directed = F, cutoff=5))[frags_nodes_v]
  clo_nodes_v=rank_scale(estimate_closeness(net1, v = frags_nodes_v,
                                            cutoff=5))[frags_nodes_v]
  tra_nodes_v=rank_scale(transitivity(net1, type = c("local"),
                                      vids = frags_nodes_v, isolates = c("zero")))
  names(tra_nodes_v)=frags_nodes_v

  #Join all information
  cat("Merging fragments details and returning new info matrix \n")
  centr_nodes_v=data.frame(ID=names(deg_nodes_v), degree=deg_nodes_v, eigenc=eig_nodes_v,
                           betweenness=btw_nodes_v, closeness=clo_nodes_v, clustering=tra_nodes_v)
  info=merge(info,centr_nodes_v,by="ID",all.x = T)
  info[is.na(info)]=0
  return(info)
}

#' Get the optimal isolation parameters for the multi-gene two-step propagation
#'
#' This function determines the optimal values for the isolation parameter in the two
#' step propagation, one value for the propagation with the gene-fragment network and one
#' with the fragment-fragment network
#'
#' @param gf_net Edge list of the chromatin interaction network such that first column are genes and second column are "FragX" fragments
#' @param ff_net Edge list of the chromatin interaction network such that first and second column are "FragX" fragments
#' @param input_m numeric matrix of cell profile, one column of a cell, rows are genes or fragments, values are expression or IAS
#' @param info 2-column dataframe, first column is character vector of node names, second colum is numeric vector of number of enhancer annotations associated to a node
#' @param n_cores number of cores to run the function with parallel computing
#' @return List of two elements: combinations and best_comb
#' combinations contains all the pair of values tested
#' best_comb contains only the optimal two values of the isolation paramter for the first and second overall propagation
#' @import data.table
#' @import Matrix
#' @import scales
#' @import parallel
#' @import doParallel
#' @import Rfast
#' @import igraph
#' @import stats
#' @import foreach
#' @export
tuning_prop_vars = function(gf_net,ff_net,input_m,info=NULL,n_cores=2){
  if(!is.null(info)){

    if(!("Bin_Enh_Anno" %in% colnames(info))){
      #Prepare info df
      info$Bin_Enh_Anno=as.numeric(info[,2]>=1)
      colnames(info)=c("ID","Num_Enh_Anno","Bin_Enh_Anno")
      info=as.data.frame(info)
    }

    max_mean_prop=mean_prop_v=0
    r1=seq(0,0.2,0.04);stop_iters=c(seq(10,40,10),50,100,200)
    r2=seq(0.7,1,0.04);
    combs=expand.grid(r1, r2, stop_iters)
    colnames(combs)=c("r1","r2","stop_iters")

    for(k_comb in 1:nrow(combs)){
      cat("\n",k_comb,"on",nrow(combs),"\n")
      comb=combs[k_comb,]
      cat("testing this propagation setting:\n");print(comb);cat("\n")
      #Propagated for the network gene-fragment
      gf_prop=rwr_OVprop(g=gf_net,input_m = input_m,no_cores=n_cores,
                         r=comb$r1, stop_step = comb$stop_iters)
      #Propagated for the network fragment-fragment
      ff_prop=rwr_OVprop(g=ff_net,input_m = gf_prop,no_cores=n_cores,
                         r=comb$r2, stop_step = comb$stop_iters)

      indxs_enhan=match(info$ID[info$Bin_Enh_Anno==1],rownames(ff_prop))
      indxs_enhan=indxs_enhan[!is.na(indxs_enhan)]
      mean_prop=median(ff_prop[indxs_enhan,1])
      mean_prop_v=c(mean_prop_v,mean_prop)
    }
    mean_prop_v=mean_prop_v[-1]
    combs$prop_enh=mean_prop_v
    best_comb=combs[which.max(combs$prop_enh),]
    res=list(combinations=combs,best_comb=best_comb)
  }else{
    best_comb=data.frame(r1=0,r2=0.7,stop_iters=200)
    res=list(combinations=NULL,best_comb=best_comb)
  }

  return(res)
}

#' Merge propagation results with meta information
#'
#' This function merges propagation results with meta information providing the number of enhancer
#' annotions associated to the fragments in the network
#'
#' @param ff_prop Propagation result obtained with the second step of the multi-gene two-step propagation
#' @param info dataframe, first column is character vector of nodes, second colum is numeric vector of number of enhancer annotations associated to a node,
#' third column is the binarized version of the second column (1 the node has enhancer annotations, 0 otherwise), the
#' left columns are the centrality measures
#' @return New info dataframe containing the input information combined to the propagation values
#' @import data.table
#' @import Matrix
#' @import scales
#' @export
merge_prop_info = function(ff_prop,info){
  #Extact fragments's propagation values
  ov_prop=data.frame(as.matrix(ff_prop));ov_prop$ID=rownames(ov_prop)
  #Merge with annotation and centrality info
  info=merge(info,ov_prop,by="ID")
  colnames(info)[ncol(info)]="propagation"
  info$propagation=rank_scale(info$propagation)
  return(info)
}

#' Enhancer classifier based on centrality measures and propagation
#'
#' This function creates a classifier able to predict if a node is an enhancer or not based
#' on centrality measures and propagation values obtained with the multi-gene two-step propagation.
#' First, it creates the task, tunes a ranger random forest, trains the model and predict
#' a randomly selected testing set in cross validation.
#'
#' @param info Info dataframe containing centrality measures and propagation values about the fragments and their corresponding nodes
#' @param proj_name name fo the project as character string
#' @param annot Define the minimum number of annotations to define an enhancer (e.g. 5)
#' @param n_cores Define number of cores to use in order to parallelize the tune, train and test tasks
#' @param perf_measure Set performace measure (e.g. classif.acc)
#' @param tune_perf Performance to reach in the tuning phase (e.g 0.7 of accuracy)
#' @param n_folds Number of folds for the cross validation (e.g. 5)
#' @return List of 7 elements: proj_name, perfs_df, data_l, task, instance, lrn_tuned, lrn_tested, best_lrn
#' proj_name is the name of the project and run
#' perfs_df is the dataframe containing the classification performances
#' data_l is the list containing the original data processed by the classifier
#' task is the mlr3 object containing the task of classification
#' instance is the mlr3 object containing the tune, train and test parameters
#' lrn_tuned is the mlr3 object containing the classifier after the tuning
#' lrn_tested is the mlr3 object containing the classifier after the testing
#' best_lrn is the mlr3 object containing the classifier which best classifiers enhancers
#' @import data.table
#' @import Matrix
#' @import scales
#' @import future
#' @import mlr3verse
#' @import mlr3tuning
#' @import mlr3learners
#' @import mlr3viz
#' @export
enhancer_classifier = function(info,proj_name="enh_pred",annot=2,n_cores=2,
                               perf_measure="classif.acc",tune_perf=0.70,n_folds=5){

  if(n_cores>1){
    options(future.fork.multithreading.enable = TRUE)
    options(future.fork.enable = TRUE)
    future::plan("multiprocess")
  }
  #Keep feature columns: the binary annotation of enhancers and the centrailty measures
  keep_col=seq(3,ncol(info))
  #Keep no-enhancer nodes and enhancer with more than annot number of annotions
  keep_row=(info$Num_Enh_Anno==0 | info$Num_Enh_Anno>=annot)
  #Subset the info dataframe
  data=info[keep_row,keep_col];
  data[,1]=as.factor(data[,1])
  target_var=colnames(data)[1]
  pred=data[,target_var]
  #Create object list with all the data saved
  data_l=list(data=data,pred=as.numeric(as.character(pred)))
  #Create the task and set the measure for eval the performance
  task = TaskClassif$new(id = proj_name, backend = data, target = target_var)
  measure = msr(perf_measure)
  #Clean
  rm(keep_row,target_var,data,pred)

  #Learner setting ----

  lrn_2tune = lrn("classif.ranger", num.threads=n_cores)

  #Params to tune: lrn_2tune$param_set
  search_space = ps(
    alpha = p_dbl(lower = 0.5, upper = 0.8),
    min.node.size = p_int(lower = 40, upper = 60),
    mtry = p_int(lower = 4, upper = 6))

  #Set the cv appraoch for the tuning
  hout_4tune = rsmp("holdout", ratio=0.9, iters=1)

  #Create instance of tuning
  perf_reached_4tune = trm("perf_reached", level = tune_perf)
  instance = TuningInstanceSingleCrit$new(
    task = task,
    learner = lrn_2tune,
    resampling = hout_4tune,
    measure = measure,
    search_space = search_space,
    terminator = perf_reached_4tune
  )
  #Criteria of tuning
  tuner = tnr("grid_search", resolution = 5)

  #Start the tuning ----
  tuner$optimize(instance)
  #Results: instance$result_learner_param_vals, instance$result_y, instance$result

  #Test the resulting setting ----
  if(is.null(instance$result_y)){
    perf_measures=c("classif.recall","classif.precision",
                    "classif.sensitivity","classif.acc",
                    "classif.ppv","classif.npv")
    perfs=rep(0,length(perf_measures))
    perfs_df=data.frame(perf_measures,perfs,stringsAsFactors = F)

    res_ml=list(proj_name=proj_name,
                perfs_df=perfs_df,
                data_l=data_l,
                task=task,
                instance=NULL,lrn_tuned=NULL,lrn_tested=NULL,best_lrn=NULL)
  }else{
    #Set tuned learner
    lrn_tuned = lrn("classif.ranger", num.threads=n_cores,
                    alpha=instance$result$alpha,
                    mtry=instance$result$mtry,
                    min.node.size=instance$result$min.node.size)

    #Set cv approach :as.data.table(mlr_resamplings)
    cv_4test = rsmp("repeated_cv", folds=n_folds)
    cv_4test$instantiate(task);

    #Test tuned classifier
    lrn_tested = resample(task, lrn_tuned, cv_4test, store_models = T)
    cat("Final check up:\n");print(lrn_tested);cat("\n")
    cat("Overall performance:",lrn_tested$aggregate(measure),"\n")
    lrn_tested_perfs=lrn_tested$score(msr(perf_measure))

    #Retrieve the best tuned, trained and tested model
    best_model_indx = which.max(lrn_tested_perfs[[perf_measure]])
    best_lrn = lrn_tested$learners[[best_model_indx]]
    # best_model = best_lrn$model

    perf_measures=c("classif.recall","classif.precision",
                    "classif.sensitivity","classif.acc",
                    "classif.ppv","classif.npv")
    perfs=rep(0,length(perf_measures))
    perfs_df=data.frame(perf_measures,perfs,stringsAsFactors = F)
    perfs_df$perfs=apply(perfs_df,1,function(x){
      perf_measure=x[1]
      lrn_tested$aggregate(msr(perf_measure))
    })

    res_ml=list(proj_name=proj_name,
                perfs_df=perfs_df,
                data_l=data_l,
                task=task,
                instance=instance,
                lrn_tuned=lrn_tuned,
                lrn_tested=lrn_tested,
                best_lrn=best_lrn)
  }
  return(res_ml)
}

#' Explain the Enhancer classifier
#'
#' This function offers an explanation about how the classifier predicts the class of a
#' node/fragment and shows how each node's feature is contribuing in the classification.
#'
#' @param res_ml Output list produced and obtained by running the function: enhancer_classifier
#' @param n_cores Define number of cores to use in order to parallelize the tune, train and test tasks
#' @return List of 8 elements: task_dalex, lrn_dalex, ranger_exp, fi_ranger, fi_ranger_df, bd_ranger_enh, bd_ranger_no, best_lrn, pr_ranger
#' task_dalex is the dalex object containing the task of classification
#' lrn_dalex is the dalex object containing the learner
#' fi_ranger is the dalex object containing the trained and tested classifier
#' fi_ranger_df is the dataframe indicating how much each node's feature contributed in the classification
#' bd_ranger_enh is the dalex object containing the classification of an enhancer
#' bd_ranger_no is the dalex object containing the classification of a non enhancer
#' @import data.table
#' @import Matrix
#' @import scales
#' @import future
#' @import mlr3verse
#' @import mlr3tuning
#' @import mlr3learners
#' @import mlr3viz
#' @import DALEX
#' @import DALEXtra
#' @export
explain_classifier = function(res_ml,n_cores=2){

  proj_name=res_ml$proj_name
  data_l=res_ml$data_l
  instance=res_ml$instance

  task_dalex = TaskClassif$new(id = "DALEX_proj", backend = data_l$data,
                               target = colnames(data_l$data)[1])
  lrn_dalex = lrn("classif.ranger", num.threads=n_cores, predict_type = "prob",
                  alpha=instance$result$alpha,
                  mtry=instance$result$mtry,
                  min.node.size=instance$result$min.node.size)
  lrn_dalex$train(task_dalex)

  ranger_exp <- explain_mlr3(lrn_dalex,
                             data = data_l$data[,-1],
                             y = data_l$pred,
                             label    = "Ranger RF",
                             colorize = FALSE,
                             verbose = TRUE)

  #Importance features
  fi_ranger <- model_parts(ranger_exp)
  fi_ranger_df=as.data.frame(fi_ranger)
  fi_ranger_df=stats::aggregate(dropout_loss ~ label+variable, fi_ranger_df, mean)
  fi_ranger_df=fi_ranger_df[order(fi_ranger_df$variable,decreasing = T),]
  fi_ranger_df=fi_ranger_df[,c(2,3)]

  best_fs=fi_ranger_df$variable[order(fi_ranger_df$dropout_loss,decreasing = T)]
  best_fs=best_fs[-grep("_",best_fs)]
  best_fs=as.character(best_fs)

  ord_data=data_l$data[order(data_l$data[,best_fs[1]],
                             data_l$data[,best_fs[2]],
                             data_l$data[,best_fs[3]],decreasing=TRUE),]

  first_enh=ord_data[which(ord_data$Bin_Enh_Anno==1)[1],]

  ord_data=data_l$data[order(data_l$data[,best_fs[1]],
                             data_l$data[,best_fs[2]],
                             data_l$data[,best_fs[3]],decreasing=FALSE),]

  first_no_enh=ord_data[which(ord_data$Bin_Enh_Anno==0)[1],]

  bd_ranger_enh <- predict_parts(ranger_exp, new_observation = first_enh)

  bd_ranger_no <- predict_parts(ranger_exp, new_observation = first_no_enh)


  #Feature profile with respect the profile of predictions
  pr_ranger = model_profile(ranger_exp)


  res_dalex=list(task_dalex=task_dalex, lrn_dalex=lrn_dalex, ranger_exp=ranger_exp,
                 fi_ranger=fi_ranger, fi_ranger_df=fi_ranger_df, bd_ranger_enh=bd_ranger_enh,
                 bd_ranger_no=bd_ranger_no, pr_ranger=pr_ranger)

  return(res_dalex)
}

#' Make the first letter a capital letter
#'
#' This function makes the first letter a capital letter while all the other letters low
#'
#' @param x character string
#' @return Make the first letter a capital letter
#' @export
#'
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
