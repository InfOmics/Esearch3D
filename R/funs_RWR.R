#' Harmonizes the matrix of cell expression profiles with respect the CIN network
#'
#' Remove all the genes that don't match between the cell expression profiles and the nodes inside the network
#'
#' Given an expression matrix of one or multiple profiles at the columns and genes located at the rows:
#' The method keeps only the genes that exist also in the network (CIN)
#' This because only existing genes in the cell profile can give information to the corresponding nodes in the network and allow the propagation
#'
#' @param m_profiles numeric matrix of profiles, rows are genes, column is a cell or sample, value is numeric indicating a gene expression in a cell column and profile
#' @param net adjacency matrix of the CIN network (it has to contain nodes with the same names of genes characterizing the input cell profiles)
#' @return matrix of profiles without the rows that are not represented by a node in the network
#' @export
#'
complete_m = function(m_profiles,net){
  #bind the genes which are not expressed in the input matrix of patient profiles but they exist in the network
  #Recover patient genes and network genes
  pat_genes=rownames(m_profiles)
  network_genes=rownames(net)

  #Find genes which exist in the patients but do no exist in the network
  network_genes_2rem=setdiff(pat_genes,network_genes)
  #Remove
  if(length(network_genes_2rem)!=0){
    m_profiles2=as.matrix(m_profiles[-match(network_genes_2rem,pat_genes),])
    colnames(m_profiles2)=colnames(m_profiles)
    m_profiles=m_profiles2
    pat_genes=rownames(m_profiles)
  }
  #Find genes which exist in the network but do no exist in the patients
  network_genes_2add=setdiff(network_genes,pat_genes)
  #Normalize
  if(length(network_genes_2add)!=0){
    m2bind=matrix(0,length(network_genes_2add),length(colnames(m_profiles)))
    m_profiles=rbind(m_profiles,m2bind)
    rownames(m_profiles)=c(pat_genes,network_genes_2add)
  }
  #Finish
  return(m_profiles)
}

#' Network based propagation with random walk for the multi-gene two-step propagation
#'
#' This function performs propagation of each individual cell expression profile associated to a column of an expression matrix.
#' The cell expression profile is characterized by the expression of genes at the rows.
#' The values of the genes indicate their expression in a specific column and cell profile
#' The propagation considers each cell profile at a time. It maps the expression values to the corresponding genes represented as nodes in the network.
#' It diffuses the gene expression values. Each node in the CIN network composed by genes, genic fragments and intergenic fragments gets a propagation score which
#' measures its activity. The imputed activity score (IAS) obtained by genes and fragments replaces their original value in the cell expression profile.
#' The result is a new matrix of cell profiles. Each cell profile is now described with genes and also fragments who own a IAS.
#'
#'
#' @param g igraph, edge list or adjacency matrix of the chromatin interaction network (CIN)
#' @param input_m numeric matrix of cell profiles, columns are cells, rows are genes or fragments, value is expression of IAS
#' @param norm Default column, options: row or column or laplacian, it sets the graph normalization performed on the network
#' @param no_cores Default 2, An integer value greater than 0 indicating the number of cores to use for computing in parallel the task
#' @param r Default 0.8, A double value lower than 1 indicating the percentage of information that a gene keeps (0.8 is 80 percentage)
#' @param stop_step Default 200, An integer value greater than 0 indicating the number of iterations of the propagation
#' @param stop_delta Default 1e-06, A double value lower than 0 indicating the threshold under which all imputed propagation values are set 0
#' @param keep_zero Default FALSE, Boolean: TRUE if to keep profiles without information
#' @return matrix of propagated profiles with imputed activity scores for the genes and fragments in the network and that replace the original expression profiles
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
#'
rwr_OVprop=function(g, input_m, norm = "column", no_cores=2, r = 0.8, stop_step=200, stop_delta = 1e-06, keep_zero=F){
  cat("**Performing two step propagation \n")
  start_time <- Sys.time()
  #Set input variables ----
  if(class(g)[[1]]=="igraph"){
    g=igraph::simplify(g)
    Isolated=which(igraph::degree(g)==0)
    g = igraph::delete.vertices(g, Isolated)
    adjM=igraph::get.adjacency(g, type = "both", attr = NULL, edges = F,
                               names = T, sparse = getIgraphOpt("sparsematrices"))

  }
  if(class(g)[[1]]=="matrix" | class(g)[[1]]=="data.frame"){
    if(ncol(g)==2){
      g=igraph::graph_from_edgelist(g, directed = F)
      g=igraph::simplify(g)
      Isolated = which(igraph::degree(g)==0)
      g = igraph::delete.vertices(g, Isolated)
      adjM=igraph::as_adjacency_matrix(g, type = "both", attr = NULL,
                                       names = T, sparse = igraph::getIgraphOpt("sparsematrices"))
    }else{
      adjM=g
    }
  }

  data=as.matrix(complete_m(input_m,adjM))
  zero_data=data[,colSums(data)==0]
  data1=as.matrix(data[,colSums(data)!=0])
  colnames(data1)=colnames(data)[colSums(data)!=0]
  data=data1;rm(data1)
  cnames=colnames(data)
  rnames=rownames(data)

  #Apply network normalization ----
  A <- adjM != 0
  if(norm == "row"){
    D <- Matrix::Diagonal(x = (Matrix::rowSums(A))^(-1))
    nadjM <- adjM %*% D
  }else if(norm == "column"){
    D <- Matrix::Diagonal(x = (Matrix::colSums(A))^(-1))
    nadjM <- D %*% adjM
  }else if(norm == "laplacian"){
    D <- Matrix::Diagonal(x = (Matrix::colSums(A))^(-0.5))
    nadjM <- D %*% adjM %*% D
  }else{
    nadjM <- adjM
  }

  ind <- match(rownames(data), rownames(adjM))
  P0matrix <- matrix(0, nrow = nrow(nadjM), ncol = ncol(data))
  P0matrix[ind[!is.na(ind)], ] <- as.matrix(data[!is.na(ind), ])

  #Prepare matrix for the propagation ----
  rownames(P0matrix) <- rownames(adjM)
  colnames(P0matrix) <- cnames
  P0matrix <- Matrix::Matrix(P0matrix, sparse = T)

  #Set parallel and apply network-based propagation ----
  if(.Platform$OS.type == "unix") {
    cat("Parallel for linux \n")
    cl <- makeCluster(no_cores,type="FORK");
  } else {
    cat("Parallel for windows \n")
    cl <- makeCluster(no_cores);
  }
  doParallel::registerDoParallel(cl)

  ind=seq(1,ncol(P0matrix))
  sets=split(ind, ceiling(seq_along(ind)/round(length(ind)/no_cores)))
  P0_l=lapply(sets,function(x){Matrix::Matrix(as.matrix(P0matrix[,x]), sparse = T)})


  PTmatrix=foreach(k = 1:length(P0_l),.inorder = T,.combine = c("cbind"),.packages=c("Matrix")) %dopar% {

    P0matrix1=P0_l[[k]]
    for(j in 1:ncol(P0matrix1)){
      P0 <-P0matrix1[, j]; step <- 0; delta <- 1; PT <- P0;
      while(step <= stop_step){
        PX <- (1 - r) * nadjM %*% PT + r * P0
        delta <- sum(abs(PX - PT))
        PT <- PX
        step <- step + 1
      }
      PT=as.matrix(PT)
      P0matrix1[,j]=PT
    }
    return(P0matrix1)

  }

  #Finish the matrix, close parallel and return the result ----
  PTmatrix[PTmatrix < stop_delta]= 0
  colnames(PTmatrix)=colnames(P0matrix)
  parallel::stopCluster(cl)

  if(keep_zero == T){
    tmp=merge(PTmatrix,zero_data,by="row.names")
    rownames(tmp)=tmp[,1];tmp=tmp[,-1]
    PTmatrix=tmp
  }

  end_time <- Sys.time()
  running_time=as.numeric(difftime(end_time, start_time, units = "mins")[[1]])
  cat("running time:",round(running_time,digits = 4),"mins", "\n")

  return(PTmatrix)
}

#' Convert matrix to diagonal matrix of row elements
#'
#' Convert matrix to diagonal matrix of row elements
#'
#' @param input_m Matrix
#' @return Return a diagonal matrix of the row elements
#'
diag_prof = function(input_m){
  input_m_sg=diag(as.vector(input_m), nrow(input_m), nrow(input_m))
  rownames(input_m_sg)=colnames(input_m_sg)=rownames(input_m)
  return(input_m_sg)
}

#' Given the seed nodes of the propagation, it provides the target nodes of each seed
#'
#' It takes the result of a run of network-based propagation and determines which and how much information
#' the target nodes received from a specific seed node
#'
#' @param ff_prop numeric matrix of propagated profiles: seed nodes at the columns and target nodes at the rows
#' @param no_cores Default 1, An integer value greater than 0 indicating the number of cores to use for computing in parallel this task
#' @return A list with contr_prop_l: For each seed node as element of the list, the vector of the targets nodes which receive information from it
#' Plus, also the original propagation profiles of the seed nodes.
#' @import data.table
#' @import Matrix
#' @import scales
#' @import parallel
#' @import doParallel
#' @import Rfast
#' @import igraph
#' @import stats
#' @import foreach
#'
get_contrib_propXorig = function(ff_prop,no_cores=1){
  ff_prop=as.matrix(ff_prop)
  ff_prop=ff_prop[order(rownames(ff_prop),decreasing = F),]
  ff_prop=t(ff_prop)
  ff_prop=ff_prop[rowsums(ff_prop,parallel = T)!=0,]
  mult_factor=100/rowsums(ff_prop,parallel = T)
  contr_prop=t(ff_prop*mult_factor)

  contr_prop_l=lapply(seq_len(ncol(contr_prop)), function(i) contr_prop[,i])
  names(contr_prop_l)=colnames(contr_prop)

  if(no_cores!=1){
    require(parallel)
    #Set parallel and apply network-based propagation ----
    if(.Platform$OS.type == "unix") {
      cat("Parallel for linux \n")
      cl <- makeCluster(no_cores,type="FORK");
    } else {
      cat("Parallel for windows \n")
      cl <- makeCluster(no_cores);
    }
    doParallel::registerDoParallel(cl)

    contr_prop_l=parLapply(cl=cl,X=contr_prop_l,fun=function(x){
      x=x[x>=1]
      x=x[order(x,decreasing = T)]
      return(x)
    })
    parallel::stopCluster(cl)
  }else{
    contr_prop_l=lapply(contr_prop_l,function(x){
      x=x[x>=1]
      x=x[order(x,decreasing = T)]
      return(x)
      })
  }

  res_contr=list(contr_prop_l=contr_prop_l,ff_prop=ff_prop)
  return(res_contr)
}

#' Given the target nodes of the propagation, it provides the seed nodes of each target
#'
#' It takes the result of a run of network-based propagation and determines which and how much information
#' the seed nodes gave to a specific target node
#'
#' @param ff_prop numeric matrix of propagated profiles: seed nodes at the columns and target nodes at the rows
#' @param no_cores Default 1, An integer value greater than 0 indicating the number of cores to use for computing in parallel this task
#' @return A list with contr_prop_l: For each target node as element of the list, the vector of the seed nodes which gave information to it
#' Plus, also the original propagation profiles of the target nodes.
#' @import data.table
#' @import Matrix
#' @import scales
#' @import parallel
#' @import doParallel
#' @import Rfast
#' @import igraph
#' @import stats
#' @import foreach
#'
get_contrib_propXdest = function(ff_prop,no_cores=1){
  ff_prop=as.matrix(ff_prop)
  ff_prop=ff_prop[order(rownames(ff_prop),decreasing = F),]
  mult_factor=100/rowsums(ff_prop,parallel = T)
  contr_prop=t(ff_prop*mult_factor)

  contr_prop_l=lapply(seq_len(ncol(contr_prop)), function(i) contr_prop[,i])
  names(contr_prop_l)=colnames(contr_prop)

  if(no_cores!=1){
    require(parallel)
    #Set parallel and apply network-based propagation ----
    if(.Platform$OS.type == "unix") {
      cat("Parallel for linux \n")
      cl <- makeCluster(no_cores,type="FORK");
    } else {
      cat("Parallel for windows \n")
      cl <- makeCluster(no_cores);
    }
    doParallel::registerDoParallel(cl)

    contr_prop_l=parLapply(cl=cl,X=contr_prop_l,fun=function(x){
      x=x[x>=1]
      x=x[order(x,decreasing = T)]
      return(x)
      })
    parallel::stopCluster(cl)
  }else{
    contr_prop_l=lapply(contr_prop_l,function(x){
      x=x[x>=1]
      x=x[order(x,decreasing = T)]
      return(x)
      })
  }

  res_contr=list(contr_prop_l=contr_prop_l,ff_prop=ff_prop)
  return(res_contr)
}

#' Determines the shortest distances from all the genes to the non genic fragments
#'
#' Given the gene-genic fragment network and the fragment-fragment network, it computes the distances from the genes to the non genic fragments
#'
#' @param gf_net Edge list representing a network such that first column are genes and second column are "FragX" fragments
#' @param ff_net Edge list representing a network such that first and second column are "FragX" fragments
#' @param frag_pattern String pattern to identify a fragment by its nomenclature in the edge list (e.g. "Frag")
#' @return Distance matrix of genes and non-genic fragments
#' @import cppRouting
#' @import Rfast
#' @import igraph
#'
get_dist_m = function(gf_net,ff_net,frag_pattern="Frag"){
  #Pattern to recognise genes in the network
  frag_pattern=paste("^",frag_pattern,sep="")
  #Reconstruct the original network
  net=rbind(gf_net,ff_net);
  #Get the genes or genic fragments
  gs=unique(as.vector(gf_net))
  gs=grep(frag_pattern,gs,value = T,invert=T)
  #Get the non genic fragments
  ff=unique(as.vector(ff_net))
  #Get the network version
  net=cbind(net,1)
  net=as.data.frame(net,stringsAsFactors = F)
  net$V3=as.numeric(net$V3)
  colnames(net)=c("from","to","cost")
  g=makegraph(net,directed = F)
  #Compute the distances from all genes and fragments
  dist_m=get_distance_matrix(g,from=gs,to=ff,allcores = T)
  return(dist_m)
}

#' Supportive function to infer the networks and the profiles needed in the single gene propagation
#'
#' @param gf_net Edge list representing a network such that first column are genes and second column are "FragX" fragments
#' @param ff_net Edge list representing a network such that first and second column are "FragX" fragments
#' @param input_m Numeric matrix of cell profiles
#' @param dist_m Distance matrix computed between the genes and the non-genic fragments in the network
#' @param gene_in Character vector of the genes of interest
#' @param frag_pattern Character string to identify the fragments in the edge list (e.g. "FragX")
#' @param degree Default 2, the degree of distance which is used to select the most important fragments connected to the genes of interest
#' @return A list with elements describing the input data which are by the method
#' @import data.table
#' @import Matrix
#' @import scales
#' @import parallel
#' @import doParallel
#' @import Rfast
#' @import igraph
#' @import stats
#' @import foreach
#' @import cppRouting
#'
get_single_prop_data = function(gf_net,ff_net,input_m,dist_m,
                                gene_in,
                                frag_pattern="Frag",
                                degree=2){
  #Pattern to recognise genes in the network
  frag_pattern=paste("^",frag_pattern,sep="")
  #Reconstruct the original network
  net=rbind(gf_net,ff_net);
  #Get the igraph version
  ig=graph_from_edgelist(net, directed = F)
  #Get the genes or genic fragments
  gs=unique(as.vector(gf_net))
  gs=grep(frag_pattern,gs,value = T,invert=T)
  #Keep only the genes of interest which are in the network
  gene_in=intersect(gene_in,gs)
  #Get the non genic fragments
  ff=unique(as.vector(ff_net))
  #Get the network version
  net=cbind(net,1)
  net=as.data.frame(net,stringsAsFactors = F)
  net$V3=as.numeric(net$V3)
  colnames(net)=c("from","to","cost")
  g=makegraph(net,directed = F)

  #Get the fragments that are close to the genes of interest
  if(length(gene_in)>1){
    frag_int=colnames(dist_m)[colSums(dist_m[gene_in,]<=degree)>0]
  }else{
    frag_int=colnames(dist_m)[dist_m[gene_in,]<=degree]
  }

  #Get all the genes that are close to the fragments
  frags_int_df=dist_m[,frag_int]
  index <- which(frags_int_df<=degree, arr.ind=TRUE)
  genes_v=rownames(frags_int_df)[index[,1]]
  frags_v=colnames(frags_int_df)[index[,2]]

  #Get the pathes among them
  pathes=get_path_pair(g,from=genes_v,to=frags_v,long=T)
  #Get all the nodes in the pathes
  nodes_in=unique(pathes$node)
  #Get the subnetwork
  sub_ig=induced_subgraph(ig, nodes_in)

  #Split in the two parts
  subnet=as_edgelist(sub_ig)
  indxs_gf=unique(c(grep(frag_pattern,subnet[,1],invert=T),
                    grep(frag_pattern,subnet[,2],invert=T)))
  indxs_ff=setdiff(seq(1,nrow(subnet)),indxs_gf)
  gf_subnet=subnet[indxs_gf,]
  ff_subnet=subnet[indxs_ff,]

  #Prepare the input matrix of profiles
  s_input_m=diag_prof(input_m)

  data_l=list(s_input_m=s_input_m,
              gf_net=gf_subnet,
              ff_net=ff_subnet)
  return(data_l)
}

#' Single gene network based propagation based on random walk with restart
#'
#' Performs the single gene propagation based on Random walk with restart.
#' This operation is suggested to apply only after having performed the standard propagation based on all the genes with the two steps.
#' First, overall propagation to diffuse expression of all genes included in a cell profile to genic fragments.
#' Second, overall propagation to diffuse the propagated expression of all genes from genic fragments to intergenic.
#'
#' This function performs the propagation of individual genes of interest belonging to a cell's expression profile.
#' It does not perform the standard propagation using all the genes of a cell profile.
#' It then returns how much each gene of interest contributed to give information to the fragments.
#' It then returns how much each fragment received information from the genes of interest.
#' This function helps to understand the contribution of the individual genes in the two-step standard propagation.
#'
#' @param gf_net Edge list of the chromatin interaction network such that first column are genes and second column are "FragX" fragments
#' @param ff_net Edge list of the chromatin interaction network such that first and second column are "FragX" fragments
#' @param gene_in Character vector of the genes of interest
#' @param input_m numeric matrix of cell profile, one column of a cell, rows are genes or fragments, values are expression or IAS
#' @param frag_pattern Character string to identify the fragments in the edge list (e.g. "FragX")
#' @param out_rda Default: sg_prop.rda, Character string to define the output rda file which will contain all the results of the analysis and intermediate products
#' @param degree Default 4, the degree of distance which is used to select the most important fragments connected to the genes of interest (e.g. degree 2 means G --> F1 --> F2)
#' @param r1 Default 0.1, A double value lower than 1 indicating the percentage of information that a gene keeps (0.1 is 10 percentage)
#' @param r2 Default 0.8, A double value lower than 1 indicating the percentage of information that a fragment keeps (0.8 is 80 percentage)
#' @param no_cores Default 2, An integer value greater than 0 indicating the number of cores to use for parallel computing
#' @return A list single_gene_prop composed of two sublists
#' contr_lxDest: For each target node as element of the list, the vector of the seed nodes which gave information to it
#' contr_lxOrig: For each seed node as element of the list, the vector of the targets nodes which receive information from it
#' @import data.table
#' @import Matrix
#' @import scales
#' @import parallel
#' @import doParallel
#' @import Rfast
#' @import igraph
#' @import stats
#' @import foreach
#' @import cppRouting
#' @export
#'
rwr_SGprop = function(gf_net,ff_net,gene_in,input_m,frag_pattern="frag",out_rda="sg_prop.rda",
                      degree=4,r1=0.1,r2=0.8,no_cores=2){
  cat("\n###################################\n")
  cat("Starting single gene propagation \n")
  single_gene_prop = list()
  #Compute distance matrix
  cat("*Computing distance matrix\n")
  dist_m=get_dist_m(gf_net, ff_net, frag_pattern = frag_pattern)
  #Consider only the genes in the network
  gs=as.vector(gf_net)
  gene_in=intersect(gene_in,gs)
  for(n_g in 1:length(gene_in)){
    g_in=gene_in[n_g]
    cat("*Obtaining subnetwork of ",g_in,"\n")
    sing_data=get_single_prop_data(gf_net, ff_net, input_m, dist_m,
                                   gene_in = g_in,
                                   degree = degree,
                                   frag_pattern = frag_pattern)
    cat("*Perfoming the propagation \n")
    sgf_prop=rwr_OVprop(g=sing_data$gf_net,input_m = sing_data$s_input_m,
                      no_cores=no_cores, r=r1)

    sff_prop=rwr_OVprop(g=sing_data$ff_net,input_m = sgf_prop,
                      no_cores=no_cores, r=r2)
    cat("*Assessing the contribution of each gene and fragment on the\n final propagation scores \n")
    #Get how much each starting informative node contributed to each node of the net
    contr_lxOrig=get_contrib_propXorig(sff_prop)
    #Get how much each node of the network got infromation from the starting nodes
    contr_lxDest=get_contrib_propXdest(sff_prop)

    cat("*Saving results")
    single_gene_prop[[g_in]]=list(contr_lxOrig=contr_lxOrig,
                                  contr_lxDest=contr_lxDest)
    save(single_gene_prop,file=out_rda)
    gc()
  }
  return(single_gene_prop)
}

