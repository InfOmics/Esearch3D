#' Create an igraph object of the chromatin interaction network (CIN) for visualization purposes
#'
#' Given the input and output of the two step network based propagation, it assembles an igraph object of the CIN with all
#' the information included in order to be visualized or analysed
#'
#' @param g_net Edge list of the chromatin interaction network such that first column are genes and second column are "FragX" fragments
#' @param input_m numeric matrix of a cell expression profile before the propagation
#' @param gf_prop numeric matrix of a cell profile after the first step of the propagation applied with only the gene-genic fragment component of the CIN
#' @param ann_net_b matrix, for each row presents the gene identifier, the chromosome in which the gene is, the starting end ending position in the sequence.
#' @param frag_pattern string, initial character of the fragments name (e.g. "F" or "Frag")
#' @param ff_net Edge list of the chromatin interaction network such that first and second column are "FragX" fragments
#' @param ff_prop numeric matrix of a cell profile after the second step of the propagation applied with the fragment-fragment component of the CIN
#' @return igraph object
#' @import igraph
#' @export
#'
create_net2plot = function(g_net, input_m, gf_prop, ann_net_b, frag_pattern="F",
                          ff_net = NULL,ff_prop = NULL){
  #Create the joint network
  if(is.null(ff_net)){
    data_net = g_net
  }
  else{
    data_net=rbind(g_net, ff_net)
    frag_pattern=paste("^",frag_pattern,sep="")
  }

  #Format and filter from duplicated rows the network edgelist
  colnames(data_net)<-c("V1","V2")
  df=as.data.frame(data_net)
  df=df[!duplicated(t(apply(df,1,sort))),]
  net_m=as.matrix(df);rm(df)

  #Create igraph obj
  net1=igraph::graph_from_edgelist(net_m,directed = F)
  namesvector<-igraph::get.vertex.attribute(net1,"name")

  #Merge the two propagations
  if(is.null(ff_net)){
    merged_prop = gf_prop
  }
  else{
    check=!(rownames(gf_prop) %in% rownames(ff_prop))
    merged_prop=rbind(gf_prop[check,],ff_prop)
    rownames(merged_prop)[1:(sum(check))]=rownames(gf_prop)[check]
  }


  #Create propagation vector for network nodes: df with name | value : str and int
  prop=data.frame(name=rownames(merged_prop),value=merged_prop[,1],stringsAsFactors = F)
  rownames(prop)=seq(1,nrow(prop))
  miss_nodes=setdiff(as.vector(data_net),prop$name)
  if(length(miss_nodes)!=0){
    prop2=data.frame(name=miss_nodes,value=min(prop$value))
    prop=rbind(prop,prop2)
  }
  prop$value=round(prop$value,digits = 2)

  #Create starting vector for network genes: df with name | expr: str and others "-"
  starting_prop=data.frame(name=rownames(input_m),expr=input_m[,1])
  if(length(setdiff(as.vector(data_net),starting_prop$name)) > 0){
    starting_prop2=data.frame(name=setdiff(as.vector(data_net),starting_prop$name),expr="-")
    starting_prop=rbind(starting_prop,starting_prop2)
  }
  rownames(starting_prop)=seq(1,nrow(starting_prop))

  #Map of the edge colors
  edge_group=rep("orange",nrow(net_m))
  if(!is.null(ff_net)){
    gene_indxs=grep(frag_pattern,V(net1)$name,invert = T)
  }else {
     gene_indxs = 1:length(V(net1))
  }
  net1=igraph::set.vertex.attribute(net1,name="shape",value="circle",index = igraph::V(net1)[gene_indxs])
  net1=igraph::set.edge.attribute(net1, name="color", value=edge_group)
  if(!is.null(ff_net)){
    edge_group[grep(frag_pattern,net_m[,1])] = "purple"
    edge_group[grep(frag_pattern,net_m[,2])] = "purple"
    net1=igraph::set.edge.attribute(net1, name="color", value=edge_group)
    #set genes as circles and fragments as squares and edges' colors
    frag_indxs=grep(frag_pattern,V(net1)$name,invert = F)
    net1=igraph::set.vertex.attribute(net1,name="shape",value="square",index = igraph::V(net1)[frag_indxs])
  }

  #assignment to network
  color_info=combine_netWprop(net1, prop, starting_prop)
  V(net1)$color=color_info[,3]
  V(net1)$propagation=color_info[,2]
  V(net1)$expr=color_info[,4]

  net1=toolkitf(net1,ann_net_b,frag_pattern,ff_net)

  return(net1)
}

#' Downloads the genes annotations for human and mouse
#'
#' @param organism ,string representing the orgnanism for which downloading the genes annotations
#' @param gene_names, vector of strings, gene names to search annotation for
#' @return list with the gene id, chromosome, start and end positions.
#' @import biomaRt
#' @export
#'
download_genes_annotation = function(organism = "human",gene_names){
  if(organism != "human" && organism != "mouse"){
    cat("Error: you are asking for an unsupported organism")
    return(NA)
  }
  if(length(gene_names) == 0){
    cat("Error: no gene names provided")
    return(NA)
  }
  dataset=""
  if(organism == "human"){
    dataset = "hsapiens_gene_ensembl"
    gene_id = "hgnc_symbol"
  }
  if(organism == "mouse"){
    dataset = "mmusculus_gene_ensembl"
    gene_id = "mgi_symbol"
  }

  ensembl <- useEnsembl(biomart = "ensembl", dataset = dataset)
  a = getBM(attributes=c(gene_id,"chromosome_name","start_position","end_position"),
    filters=c("chromosome_name",gene_id),values=list(c(1:19,"X","Y"),gene_names),mart=ensembl,useCache = FALSE)
  a = a[a[,1] != "",]
  a[,2] = paste("chr",a[,2],sep="")
  colnames(a) = c("ID","chr","start","end")
  differece_length = length(gene_names) - dim(a)[1]
  if(differece_length > 0 ){
    cat(differece_length, " gene(s) annotation(s) has/have not been found\n",sep="")
  }
  return(a)
}

#' Combine network data with propagation values
#'
#' Support function for the network visualization
#'
#' @param net igraph object
#' @param propag propagation matrix
#' @param totaldataframe metainformation matrix
#' @import igraph
#' @import grDevices
#'
combine_netWprop = function(net, propag, totaldataframe) {
  totaldataframe[,2]=as.character(totaldataframe[,2])

  dtnames<-data.frame(names = names(V(net)))
  dtinfo=merge(dtnames,propag, by.x="names", by.y="name", sort=F)

  dtinfo[,2]=as.double(as.character(dtinfo[,2]))
  miniprop=dtinfo[order(dtinfo$value),]

  my_resolution = 5000
  my_palette    = colorRampPalette(c('white','orange1'))
  my_colors = my_palette(my_resolution)[as.numeric(cut(miniprop[,2], breaks=my_resolution))]

  miniprop[,3]=my_colors
  colnames(miniprop)<-c("Name","Propagation","Color")

  color_info=merge(dtnames,miniprop, by.x="names", by.y="Name", sort=F)

  color_info=merge(color_info, totaldataframe, by.x="names", by.y="name", sort=F )

  colnames(color_info)<-c("Name","Propagation","Color","Starting expression")

  return (color_info)

}

#' Assign colours to node types in the igraph object
#'
#' Support function for the network visualization
#'
#' @param subnet igraph object
#' @import igraph
#' @import grDevices
#'
fun2 <- function(subnet) {
  dtinfo=get.data.frame(subnet, what="vertices")
  miniprop=dtinfo[order(dtinfo$prop),]
  my_resolution = 5000
  my_palette = colorRampPalette(c('white','orange1'))
  my_colors = my_palette(my_resolution)[as.numeric(cut(miniprop[,4], breaks=my_resolution))]
  miniprop[,6]<-my_colors

  dtnames<-data.frame(names = names(V(subnet)))
  colordata=merge(dtnames, miniprop, by.x="names", by.y="name", sort=F)


  V(subnet)$color_rel<-colordata[,6]

  return(subnet)
}


#' Starts a shiny gui for the visualization of propagation results
#'
#' Given the graph obtained with create_net2plot function, the annotation of the genes included in the graph with
#' download_genes_annotation and the length of the chromosomes, it starts a graphical user interface to visualize
#' the propagation values on the network including genes and fragments.
#' This is useful to visualize the results of the two step propagation or the single gene propagation.
#'
#' @param net1 igraph object to visualize
#' @param ann_net_b dataframe, specifies the genes positions (chromosome, start and end position)
#' @param chr_len dataframe, specifies the chromosomes (first column) and their length  (second column)
#' @param example boolean, show 5 genes to select as an example
#' @return NULL
#' @import shiny
#' @import shinydashboard
#' @import shinyjs
#' @import RCy3
#' @import visNetwork
#' @import plotfunctions
#' @import fields
#' @import DT
#' @import spam
#' @import dplyr
#' @import igraph
#' @import plyr
#' @import grDevices
#' @import graphics
#' @import stats
#' @export
start_GUI = function(net1, ann_net_b, chr_len, example=FALSE){

  #index annotated genes to chromosomes
  rownames(chr_len)=chr_len$chrom
  index_all=list()
  for (i in c(1:19,"X","Y")){
    index_all[[i]]=which(ann_net_b[,"chr"] %in% paste("chr",i,sep=""))
  }
  names(index_all)=paste("chr",names(index_all),sep="")

  if(example){
    for(i in 1:length(index_all)){
      index_all[[i]]=seq(1,min(dim(ann_net_b)[1],5))
    }
  }

  #Get connected components
  dg<-decompose.graph(net1)
  length_dg=1:length(dg)
  name_in_dg=list()
  for (i in 1:length(dg)){
    length_dg[i]=length(V(dg[[i]]))
    name_in_dg[[i]]=igraph::get.vertex.attribute(dg[[i]],"name")
  }
  if(sum(igraph::get.vertex.attribute(net1)$shape == "square") > 0){
    wp <- wellPanel(
      h3("Legend",style="color:#0066CC"),
      fluidRow(
        column(12,h5(tags$b("Node's shape: "))),
        column(6,h5(tags$b("Genes")),icon("fas fa-circle")),
        column(6,h5(tags$b("Fragments")),icon("fas fa-square"),style="height:100px")
      ),
      fluidRow(
        column(12,h5(tags$b("Edge's color: "))),
        column(6,h5(tags$b("Gene-fragment or gene-gene connection: orange")), style="color:#E96928"),
        column(6,h5(tags$b("Fragment-fragment connection: purple")),style="color:#B616DE")
      ),
      br(),
      dataTableOutput("legend")
    )
  }
  else {
    wp <- wellPanel(
      h3("Legend",style="color:#0066CC"),
      fluidRow(
        column(12,h5(tags$b("Node's shape: "))),
        column(6,h5(tags$b("Genes")),icon("fas fa-circle")),
      ),
      fluidRow(
        column(12,h5(tags$b("Edge's color: "))),
        column(6,h5(tags$b("Gene-gene connection: orange")), style="color:#E96928"),
      ),
      br(),
      dataTableOutput("legend")
    )
  }

  ui <- dashboardPage(
    dashboardHeader(
      title=span("Visualization tool"),
      tags$li(class = "dropdown",
              tags$style(".main-header {max-height: 80px}"),
              tags$style(".main-header .logo {height: 60px;}"))
    ),

    dashboardSidebar(
      fluidRow(
        column(
          width = 12,
          selectInput("chr", shiny::HTML("<p><span style='color: blue'>Select nodes by chromosome</span></p>"),choices =names(index_all), multiple = FALSE),
          sliderInput( "region", shiny::HTML("<p><span style='color: blue'>Select nodes by genome region</span></p>"),min=1,max=195456987,step=10000,value=c(1,195456987)),
          selectizeInput(
            'node', shiny::HTML("<p><span style='color: blue'>Select or type node</span></p>"), choices = c("node1","node2"),
            options = list(maxOptions = 2,placeholder = 'Please select an option below',
                           onInitialize = I('function() { this.setValue(""); }'))
          ),
          numericInput("dist", shiny::HTML("<p><span style='color: blue'>Distance from selected node</span></p>"), 0, min = 0, max = 100),
          br(),
          actionButton("button","Search", style="background-color:#0066CC; color:#FFF; border-color: #004C99")
        )
      )

    ),
    dashboardBody(
      tags$head(tags$style(HTML('.content-wrapper, .right-side {
      background-color: #E2ECEB}
      .skin-blue .main-sidebar {
      background-color: #C2CEFF;
      }
'))),


      useShinyjs(),
      fluidRow(class="myRow",
               column(
                 width = 9,
                 visNetworkOutput("network")#,
                 #div(style = "height:1000px")
               ),
               column(width=3, fluid=FALSE,
                      plotOutput("legendG")#,
                      #div(style = "height:1000px")
               )
      ),

      checkboxInput("relativebox","Scale colours based by only this subnetwork", value = FALSE, width = NULL),
      fluidRow(
        column(1,actionButton("morebutton","More", style="background-color:#0066CC; color:#FFF; border-color: #004C99", offset=0)),
        column(1,actionButton("lessbutton","Less", style="background-color:#0066CC; color:#FFF; border-color: #004C99")),
        column(3,actionButton("cys","Open in cytoscape", style="background-color:#0066CC; color:#FFF; border-color: #004C99"),offset = 5),
        column(2,downloadButton("export",style="background-color:#0066CC; color:#FFF; border-color: #004C99"))
      ),
      br(),
      #legend space
      wp
    )
  )


  server <- function(input, output, session) {

    value<-reactiveVal(NULL)

    observeEvent(input$chr,{
      cat("observeEvent chr \n")
      value_max=chr_len[input$chr,]$length
      element<<-ann_net_b[index_all[[input$chr]],]
      updateSliderInput(session,"region",max=value_max,value = c(1,value_max),min = 1)

    })

    observeEvent(input$region,{
      cat("observeEvent region \n")
      element_pos<-ifelse((element[,3]>=input$region[1] & element[,3]<=input$region[2]) |
                            (element[,4]>=input$region[1] & element[,4]<=input$region[2]),TRUE,FALSE)
      element_select<-element[element_pos,1]

      updateSelectizeInput(session, 'node', choices = element_select, server = TRUE,options = list(maxOptions = length(element_select)))

    })

    #find the connected component to which the node belongs to...
    observeEvent(input$node, {
      cat("observeEvent node \n")
      if (input$node == ""){
        return()
      }
      print(input$node)
      for (i in 1:length(dg)) {
        if (input$node %in% name_in_dg[[i]]){
          comp<-dg[[i]]
          break
        }
      }
      print(comp)
      #...and set the max distance based on the result
      max=max(distances(comp,v=input$node, V(comp), weights=NULL))
      cat("node max:");print(max);cat("\n")
      updateNumericInput(session, "dist",value = 0, max = max)

    })


    #reaction to search button
    observeEvent(input$button, {
      cat("observeEvent search button \n")
      value_dist<<-input$dist
      #generate the graph
      subnet=make_ego_graph(net1, order = input$dist, nodes = input$node, mode = "all", mindist = 0)
      plot_gr<<-subnet[[1]]
      plot_gr<<-fun2(plot_gr)
      #vertex_attr(plot_gr,"color",which(V(plot_gr)$name==input$node))<-"#00FF00"
      dati<-toVisNetworkData(plot_gr,idToLabel = TRUE)
      newValue=dati
      value(newValue)
    })

    #reaction to more button
    observeEvent(input$morebutton, {
      cat("observeEvent more button \n")
      value_dist<<-value_dist+1
      cat("value_dist:",value_dist,"\n")
      cat("input node:");print(input$node);cat("\n")
      #generate the graph
      subnet=make_ego_graph(net1, order = value_dist, nodes = input$node, mode = "all", mindist = 0)
      plot_gr<<-subnet[[1]]
      plot_gr<<-fun2(plot_gr)
      dati<-toVisNetworkData(plot_gr,idToLabel = TRUE)
      newValue=dati
      value(newValue)
    })

    #reaction to less button
    observeEvent(input$lessbutton, {
      cat("observeEvent less button \n")
      value_dist<<-value_dist-1
      if(value_dist<0){
        value_dist<<-0
      }

      #generate the graph
      subnet=make_ego_graph(net1, order = value_dist, nodes = input$node, mode = "all", mindist = 0)
      plot_gr<<-subnet[[1]]
      plot_gr<<-fun2(plot_gr)
      dati<-toVisNetworkData(plot_gr,idToLabel = TRUE)
      newValue=dati
      value(newValue)
    })

    #visualization of the visNetowrk
    output$network <- renderVisNetwork({
      cat("observeEvent renderVisNetwork \n")
      if(is.null(value())){
        return()
      }

      if (input$relativebox==TRUE) {
        data=value()
        data$nodes$color=value()$nodes$color_rel
        data$nodes$color_rel=value()$nodes$color
      } else {
        data=value()
      }

      data$nodes$color=replace(data$nodes$color,which(data$nodes$id==input$node),"#00FF00")
      prv=data$nodes$propagation
      qs=quantile(prv,probs = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
      for(k in 1:length(qs)){
        if(k==1){
          qs_names=paste(qs[k],sep="-")
        }else{
          new_name=paste(qs[k-1],qs[k],sep="-")
          qs_names=c(qs_names,new_name)
        }
      }
      map=data.frame(indx=seq(1,11),qs_names,stringsAsFactors = F)
      propagation_ranges=findInterval(prv, qs)
      propagation_ranges=plyr::mapvalues(propagation_ranges,from=map$indx,to=map$qs_names,warn_missing = F)
      data$nodes$propagation_ranges=propagation_ranges

      #creation of the visNetwork object
      save(data,file="data.rda")
      vis<<-visNetwork(nodes = data$nodes, edges =data$edges, width = "auto", height = "auto")%>%
        visNodes( font = list("color"="#000",size=15),size=20,) %>%
        visEdges(length=150, width = 2, smooth = T) %>%
        visLayout(randomSeed = 250) %>%
        visEvents(stabilizationIterationsDone="function () {this.setOptions( { physics: false } );}") %>%
        visPhysics(maxVelocity = 50, minVelocity = 1) %>%
        visOptions(highlightNearest = list(enabled = T, hover = T),
                   selectedBy = list(variable = "propagation_ranges" , multiple=T, sort=T, main="Select by propagation ranges",
                                     style = 'width: 220px; height: 26px;
                                   background: #f8f8f8;
                                   color: darkblue;
                                   border:none;
                                   outline:none;'))


    })

    #visualization of the reactive sidebar legend
    output$legendG=renderPlot({
      cat("observeEvent renderPlot \n")
      par(bg="#E2ECEB")
      emptyPlot(0,0,axes=F)
      if(is.null(value())){
        return()
      }

      if(input$relativebox==TRUE) {
        if (length(value()$nodes$propagation)==1) {
          legend("bottomright", bty = "n", fill=value()$nodes$color_rel,legend=round(value()$nodes$propagation, digits=2), border = 'black',cex = 1.2)
        } else {
          zlim=c(round(min(value()$nodes$propagation),digits=2),round(max(value()$nodes$propagation), digits=2))
          my_palette    = colorRampPalette(c('white','orange1'))
          image.plot(legend.only=T, horizontal=F, legend.shrink=0.5, legend.width=0.65, col= my_palette(5000), zlim=zlim, axes=F,
                     axis.args=list(at=zlim, labels=zlim),smallplot=c(0.1,0.2,0.1,0.9))
        }

      } else {
        zlim=c(round(min(V(net1)$propagation),digits=2),round(max(V(net1)$propagation),digits=2))
        my_palette    = colorRampPalette(c('white','orange1'))
        image.plot(legend.only=T, horizontal=F, legend.shrink=0.5, legend.width=0.65, col= my_palette(5000), zlim=zlim, axes=F,
                   axis.args=list(at=zlim, labels=zlim),smallplot=c(0.1,0.2,0.1,0.9))

      }

    })

    #button to open the graph in cytoscape
    observeEvent(input$cys,{
      cat("observeEvent cys \n")
      cytoscapePing()
      createNetworkFromIgraph(plot_gr,"myIgraph")
    })


    #download network
    output$export <- downloadHandler(
      filename = "Network.html",
      content = function(con) {
        visSave(vis, con)
      })

    #output of the complete legend
    output$legend <- renderDataTable({
      if(input$button) {
        data_leg<-value()$nodes
        cols=c("label","shape","expr","propagation","color","color_rel","title","id")
        datatable(data_leg[,cols],options=list(columnDefs=list(list(visible=F,targets=c(1,2,7,8)))),colnames =c("Starting expression"="expr", "Propagation"="propagation", "Color based on the global network"="color", "Color based on the subnetwork"="color_rel")) %>%
          formatStyle (columns = "Color based on the global network", backgroundColor=styleEqual(data_leg$color ,as.character(data_leg$color)), color=styleEqual(data_leg$color ,as.character(data_leg$color))) %>%
          formatStyle (columns = "Color based on the subnetwork", backgroundColor=styleEqual(data_leg$color_rel ,as.character(data_leg$color_rel)), color=styleEqual(data_leg$color_rel ,as.character(data_leg$color_rel)))

      }

    })

  }

  shiny::shinyApp(ui, server)

}

#' Connects network nodes to metainformation
#'
#' Support function for the network visualization
#'
#' @param net1 igraph object
#' @param ann_net_b annotation dataset of the nodes included in the igraph object
#' @param frag_pattern character string to identify the fragment nodes by their name in the igraph object
#' @param ff_net edge list of the gene fragment network
#' @import igraph
#'
toolkitf <- function(net1, ann_net_b, frag_pattern = "F", ff_net = NULL) {
  toolnames=V(net1)$name
  if(is.null(ff_net)){
    symbol_position = 1:length(V(net1))
  }else {
    symbol_position = grep(frag_pattern,toolnames, invert=T)
  }
  symbol_names=toolnames[symbol_position]

  if(is.null(ff_net)){
    frag_position = 0
  }else {
    frag_position = grep(frag_pattern,toolnames, invert = F)
  }
  frag_names=data.frame(names=toolnames[frag_position])

  data_annotation=merge(frag_names, ann_net_b, by.x="names", by.y="ID", sort=F)

  #correct indices for annotation
  ann_frag_position=which(toolnames %in% data_annotation[,1])
  frag_link=paste("<a target='_blank' href=https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
    data_annotation$chr,"%3A",data_annotation$start,"%2D",data_annotation$end,
    "&hgsid=981993861_VmexClFuvk6xFB6uwwbxkArgAkC9> Annotation of ",
    data_annotation$names,"</a>",sep="")
  vector_link=paste("<a target='_blank' href=https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
    sep = "",symbol_names, "> GeneCards of ", symbol_names,"</a>")

  nolink_vector=c(rep("No link available", length(toolnames)))
  toolkit_title=replace(nolink_vector,ann_frag_position,frag_link)
  toolkit_title=replace(toolkit_title,symbol_position,vector_link)

  V(net1)$title=toolkit_title

  return(net1)
}
