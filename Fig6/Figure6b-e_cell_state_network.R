
####################################################Figure 6b-6e network graph#########################################


library(igraph)
library(plotrix)
library(s2dverification)
library(RColorBrewer)
library(dvmisc)


#This plot needs a dataset summarizing the total number of clonotypes shared by two cell state from at baseline to on induction.
###########define functions################
get_adjacency <- function(data,response){
  cells=unique(data$node.baseline)[c(2:8,10,11)]#only include TRM, TEM, effector and proliferating subsets
  #cells=unique(data$node.baseline)
  adjacency_matrix_raw=matrix(nrow = length(cells),ncol = length(cells))
  group=responses[response]
  
  for( i in 1:length(cells)){
    for(j in 1:length(cells)){
      ind=which(data$node.baseline==cells[i]&data$node.oninduction==cells[j])
      if(data[ind,group]<1){
        adjacency_matrix_raw[i,j]=0
      }
      else{
        adjacency_matrix_raw[i,j]=data[ind,group]
      }
    }
  }
  rownames(adjacency_matrix_raw)=cells
  colnames(adjacency_matrix_raw)=cells
  adjacency_matrix_raw=adjacency_matrix_raw/normalization_factor[response]
  return(adjacency_matrix_raw)
}
plotNetwork<-function(adjacency_matrix_raw,response,weight_all){
  adjacency_matrix=adjacency_matrix_raw
  adjacency_matrix[(adjacency_matrix>0)]<-1
  
  graph_object <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
  graph_object$layout=layout.circle
  V(graph_object)$color=colors[match(names(V(graph_object)),unique(dat$node.1))]
  E(graph_object)$color='grey'
  weight=as.vector(t(adjacency_matrix_raw))[which(as.vector(t(adjacency_matrix_raw))!=0)]
  E(graph_object)$width=weight*30
  E(graph_object)$color=colorRampPalette(c('grey90','grey75','moccasin','gold','orange','darkorange3','firebrick'))(length(weight_all))[match(weight,weight_all)]
  p<-plot(graph_object, 
          vertex.size = 30,vertex.label.color='black',
          edge.arrow.size = 0.25, edge.arrow.width = 2,main=paste0(responses[response]),
          edge.curved = TRUE)
  return(p)
  
}
colors=c("#A6CEE3", "#4F96C4","royalblue" ,"#CC3311", "#EE3377", "#FB9A99", "#A7D78D","#69BB54","#5D9E43" ,"#117733" , "#DDCC77", "#D9A398", "#CAB2D6" ,"#9A77B8", "#6A3D9A"  )   


###########read data#################
#if plotting expaned clonotypes
dat=read_xlsx('../DataShare/No of expanded clonotypes shared between two states at two timepoint.xlsx')
#or if plotting contracted clonotypes
dat=read_xlsx('../No of contracted clonotypes shared between two states at two timepoint.xlsx')

responses=c('MHR','nMHR')

normalization_factor=c(47,31)#if plotting expaned clonotypes
normalization_factor=c(44,13)#if plotting contracted clonotypes


adjacency_matrix_raw.MHR=get_adjacency(dat.1,1)
adjacency_matrix_raw.nonMHR=get_adjacency(dat.1,2)

weight.MHR=as.vector(t(adjacency_matrix_raw.MHR))[which(as.vector(t(adjacency_matrix_raw.MHR))!=0)]
weight.nonMHR=as.vector(t(adjacency_matrix_raw.nonMHR))[which(as.vector(t(adjacency_matrix_raw.nonMHR))!=0)]
weight_all=sort(unique(c(weight.nonMHR,weight.MHR)))

plotNetwork(adjacency_matrix_raw.MHR,1,weight_all)
plotNetwork(adjacency_matrix_raw.nonMHR,2,weight_all)

pdf('./colorbar.pdf',width = 1,height = 2)
ColorBar(round(c(0,weight_all),digits = 2),
         colorRampPalette(c('grey90','grey75','moccasin','gold','orange','darkorange3','firebrick'))(length(weight_all)))
dev.off()
