library(shinyjs)
library(dplyr)
library(plyr)
library(stringr)
library(piano)
library(shinycssloaders)
require(visNetwork)
library(qvalue)
library(DT)
library(shinyWidgets)
library(ggwordcloud)
library(shinythemes)
library(igraph)
library(shinyjs)


#read data

TF.tbl=read.table("../data/TF_DS.full.tbl", header=T, sep="\t")
TF.tbl=TF.tbl[,c("TF","Rank","Drought_score","Description","Family","Label")]
colnames(TF.tbl)=gsub("_"," ",colnames(TF.tbl))


#MSU to RAP mappings
rap.msu=read.table("../data/msu_rap_mapings", header=T)

#consensus net
con_net=read.table("../data/consensus.network", header=T)

#coreg net
net=read.table("../data/consensus.gene-TF_overlap.JI_mt0pt1.top1perc.dat", header=T)
colnames(net)=c("src","target","ji")

#hetnet
hetnet=read.table("../data/hetnet.dat", header=T)
all.nodes=read.table("../data/hetnet.nodes", header=T)
all.nodes.titles=read.table("../data/hetnet.nodes.titles")
colnames(all.nodes.titles)=c("Node","title")
all.nodes=all.nodes %>% left_join(all.nodes.titles, by = c("Node" = "Node"))
Node1=as.character(unique(hetnet$from))
Node2=as.character(unique(hetnet$to))
nodes = data.frame(name = unique(c(Node1, Node2)), stringsAsFactors = FALSE)
indx=match(nodes$name, all.nodes$Node)
hetnet.nodes=all.nodes[indx,]
hetnet.nodes$color=gsub("olivegreen","#009966",hetnet.nodes$color)
hetnet.nodes$color=gsub("lightblue","#0066FF",hetnet.nodes$color)
hetnet.igraph <- graph_from_data_frame(hetnet, directed=FALSE, vertices=hetnet.nodes)


#module net
modnet=read.table("../data/TF-module.JI_mt0.1.dat")
colnames(modnet)=c("TF","Module","Jaccard")

#modules
modules=read.table("../data/module_attribute.tbl", header=T)
module.gs=modules[,c("Gene","Module")]
module.gs=module.gs[!grepl("M0001",module.gs$Module),] #remove M0001
mod_to_use=as.data.frame(table(unlist(module.gs$Module)))
module.gsets=mod_to_use %>% filter(Freq > 10) %>% filter(Freq < 500)
module.gs=module.gs[module.gs$Module %in% module.gsets$Var1,]
module.gs=loadGSC(module.gs)


#GO enrich mat
GO.mat=read.table("../data/GO.mcl.mat", header=T)
colnames(GO.mat)=c("Module","GO")
GO.mat$GO=gsub("_"," ",GO.mat$GO)
GO.mat$GO=str_to_lower(GO.mat$GO)
GO.mat=plyr::ddply(GO.mat, .(Module), colwise(paste), collapse = ",  ")

#mapman
mapman.mat=read.table("../data/mapman.mcl.mat", header=T)
colnames(mapman.mat)=c("Module","mapman")
mapman.mat$mapman=gsub("_"," ",mapman.mat$mapman)
mapman.mat$mapman=str_to_lower(mapman.mat$mapman)
mapman.mat=plyr::ddply(mapman.mat, .(Module), colwise(paste), collapse = ",  ")

#fire
fire.mat=read.table("../data/fire.mcl.mat", header=T)
colnames(fire.mat)=c("Module","fire")
fire.mat$fire=gsub("_"," ",fire.mat$fire)
fire.mat$fire=str_to_lower(fire.mat$fire)
fire.mat=plyr::ddply(fire.mat, .(Module), colwise(paste), collapse = ",  ")


#CREs
cre=read.table("../data/module-motif.mat", header=T)

#Example genelist 
gene_annot=read.table("../data/msu_gene_annotations.txt")
colnames(gene_annot)=c("Gene","Annotation")
gene_annot$Annotation=gsub("(.*),.*", "\\1", gene_annot$Annotation)

GOBP=read.table("../data/gene_association.rice.agbase.propagated.BP.for_overlap_calculation.lt1000mt10.JI0pt9.setdiff5.filtered.genesets")
GOBP.desc=read.table("../data/gene_association.rice.agbase.propagated.BP.for_overlap_calculation.lt1000mt10.JI0pt9.setdiff5.filtered.genesets.desc")
colnames(GOBP.desc)=c("GOTerm","Description")
GOBP.gs=loadGSC(GOBP)


