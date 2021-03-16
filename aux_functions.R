

parse.textarea.input <- function(text){
  sep <- NULL
  if(grepl(";",text)) sep <- ";"
  if(grepl(",",text)) sep <- ","
  if(grepl("\n",text)) sep <- "\n"
  if(is.null(sep)) {
    text <- text
  } else {
    text <- unlist(stringr::str_split(text,sep))
  }
  text=text[text != ""]
  return (text)
}
    
    
# function to check if input locus IDs are MSU format
	not_msu <- function(input){
  		check=startsWith(as.character(input), "LOC")
  		#print(check)
		if ('FALSE' %in% check)
		{
			showModal(modalDialog(
        		title = "Gene not found!",style = "color: red;",
        		"The uploaded file has incorrect locus IDs, Please check for MSU format. Syngenta chromosomes are not present in the network.",
        		easyClose = TRUE,
        		footer = NULL
      		))
 		   	showNotification("Stopping progress...", duration = 5,type="error")
 		  	return("Error! File contains incorrect locus IDs...MSU format only")
  		}
  		else
  		{
	    	NULL 
	    }
	}
	
	 TextareaData <- eventReactive(input$list,{ 
    	showNotification("Parsing input...",type="message")
    	 validate(
     	 	need(input$inputGeneList, 'Enter at least one gene!')
    	)
    	inputGeneList=input$inputGeneList
    	if(is.null(inputGeneList))
    	{
        	showModal(modalDialog(
        		title = "Empty input box!",style = "color: red;",
        		"Enter at least one gene in MSU format.",
        		easyClose = TRUE,
        		footer = NULL
      		))
        }
    	if(!is.null(inputGeneList))
    	{
        	geneList <- parse.textarea.input(inputGeneList)
        	validate(not_msu(geneList))
        	geneList
        }
      })
	

	# function to prep enrichment table 
	get_enrichment_results <- function(input){
	usrgenes=input
	showNotification("Starting enrichment analysis", duration = 5)
	enriched.gs.pval=enrich_hyper(usrgenes)
	enriched.gs.pval=enriched.gs.pval[enriched.gs.pval$adjusted_pvalue<0.1,]
	if (nrow(enriched.gs.pval)>0)
		{	
			enrich.tbl=join_enrich_tbl(enriched.gs.pval)
			enrich.tbl
		}
		else
		{
			showNotification("Warning: no enrichment found with input genes. Expanding the input using first degree neighbors", duration = 5, type="warning")
			f=net[net$target %in% usrgenes | net$src %in% usrgenes,] 
			firstNei=unique(as.character(f$src,f$target))
			enriched.gs.pval=enrich_hyper(firstNei)		
			enriched.gs.pval=enriched.gs.pval[enriched.gs.pval$adjusted_pvalue<0.1,]
			if (nrow(enriched.gs.pval)>0)
			{
				enrich.tbl=join_enrich_tbl(enriched.gs.pval)
				enrich.tbl
			}
			else
			{
				showNotification("No significant overlap found.", duration = 5,type="error")
				enriched.gs.pval=data.frame("Gene"="","Module"="No enrichment found! Try another geneset")
				enriched.gs.pval
			}
			
		}
	}
	
	get_regulon_enrichment_results <- function(input){
	usrgenes=input
	enriched.gs.pval=enrich_hyper_regulon(usrgenes)
	enriched.gs.pval=enriched.gs.pval[enriched.gs.pval$adjusted_pvalue<0.1,]
	if (nrow(enriched.gs.pval)>0)
		{	
			enrich.tbl=join_Regulon_enrich_tbl(enriched.gs.pval)
			enrich.tbl
		}
		else
		{
			
				showNotification("No significant overlap found.", duration = 5,type="error")
				enriched.gs.pval=data.frame("Gene"="","GO Term"="No enrichment found! Try another geneset")
				enriched.gs.pval
			
			
		}
	}
	
	
	
  # function to perform enrichment analysis using Hypergeo
	enrich_hyper <- function(input,output,session){
		#showNotification("Enrichment network ...", duration = 5,type="message")
		hyper_query=unique(input)
		enriched.gs=runGSAhyper(hyper_query,gsc=module.gs)
		resTab=as.data.frame(enriched.gs$resTab)
		resTab=resTab[,1:4]
		colnames(resTab)=c("pvalue","adjusted_pvalue","Overlap","Module_Size")
		resTab$Module=rownames(resTab)
		resTab=resTab[,c("Module","Module_Size","Overlap","adjusted_pvalue","pvalue")]
		resTab=resTab[,1:4]
		resTab$adjusted_pvalue=round(resTab$adjusted_pvalue,5)
		return(resTab)
	}
      
	enrich_hyper_regulon <- function(input,output,session){
		hyper_query=unique(input)
		enriched.gs=runGSAhyper(hyper_query,gsc=GOBP.gs)
		resTab=as.data.frame(enriched.gs$resTab)
		resTab=resTab[,1:2]
		colnames(resTab)=c("pvalue","adjusted_pvalue")
		resTab$GOTerm=rownames(resTab)
		resTab=resTab[,c("GOTerm","pvalue","adjusted_pvalue")]
		resTab$adjusted_pvalue=round(resTab$adjusted_pvalue,5)
		resTab$pvalue=round(resTab$pvalue,5)
		
		return(resTab)
	}
	
	# function to add attribute columns to enriched module table
	join_enrich_tbl <- function(input, output, session){
		tbl=input
		#merged.tbl=merge(tbl, GO.mat)
		#merged.tbl=merge(merged.tbl,mapman.mat)
		#merged.tbl=merge(merged.tbl,fire.mat)
		merged.tbl=tbl %>% left_join(GO.mat, by=c("Module"))
		merged.tbl=merged.tbl %>% left_join(mapman.mat, by=c("Module"))
		merged.tbl=merged.tbl %>% left_join(fire.mat, by=c("Module"))
		mod.list = list()
  		for (i in 1:length(unique(merged.tbl$Module)))
  		{
  			mod.regulators=modnet[modnet$Module %in% merged.tbl$Module[i],]
  			mod.regulators=mod.regulators[order(-mod.regulators$Jaccard),]
  			mod.regulators=mod.regulators[1:10,] #select top 10 regulators
			mod.regulators$TF=paste0("<a href='http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?&orf=",mod.regulators$TF,"' target='_blank'>",mod.regulators$TF,"</a>")
			mod.regulators=plyr::ddply(mod.regulators[,1:2], .(Module), colwise(paste), collapse = ",  ")
  			mod.list[[i]]= mod.regulators	
  		}
  		mod.TF = do.call(rbind, mod.list)
  		merged.tbl=merged.tbl %>% left_join(mod.TF, by = c("Module" = "Module"))
		merged.tbl=merged.tbl[order(merged.tbl$adjusted_pvalue),]
		colnames(merged.tbl)=c("Module ID","Size","Overlap","adjusted p-value","GO process","Mapman pathway","CRE","Top Regulators")
		return(merged.tbl)		
	}
	
	join_Regulon_enrich_tbl <- function(input, output, session){
		tbl=input
		merged.tbl=tbl %>% left_join(GOBP.desc, by=c("GOTerm"))
		#print(head(merged.tbl))
		merged.tbl=merged.tbl[order(merged.tbl$adjusted_pvalue),]
		colnames(merged.tbl)=c("GO Term","p-value","adjusted p-value","Description")
		merged.tbl$Description=gsub("_"," ",merged.tbl$Description)
		return(merged.tbl)		
	}
	
	#function to get neighbors of enriched modules
	get_module_neighbors <- function(input,output,session){	
		queryMods =input
		#queryMods =c("M0025")
		net1=hetnet[hetnet$from %in% queryMods | hetnet$to %in% queryMods,]
		Node1=as.character(unique(net1$from))
		Node2=as.character(unique(net1$to))
		nodes = data.frame(name = unique(c(Node1, Node2)), stringsAsFactors = FALSE)
		net2.igraph=induced_subgraph(hetnet.igraph,nodes$name)
		net2.visnet=toVisNetworkData(net2.igraph)
		net2.visnet$nodes$title=gsub("_"," ",net2.visnet$nodes$title)
		net2.visnet$nodes$label=gsub("_"," ",net2.visnet$nodes$label)
		net2.visnet$nodes$shadow=TRUE
		net2.visnet$nodes$borderWidth=3

		visnet=visNetwork(nodes = net2.visnet$nodes, edges = net2.visnet$edges) %>% 
		visGroups(groupname = "module", shape = "circle", color="red") %>%
  		visGroups(groupname = "Regulator", shape = "triangle",color="grey") %>%
		visGroups(groupname = "mapman", shape = "square",color="#0066FF") %>%
		visGroups(groupname = "CRE", shape = "diamond",color="#009966") %>%
		visGroups(groupname = "GOBP", shape = "square", color="cornflowerblue") %>%
		visInteraction(navigationButtons = TRUE) %>% 
	#	visPhysics(stabilization = TRUE) %>%
		visEdges(smooth = TRUE) %>%
  		visIgraphLayout() %>% visLegend(width = 0.1, position = "right", main = "Group",zoom=FALSE)
		return(visnet)	
	}	
	
	
	