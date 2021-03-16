source('aux_functions.R')

server = function(input, output, session) { 	
	
    textLog <- reactiveVal("")
    
    observeEvent(input$demo, {
    	shinyjs::show(id = "table")
    	shinyjs::show(id = "demonetworks")
    	shinyjs::hide(id = "download-btns")
    	shinyjs::hide(id = "SingleRegulon")
    	shinyjs::hide(id = "download-Regulon-btns")
    	shinyjs::hide(id="intro")
    	shinyjs::hide(id ="selectPathnetworks")
    	shinyjs::show(id ="netColumn")
    	shinyjs::show(id = "download-btns")

    })
     
    observeEvent(input$SearchRegulon, {
    	shinyjs::show(id = "SingleRegulon")
    	shinyjs::show(id = "table")
		shinyjs::show(id = "download-Regulon-btns")
    	shinyjs::hide(id="networks")
		shinyjs::hide(id="demonetworks")
		shinyjs::hide(id="download-btns")
		shinyjs::hide(id="gsenrichtbl")
		shinyjs::hide(id="demorestbl")

    })
     
    observeEvent(input$list, {
   		shinyjs::show(id = "logs")
    	shinyjs::show(id = "networks")
    	shinyjs::hide(id="intro")
    	shinyjs::show(id ="netColumn")
    	shinyjs::show(id = "download-btns")
    	shinyjs::hide(id = "SingleRegulon")
    	shinyjs::hide(id = "download-Regulon-btns")
    	genes=TextareaData()
    	qurygenesTbl=gene_annot[gene_annot$Gene %in% genes,]
		message=paste("Total genes in input:",length(qurygenesTbl$Gene),sep=" ");
    	textLog(paste(textLog(), message))
    	
    	enrichtbl=enrichTblData()
		message=paste("Total enriched modules:",length(unique(enrichtbl[,1])))
		textLog(paste(textLog(), message, sep="\n"))
		
		enrichedBP=unlist(strsplit(enrichtbl[,5],","))
		enrichedBP=enrichedBP[!is.na(enrichedBP)]
		totalBP=length(unique(enrichedBP))
		message=paste("\nTotal enriched GO BP:",totalBP,sep=" ");
		textLog(paste(textLog(), message, sep="\n"))
		
		enrichedmm=unlist(strsplit(enrichtbl[,6],","))
		enrichedmm=enrichedmm[!is.na(enrichedmm)]
		totalmm=length(unique(enrichedmm))
		message=paste("\nTotal enriched Mapman:",totalmm,sep=" ");
		textLog(paste(textLog(), message, sep="\n"))
		
		enrichedcre=unlist(strsplit(enrichtbl[,7],","))
		enrichedcre=enrichedcre[!is.na(enrichedcre)]
		totalcre=length(unique(enrichedcre))
		message=paste("\nTotal enriched CREs:",totalcre,sep=" ");
		textLog(paste(textLog(), message, sep="\n"))
		
		enrichedreg=unlist(strsplit(enrichtbl[,8],","))
		enrichedreg=enrichedreg[!is.na(enrichedreg)]
		totalreg=length(unique(enrichedreg))
		message=paste("\nTotal unique regulators:",totalreg,sep=" ");
		textLog(paste(textLog(), message, sep="\n"))

    })
    
    observeEvent(input$refresh, {
    	session$reload()
    }) 
    
	output$Log = renderText({
		textLog()	
	})

    
     demoData <- eventReactive(input$demo,{ 
     showNotification("Starting demo...",type="message")
    	 data=read.table("../data/sucrose.txt")
    	 geneList=data$V1
    	 geneList
      })
    
    	enrichTblData = reactive({
		genes=TextareaData()
		usrgenes=unique(genes)
		enrichment_table=get_enrichment_results(usrgenes)
		enrichment_table
  	})
	
    enrichRegulonData = reactive({
		TFNet=TFdata()
		src=as.character(unique(TFNet$src))
		target=as.character(unique(TFNet$target))
		nodes = data.frame(name = unique(c(src, target)), stringsAsFactors = FALSE)
		usrgenes=unique(nodes)
		enrichment_table=get_regulon_enrichment_results(usrgenes)
		enrichment_table
  	})
  	
	output$RegulonEnrichtbl = renderDT(rownames=FALSE,escape = FALSE,class = "cell-border stripe",caption="Table 2: GO Enrichment of targets",{
		#showNotification("Calculating enrichment..",type="message",duration = 5)
		enrichRegulonData()
	})
    
    output$demorestbl = renderDT(rownames=FALSE,escape = FALSE, class = "cell-border stripe",{
		genes=demoData()
		usrgenes=genes
		enrichment_table=get_enrichment_results(usrgenes)
		enrichment_table
	})

	output$demonetplot = renderVisNetwork({
		genes=demoData()
		usrgenes=genes
		enrichment_table=get_enrichment_results(usrgenes)
		visnet=get_module_neighbors(enrichment_table[,1])   
		visOptions(visnet, highlightNearest = TRUE, selectedBy = "group", nodesIdSelection = TRUE) 
	})	
	
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
      
  	
    
	output$gsenrichtbl = renderDT(rownames=FALSE,escape = FALSE,class = "cell-border stripe",caption="Table 1: Clusters enriched in the query gene set",{
		#showNotification("Calculating enrichment..",type="message",duration = 5)
		enrichTblData()
	})
	
	output$netplot = renderVisNetwork({
		genes=TextareaData()
		usrgenes=genes
		showNotification(paste("Total unique genes matching MSU annotations",length(unique(usrgenes)),sep=": "), duration = 5)
		enrichment_table=enrichTblData()
		
		if (enrichment_table$Module %in% 'No enrichment found! Try another geneset')
		{
			showModal(modalDialog(
        		title = "No enrichment found! Try another geneset!",style = "color: red;",
        		"Maybe too few genes were input. Consider add more genes to input and try again.",
        		easyClose = TRUE,
        		footer = NULL
      		))
 		   	showNotification("Stopping progress...", duration = 5,type="error")
  		}
  		else
  		{
  			showNotification("Plotting network", duration = 5)
			visnet=get_module_neighbors(enrichment_table[,1])   
			visOptions(visnet, highlightNearest = TRUE, selectedBy = "group", nodesIdSelection = TRUE) 
		}
	})	
	
	

	output$downloadTblbtn=downloadHandler(
    	filename = function() {
      		paste("GRAiNS.modules.output","txt", sep = ".")
    	},
    	content = function(file) {
      		write.table(enrichTblData(), file, row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    	}
  	)

	##Regulon searches
		
	TFdata = eventReactive(input$SearchRegulon, {
		validate(need(input$TF != "", "Please enter a gene"))
		gene = input$TF
		req(gene)
		
		gene=trimws(gene, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
		check=startsWith(gene, "LOC_Os")
		#israp=startsWith(gene, "Os")
		if (check=='TRUE')
		{
 		  	TFnet=con_net[con_net$src %in% gene,]
  		}
  		else 
  		{
  			showModal(modalDialog(
        		title = "Gene not found!",style = "color: red;",
        		"Either the ID format is incorrect, or it is not present in our network. Try another gene with format LOC_OsXXgXXXXX",
        		easyClose = TRUE,
        		footer = NULL
      		))
  		}
  		TFnet
  	})
  	
  	output$RegulonPlot = renderVisNetwork({
		showNotification("Plotting network...", duration = 5,type="message")
		TFNet = TFdata()
		if(dim(TFNet)[1]==0)
			{
      		showModal(modalDialog(
        		title = "Network not found!",style = "color: red;",
        		"There are no genes predicted as targets for this regulator!!",
        		easyClose = TRUE,
        		footer = NULL
      		))
    		}
		TFNet=TFNet[1:input$nedges,]
		src=as.character(unique(TFNet$src))
		target=as.character(unique(TFNet$target))
		nodes = data.frame(name = unique(c(src, target)), stringsAsFactors = FALSE)
		nodes$id = 0:(nrow(nodes) - 1)
		nodes$group = ifelse(nodes$name %in% src, "TF", "targets")
		indx=match(nodes$name, gene_annot$Gene)
		gene_desc.tmp=gene_annot[indx,]
		nodes=cbind(nodes, gene_desc.tmp)
		nodes=nodes[,2:5]
		modules=modules[,c("Module","Gene")]
		colnames(modules)=c("Module","Gene")
		nodes=left_join(nodes,modules)
		edges = TFNet %>%
		left_join(nodes, by = c("src" = "Gene")) %>%
		select(-src) %>%
		left_join(nodes, by = c("target" = "Gene")) %>%
		select(-target) 
		vis.nodes=nodes
		vis.links=edges[,c("id.x","id.y")]
		colnames(vis.links)=c("from","to")
		vis.links$weight=20
		vis.nodes$shadow <- TRUE # Nodes will drop shadow
		vis.nodes$title  <- vis.nodes$Annotation # Text on click
		vis.nodes$title <- gsub("_"," ",vis.nodes$title)
		vis.nodes$label  <- vis.nodes$Gene # Node label		
		vis.nodes$size   <- 30 # Node size
		vis.nodes$borderWidth <- 2 # Node border width
		#vis.nodes$shape="dot"
		vis.nodes$shape=ifelse(nodes$group %in% "TF", "triangle", "dot")
		vis.nodes$color.background <- c("slategrey")
		vis.nodes$color.border <- "black"
		vis.nodes$color.highlight.background <- "orange"
		vis.nodes$color.highlight.border <- "darkred"
		visnet=visNetwork(vis.nodes, vis.links) %>% 
		visGroups(groupname = "targets", shape = "circle", color="red") %>%
  		visGroups(groupname = "TF", shape = "triangle",color="grey")
		visOptions(visnet)		
		 
	})

	output$Regulontbl= renderDT(rownames=FALSE,escape = FALSE,class = "cell-border stripe",caption="Table 1: Predicted targets of query TF",{
		TFNet = TFdata()
		#print(head(TFNet))
		TFNet=TFNet[,3:4]
		colnames(TFNet)=c("Target","Annotation")
		TFNet$Annotation=gsub("_"," ",TFNet$Annotation)
		TFNet
	})
  	
	output$downloadRegulonbtn=downloadHandler(
    	filename = function() {
      		paste("GRAiNS.regulon.output","txt", sep = ".")
    	},
    	content = function(file) {
      		write.table(TFdata(), file, row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    	}
  	)

	output$downloadRegulonEnrichmentbtn=downloadHandler(
    	filename = function() {
      		paste("GRAiNS.regulon.GOEnrichment.output","txt", sep = ".")
    	},
    	content = function(file) {
      		write.table(enrichRegulonData(), file, row.names = FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    	}
  	)
}