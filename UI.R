#source('global.R')
#Shiny 
UI=shinyUI(
	fluidPage(
		shinyjs::useShinyjs(),
		id="main",
		title="GRAiN",
		theme = shinytheme("simplex"),
		tags$head(tags$style(
                  HTML('
                       body, label, input, button, select { 
                       font-family: "Helvetica","Arial";
                       
                       } 
                       .btn {
                        border-color: grey; 
                        background: grey;
                       }
                       
                      .bttn-bordered.bttn-sm {
                          width: 200px;
                          text-align: left;
                          margin-bottom : 0px;
                          margin-top : 20px;
                       }
                       #refresh{
                       	
                       	margin-top: 40px;
                       }
                       .div{
                       	margin-top : 20px;
                       }
    					table.dataTable th {
    						background-color: #555555 !important;
    						color: white;
    					}
    					#netColumn{
    						border: 1.5px solid #F0F0F0;
    						border-top:0px;
    						margin-left:5px;
    					}
                       '
                  )
                )),
		titlePanel("GRAiN"),
		
		tabsetPanel(
			tabPanel("Home",  
                icon = icon("object-group"),
                fluidRow(
                	column(2,
               			dropdown(id="sidebar",label="Gene list upload",
               			icon = icon("list"),
                    		style = "bordered", 
                    		status = "primary", 
                    		width = "200px",
                    		size =  "sm",
                    		animate = animateOptions(
                    			enter = animations$fading_entrances$fadeInDown,
                   	     		exit = animations$fading_exits$fadeOutUp
                    		),                    	
               				textAreaInput(inputId = "inputGeneList",
                        		label = "Gene list: separate gene names by , or ; or newline",
                            	value =  NULL, 
                            	width = "100%",height="200px"
                        	),
                        	actionButton("list", "Submit",icon=icon("upload")), 
                        	actionButton("demo", "Demo"), hr(),
                        ),
                        
               			dropdown(id="TargetSelect",label="Search a TF",
                        icon = icon("search"),
                    		style = "bordered", 
                    		status = "primary", 
                    		width = "300px",
                    		size =  "sm",
                    		animate = animateOptions(
                    			enter = animations$fading_entrances$fadeInDown,
                   	     		exit = animations$fading_exits$fadeOutUp
                    		),         	
                        	textInput("TF", "Enter a locus ID",value="LOC_Os03g53020"),
                        	actionButton("SearchRegulon", "Submit",icon=icon("search")),
      						hr()
               			),#end dropdown
               			actionButton("refresh", "Refresh",icon=icon("refresh")),
                		hr(),	           
                		hidden(
                			div(id="logs",
                				textOutput("Log"),
                				                		
 							)	
                		)
                	),#end column
                	column(9,id="netColumn",
                		div(id="landingpagetbl",
                			DT::dataTableOutput("TFDSmaintbl")
                		),
                		hidden(
                			div(id="networks",
                				visNetworkOutput("netplot", width = "100%", height = "250px")%>% withSpinner(color="black", type=6, size=1)
                			),
                			div(id="demonetworks",
                 				visNetworkOutput("demonetplot", width = "100%", height = "250px")%>% withSpinner(color="black", type=6, size=1)               	
                			),
                			div(id="SingleRegulon",
                				sliderInput("nedges","# of edges to display",min = 5,max = 100,value = 10),
                 				visNetworkOutput("RegulonPlot", width = "100%", height = "250px")%>% withSpinner(color="black", type=6, size=1)               	
                			)
                		)
                	),
                	column(1)
                	
                ),#fluidRow
                fluidRow(	
                	div(id="table",hr(),
                		hidden(
                		column(2,id="download-btns",
                			
                				dropdown(label="Download reports",
                					icon = icon("download"),
                    				style = "bordered", 
                    				status = "primary", 
                    				width = "200px",
                    				size =  "sm",
                    				animate = animateOptions(
                    					enter = animations$fading_entrances$fadeInDown,
                   	     				exit = animations$fading_exits$fadeOutUp
                    				),             
                					downloadButton("downloadTblbtn","Download this table")
                				)
                			)
                		),
                		hidden(
                		column(2,id="download-Regulon-btns",
                			
                				dropdown(label="Download reports",
                					icon = icon("download"),
                    				style = "bordered", 
                    				status = "primary", 
                    				width = "300px",
                    				size =  "sm",
                    				animate = animateOptions(
                    					enter = animations$fading_entrances$fadeInDown,
                   	     				exit = animations$fading_exits$fadeOutUp
                    				),             
                					downloadButton("downloadRegulonbtn","Download target table"),hr(),
                					downloadButton("downloadRegulonEnrichmentbtn","Download enrichment table")
                					
                				)
                			)
                		),
                		column(9,
                			DT::dataTableOutput("gsenrichtbl"),	           
 							DT::dataTableOutput("demorestbl"),
 							DT::dataTableOutput("Regulontbl"),
 							hr(),

 							DT::dataTableOutput("RegulonEnrichtbl")
 						),
 						column(1)
 					)	
                )#fluidRow
			),#end tab panel
			
			tabPanel("About", icon = icon("info"),
				includeMarkdown("about.Rmd")
            ),
            tabPanel("Help", icon = icon("book"),
				includeMarkdown("docs.Rmd")
            )
		) #end of tabset panel   
	)#end of fluid page		
)#end of shinyUI
