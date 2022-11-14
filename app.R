##
## CRISPR web app - app version 41 
## compatible with R - 4.2.1
##

library(tidyverse)
library(data.table)
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinyalert)
library(sjmisc)
library(DT)
library(reshape2)
library(plotly)
library(ggbeeswarm)

setwd("/opt/shiny-server/samples/gftapps/analysis/crispr_data_gftapp/")

##Define input:-
#1) read metadata and gene-level scores
parpi_genez <- fread("combined_gene_level_scores_annotated_260922.txt", colClasses=(`CISPLATIN_PEMETREXED-MEAN-LFC:pub::M1`="numeric"))
parpi_genez <- as.data.frame(parpi_genez)
selectInput_data <- select(parpi_genez,matches(c(":gft-C",":pub::","^VE.")))
parpi_metadata <- read_tsv("v11.metadata.txt")
parpi_metadata <- as.data.frame(parpi_metadata)

#2) gene annotations
annotation_columns <- c("^DDR$","Process","CGC","CEG")

#3) drugs
parp_inhibitors <- c("Olaparib","Talazoparib", "Niraparib", "Rucaparib", "Veliparib")
atr_inhibitors <- c("AZD6738","VE821","VE-821","VX970","VX-970","M4344","M6620","BAY1895344")
platinum_inhibitors <- c("CISPLATIN","CARBOPLATIN","OXALIPLATIN")

polq_inhibitors <- c("ART558","ART182","ART693")

#4) Isogenic screens with no drug arm
ve_only <- "nodrug"

#4) metadata columns to display by default - 15 columns
cols_metadata_default_display <- c("arm_ID","cell_line","cas9_type","modification","timepoint","drug","target","DNA_repair","drug_conc_nM","SF","library","library_type","analysis_pipeline","analysis_method","short_arm_key")

#5) thresholds
pos.hit.threshold <- 2
neg.hit.threshold <- -2

# shiny UI
ui <- navbarPage(title=div(a(href="http://gft.icr.ac.uk/",img(src="icr-3-1.jpg",height='40',width='65'),img(src="BCN_Logo_LargeStrap_Pink_RGB.PNG",height='35',width='70')),"CRISPR screens",titleWidth='100%', style="color: dimgrey; font-size: 28px; font-family: sans-serif;"),
                 windowTitle="GFT CRISPR app",
                 #Select datasets tab
                 tabPanel("Select datasets",
                          shinyjs::useShinyjs(debug=TRUE),
                          div(id="Sidebar",sidebarPanel(
                            #style = "position:fixed;width:27%;",
                            style = "width:88%;",
                            conditionalPanel(
                              'input.dataset == "parpi_metadata',
                              checkboxGroupInput("show_vars", "Select metadata columns to show:",
                                                 names(parpi_metadata), selected = cols_metadata_default_display) #14 metadata columns to display by default
                            ))),
                          wellPanel(
                            p("Note:- 'Select datasets' tab displays a table with all the information about CRISPR screen data. More columns can be added to the metadata table by selecting column names listed in the sidebar menu on the left. Table can be filtered to view only the drug arm comparisons/screens of interest in two ways: 1) by performing keyword search using column filters provided below each column header (Please note that it only allows case sensitive search). 2) by using regular expression based filtering of multiple arm IDs using the search box provided in the top right corner of the table. For example, arm IDs - gft-C5, gft-C50 and pub::B1 can be searched by entering \"^gft-C5$|^gft-C50$|^pub::B1$\" in the search box. Here, '^' and '$' represent beginning and end of arm ID, respectively. The operator, '|' means 'or' and searches the table for arm IDs matching gft-C5 or gft-C50 or pub::B1. Metadata table is linked to the table with gene essentiality scores that displays data for only the selected screens/arm comparisons. The filtered table with scores can be viewed by clicking on 'Gene-level scores' tab on the top.", style = "font-family: Ariel; font-size: 14px; color: DarkBlue")),
                          mainPanel(
                            id ="Main",
                            column(5,offset = 5,bsButton("showpanel", "Show/Hide sidebar", type = "toggle", value = TRUE)),
                            DT::dataTableOutput("parpi_metadata_tbl")
                          )
                 ),
                 # Gene-level scores tab
                 tabPanel("Gene-level scores",
                          wellPanel(
                            p("Note:- 'Gene-level scores' tab displays table with gene essentiality scores for CRISPR screens/arm comparisons selected in the metadata table. Analyses methods used for identifying essential genes for CRISPR screens displayed in this table include z-normalisation (median z-scores & average z-scores), log fold change, normz, MAGeCK rank-based scores and STARS rank-based scores. Values in the columns with scores are coloured based on arbitrary thresholds (scores < -2 are in red and scores > 2 are in blue). This threshold only applies to arm comparisons with gene essentiality represented in terms of normz, zscore-avg and log foldchange scores. Based on these thresholds, total number of resistant or sensitising gene hits for each drug target have been summarised in  columns 7-12 of the table. Clicking on a column displays a pop-up window that contains information about that column. Column filters for columns with scores display a slider window for selecting rows with specific range of scores. Files with scores and corresponding metadata can be downloaded by clicking on the download buttons provided in the bottom left corner of the table. Based on this table, genes of interest can be visualised as dot plots in the 'Gene plots' tab. Users are required to enter a comma-separated list of gene symbols in the text box provided. Two separate plots are generated for the same gene depending on the data: one for gene essentiality based on z-score, log foldchange and normz methods (where a negative score represents a negative hit and a positive score represents a positive hit), and the other for gene essentiality based on other scoring methods like MAGeCK and STARS (where the scores are only positive). There is also an option of downloading plots and a table of scores and metadata for selected genes.", style = "font-family: Ariel; font-size: 14px; color: DarkBlue")),
                          DT::dataTableOutput(outputId = "parpi_genez_tbl"),
                          fluidRow(
                            div(style="display:inline-block",downloadButton("downloadGenez", "Download scores")),
                            div(style="display:inline-block",downloadButton("downloadMetadata", "Download metadata"))
                          )),
                 # Gene-level dot plots tab
                 tabPanel("Gene dotplots",
                          fluidRow(
                            #useShinyalert(),			   
                            wellPanel(
                              textInput("geneIDsVEv2", label="Enter gene symbols:", value = "", width = NULL,
                                        placeholder = "e.g., ATM,CDK1,EME1,POLK"))),
                          fluidRow(
                            actionButton(inputId="addGenesVEv2",label="Make dotplots"),
                            downloadButton("downloadDotplotTableVEv2", "Download table with scores and metadata for selected genes")),
                          br(),
                          plotlyOutput("genePlotVEv2",height='600px'),
                          tableOutput("geneTableVEv2")
                 ),
                 # Scatterplots tab
                 tabPanel("Scatterplot",
                          fluidRow(
                            wellPanel(
                              uiOutput("colnames1"),
                              uiOutput("colnames2"))),
                          fluidRow(
                            actionButton(inputId="scatter_plot",label="Make scatterplot")),
                          #downloadButton("downloadDScatterTableVEv2", "Download table with scores and metadata for selected genes")),
                          br(),
                          plotlyOutput("scatterPlot",height=600,width=600)),
                 #tableOutput("scatterplotTable")),
                 tags$head(
                   tags$style(type = 'text/css', 
                              HTML('.navbar { background-color: white;
                              font-weight: normal;
                              font-size: large;
                              font-style: normal;
                              color: grey;
                              font-family: sans-serif;}
                          .navbar-default .navbar-brand{color: black;}
                          .tab-panel{ background-color: black; color: white}
                          .navbar-default .navbar-nav > .active > a, 
                          .navbar-default .navbar-nav > .active > a:focus, 
                          .navbar-default .navbar-nav > .active > a:hover {
                                color: #555;
                                background-color: #fad9e7;
                            }')
                   )
                 )
)

# Server logic
server <- function(input, output, session) {
  #show <- methods::show
  #hide/show metadata sidebar menu
  observeEvent(input$showpanel, {
    
    if(input$showpanel == TRUE) {
      removeCssClass("Main", "col-sm-12")
      addCssClass("Main", "col-sm-8")
      shinyjs::show(id = "Sidebar")
      shinyjs::enable(id = "Sidebar")
    }
    else {
      removeCssClass("Main", "col-sm-8")
      addCssClass("Main", "col-sm-12")
      shinyjs::hide(id = "Sidebar")
    }
  })
  
  #display metadata table
  output$parpi_metadata_tbl  <- DT::renderDataTable({
    DT::datatable(parpi_metadata[, input$show_vars,drop=FALSE],filter="top", selection='none', escape=FALSE, 
                  options = list(
                    autoWidth=TRUE,
                    searching = TRUE,
                    scrollX = TRUE,
                    search = list(regex = TRUE, caseInsensitive = FALSE, search = '(^gft-C[0-9]+$)|(^pub::[A-Z0-9]+$)'),
                    columnDefs = list(list(width = '125px', targets = c(1:15))), #adjust width of all metadata column filters displayed by default
                    initComplete = JS(
                      "function(settings, json) {",
                      "$(this.api().table().header()).css({'background-color': '#e5e7e9', 'color': '#000'});",
                      "}")))
  })
  
  #save user selected score table in an object and use this in required functions later in the code
  filtered_parpi_genez <- reactive({
    sel <- input$parpi_metadata_tbl_rows_all
    drug_to_show <- parpi_metadata[sel,"drug"]
    short_name_to_show <- parpi_metadata[sel,"short_arm_key"]
    f_parpi_genez <- NULL
    f_parpi_genez <- select(parpi_genez, matches(c("gene",annotation_columns,paste0("^",drug_to_show,".*",short_name_to_show,".*|^VE.*",short_name_to_show,".*$|^",short_name_to_show,"$"),"shRNA_hits","siRNA_hits")))
    f_parpi_genez <- select(f_parpi_genez, matches(c("gene",annotation_columns,".*normz.*",".*lfc.*",".*zscore-avg.*",".*MAGECK-MLE.*",".*shRNA_hits",".*siRNA_hits"),ignore.case=TRUE))
    f_parpi_genez <- f_parpi_genez[,!grepl("zscore_med", colnames(f_parpi_genez))]
    
    ## total hits
    total_hits <- NULL
    total_cols <- select(f_parpi_genez,matches(c("Gene",drug_to_show),ignore.case=TRUE))
    total_cols <- select(total_cols,matches(c("Gene",".*normz.*",".*lfc.*",".*zscore-avg.*",".*MAGECK-MLE.*",".*shRNA_hits.*",".*siRNA_hits.*"),ignore.case=TRUE))
    total_cols <- as.data.frame(total_cols)
    total_cols <- column_to_rownames(total_cols,var="gene")
    col_names <- colnames(total_cols)[3:ncol(total_cols)]
    total_cols <- total_cols[,-c(1,2)] 
    total_cols <- as.data.frame(total_cols)
    
    # calculate total hits based on thresholds defined at the beginning of the code
    if(ncol(total_cols)==1){
      colnames(total_cols) <- col_names
      total_cols$total_neg_hits <- ifelse(total_cols[,col_names] < neg.hit.threshold, 1,0)
      total_cols$total_pos_hits <- ifelse(total_cols[,col_names] > pos.hit.threshold, 1,0)
    } else if(ncol(total_cols)>1){
      colnames(total_cols) <- col_names
      total_cols <- total_cols %>% mutate_all(funs(. < neg.hit.threshold)) %>% transmute(total_neg_hits = rowSums(., na.rm=T)) %>%
        bind_cols(total_cols, .)
      
      total_cols <- total_cols <- total_cols %>% mutate_all(funs(. > pos.hit.threshold)) %>% transmute(total_pos_hits = rowSums(., na.rm=T)) %>%
        bind_cols(total_cols, .)
      
      total_hits <- total_cols[,c(ncol(total_cols)-1,ncol(total_cols))]
      
    } else {
      total_hits <- setNames(data.frame(matrix(data=rep(NA, times=2*nrow(f_parpi_genez)),ncol = 2, nrow = nrow(f_parpi_genez))), c("total_neg_hits", "total_pos_hits"))
    } 
    
    #total no. of PARPi hits
    parp_hits <- NULL
    parp_inhibitors_present <- paste0("'",paste(toupper(parp_inhibitors),collapse="' %in% toupper(drug_to_show) | '"),"' %in% toupper(drug_to_show)")
    if(as.factor(parp_inhibitors_present)){
      parp_cols <- select(f_parpi_genez,matches(c("Gene",parp_inhibitors),ignore.case=TRUE))
      parp_cols <- select(parp_cols,matches(c("Gene",".*normz.*",".*lfc.*",".*zscore-avg.*",".*MAGECK-MLE.*"),ignore.case=TRUE))
      parp_cols <- as.data.frame(parp_cols)
      parp_cols <- column_to_rownames(parp_cols,var="gene")
      col_names <- colnames(parp_cols)[3:ncol(parp_cols)]
      parp_cols <- parp_cols[,-c(1,2)] 
      parp_cols <- as.data.frame(parp_cols)
      
      if(ncol(parp_cols)==1){
        colnames(parp_cols) <- col_names
        parp_cols$parp_neg_hits <- ifelse(parp_cols[,col_names] < neg.hit.threshold, 1,0)
        parp_cols$parp_pos_hits <- ifelse(parp_cols[,col_names] > pos.hit.threshold, 1,0)
      } 
      else if(ncol(parp_cols)>1){
        colnames(parp_cols) <- col_names
        parp_cols <- parp_cols %>% mutate_all(funs(. < neg.hit.threshold)) %>% transmute(parp_neg_hits = rowSums(., na.rm=T)) %>%
          bind_cols(parp_cols, .)
        
        parp_cols <- parp_cols <- parp_cols %>% mutate_all(funs(. > pos.hit.threshold)) %>% transmute(parp_pos_hits = rowSums(., na.rm=T)) %>%
          bind_cols(parp_cols, .)
        
        parp_hits <- parp_cols[,c(ncol(parp_cols)-1,ncol(parp_cols))]
        
      } else {
        parp_hits <- setNames(data.frame(matrix(data=rep(NA, times=2*nrow(f_parpi_genez)),ncol = 2, nrow = nrow(f_parpi_genez))), c("parp_neg_hits", "parp_pos_hits"))
      } 
    } else {
      parp_hits <- setNames(data.frame(matrix(data=rep(NA, times=2*nrow(f_parpi_genez)),ncol = 2, nrow = nrow(f_parpi_genez))), c("parp_neg_hits", "parp_pos_hits"))
    }
    
    #total no. of ATRi hits
    atr_hits <- NULL
    atr_inhibitors_present <- paste0("'",paste(toupper(atr_inhibitors),collapse="' %in% toupper(drug_to_show) | '"),"' %in% toupper(drug_to_show)")
    if(as.factor(atr_inhibitors_present)){
      atr_cols <- select(f_parpi_genez,matches(c("Gene",atr_inhibitors),ignore.case=TRUE))
      atr_cols <- select(atr_cols,matches(c("Gene",".*normz.*",".*lfc.*",".*zscore-avg.*",".*MAGECK-MLE.*"),ignore.case=TRUE))
      atr_cols <- as.data.frame(atr_cols)
      atr_cols <- column_to_rownames(atr_cols,var="gene")
      col_names <- colnames(atr_cols)[3:ncol(atr_cols)]
      atr_cols <- atr_cols[,-c(1,2)]
      atr_cols <- as.data.frame(atr_cols)
      
      if(ncol(atr_cols)==1){
        colnames(atr_cols) <- col_names
        atr_cols$atr_neg_hits <- ifelse(atr_cols[,col_names] < neg.hit.threshold, 1,0)
        atr_cols$atr_pos_hits <- ifelse(atr_cols[,col_names] > pos.hit.threshold, 1,0)
        atr_hits <- atr_cols[,c(ncol(atr_cols)-1,ncol(atr_cols))]
      }
      else if(ncol(atr_cols)>1){
        colnames(atr_cols) <- col_names
        atr_cols <- atr_cols %>% mutate_all(funs(. < neg.hit.threshold)) %>% transmute(atr_neg_hits = rowSums(., na.rm=T)) %>%
          bind_cols(atr_cols, .)
        
        atr_cols <- atr_cols %>% mutate_all(funs(. > pos.hit.threshold)) %>% transmute(atr_pos_hits = rowSums(., na.rm=T)) %>%
          bind_cols(atr_cols, .)
        
        atr_hits <- atr_cols[,c(ncol(atr_cols)-1,ncol(atr_cols))]
      } else {
        atr_hits <- setNames(data.frame(matrix(data=rep(NA, times=2*nrow(f_parpi_genez)),ncol = 2, nrow = nrow(f_parpi_genez))), c("atr_neg_hits", "atr_pos_hits"))
      }
    } else {
      atr_hits <- setNames(data.frame(matrix(data=rep(NA, times=2*nrow(f_parpi_genez)),ncol = 2, nrow = nrow(f_parpi_genez))), c("atr_neg_hits", "atr_pos_hits"))
    }
    
    #total no. of Platinum hits
    platinum_hits=NULL
    platinum_inhibitors_present <- paste0("'",paste(toupper(platinum_inhibitors),collapse="' %in% toupper(drug_to_show) | '"),"' %in% toupper(drug_to_show)")
    if(as.factor(platinum_inhibitors_present)){
      platinum_cols <- select(f_parpi_genez,matches(c("Gene",platinum_inhibitors),ignore.case=TRUE))
      platinum_cols <- select(platinum_cols,matches(c("Gene",".*normz.*",".*lfc.*",".*zscore-avg.*",".*MAGECK-MLE.*"),ignore.case=TRUE))
      platinum_cols <- as.data.frame(platinum_cols)
      platinum_cols <- column_to_rownames(platinum_cols,var="gene")
      col_names <- colnames(platinum_cols)[3:ncol(platinum_cols)]
      platinum_cols <- platinum_cols[,-c(1,2)]
      platinum_cols <- as.data.frame(platinum_cols)
      
      if(ncol(platinum_cols)==1){
        colnames(platinum_cols) <- col_names
        platinum_cols$platinum_neg_hits <- ifelse(platinum_cols[,col_names] < neg.hit.threshold, 1,0)
        platinum_cols$platinum_pos_hits <- ifelse(platinum_cols[,col_names] > pos.hit.threshold, 1,0)
        
        platinum_hits <- platinum_cols[,c(ncol(platinum_cols)-1,ncol(platinum_cols))]
      }
      else if(ncol(platinum_cols)>1){
        colnames(platinum_cols) <- col_names
        platinum_cols <- platinum_cols %>% mutate_all(funs(. < neg.hit.threshold)) %>% transmute(platinum_neg_hits = rowSums(., na.rm=T)) %>%
          bind_cols(platinum_cols, .)
        
        platinum_cols <- platinum_cols %>% mutate_all(funs(. > pos.hit.threshold)) %>% transmute(platinum_pos_hits = rowSums(., na.rm=T)) %>%
          bind_cols(platinum_cols, .)   
        
        platinum_hits <- platinum_cols[,c(ncol(platinum_cols)-1,ncol(platinum_cols))]
      } else {
        platinum_hits <- setNames(data.frame(matrix(data=rep(NA, times=2*nrow(f_parpi_genez)),ncol = 2, nrow = nrow(f_parpi_genez))), c("platinum_neg_hits", "platinum_pos_hits"))
      } 
    }else {
      platinum_hits <- setNames(data.frame(matrix(data=rep(NA, times=2*nrow(f_parpi_genez)),ncol = 2, nrow = nrow(f_parpi_genez))), c("platinum_neg_hits", "platinum_pos_hits"))
    }
    
    #add columns with hits to the main gene-level score (genez) table
    filtered_parpi_genez <- NULL
    filtered_parpi_genez <- cbind(
      f_parpi_genez[,1:7],
      total_hits,
      parp_hits,
      atr_hits,
      platinum_hits,
      f_parpi_genez[,8:ncol(f_parpi_genez)]
    )
    rownames(filtered_parpi_genez) <- NULL
    
    filtered_parpi_genez$gene <- paste0("<a href=https://pubmed.ncbi.nlm.nih.gov/?term=", filtered_parpi_genez$gene,"&sort=date target=\"_blank\">", filtered_parpi_genez$gene, "</a>" )
    
    return(filtered_parpi_genez)
  })
  
  #display genez table
  output$parpi_genez_tbl <- DT::renderDataTable({
    filtered_parpi_genez <- filtered_parpi_genez()
    colnames(filtered_parpi_genez) <- gsub(":", " ",colnames(filtered_parpi_genez))
    score_data <- select(filtered_parpi_genez,matches(c("normz","LFC","ZSCORE-AVG","MAGECK-MLE"),ignore.case=FALSE))
    score_cols <- which(colnames(filtered_parpi_genez) %in% colnames(score_data))
    
    DT::datatable(filtered_parpi_genez,
                  filter="top",
                  escape=FALSE, 
                  selection = list(mode='single',target = 'column',selectable=-2),
                  extensions = 'FixedColumns', 
                  options=list(
                    searching=TRUE,
                    autoWidth=TRUE,
                    scrollX = TRUE,
                    columnDefs = list(list(width = '125px', targets = c(2,3,4,5,6))), #adjust width of columns with annotations 
                    fixedColumns = list(leftColumns = 2),
                    initComplete = JS(
                      "function(settings, json) {",
                      "$(this.api().table().header()).css({'background-color': '#e5e7e9', 'color': '#000'});",
                      "}"))) %>%
      formatStyle(columns=score_cols, valueColumns=score_cols,
                  color = JS("value < -2 ? 'red' : value > 2 ? 'blue' : 'black'"))
  })
  
  #download filtered score table
  output$downloadGenez <- downloadHandler(
    filename = function() {
      paste('filtered_genez_data-', Sys.time(), '.csv', sep = '')
    },
    content = function(file){
      filtered_parpi_genez <- filtered_parpi_genez()
      filtered_parpi_genez$gene <- gsub("</a>$","",filtered_parpi_genez$gene)
      filtered_parpi_genez$gene <- gsub("^.*>","",filtered_parpi_genez$gene)
      write.csv(filtered_parpi_genez[input[["parpi_genez_tbl_rows_all"]], ], file,row.names=FALSE)
    }
  )
  
  #download metadata for filtered genez table
  output$downloadMetadata <- downloadHandler(
    filename = function() {
      paste('metadata-', Sys.time(), '.csv', sep = '')
    },
    content = function(file){
      filtered_parpi_genez <- filtered_parpi_genez()   
      i=NULL
      short_colnames <- NULL
      for(i in 1:ncol(filtered_parpi_genez)){
        if(str_contains(colnames(filtered_parpi_genez)[i], "gft-C")){
          shortCols <- gsub("^.*:","",colnames(filtered_parpi_genez)[i])
        } else{
          shortCols <- colnames(filtered_parpi_genez)[i]
        }
        short_colnames <- rbind(short_colnames,shortCols)
      }
      columns_to_drop <- c("data","short_arm_key","screen_ID","analysis_workflow","analysis_report","geneZ_file","guideZ_file")
      metadata_table <- parpi_metadata[which(parpi_metadata$arm_ID %in% short_colnames),!(names(parpi_metadata) %in% columns_to_drop)]
      write.csv(metadata_table, file, row.names=FALSE)
    }
  )
  
  #popup genez column details - show column info when each column is clicked
  observeEvent(input$parpi_genez_tbl_columns_selected, {
    selected_cols <- input$parpi_genez_tbl_columns_selected
    filtered_parpi_genez <- filtered_parpi_genez()
    drug_columns <- colnames(filtered_parpi_genez)[!grepl(c("^VE."), colnames(filtered_parpi_genez))]
    drug_columns <- drug_columns[!grepl(c("^gene"), drug_columns)]
    drug_columns <- drug_columns[!grepl(c("^DDR"), drug_columns)]
    drug_columns <- drug_columns[!grepl(c("^Process"), drug_columns)]
    drug_columns <- drug_columns[!grepl(c("^CGC"), drug_columns)]
    drug_columns <- drug_columns[!grepl(c("^CEG"), drug_columns)]
    drug_columns <- drug_columns[!grepl(c("hits"), drug_columns)]
    if (length(selected_cols) > 0){
      filtered_data <- filtered_parpi_genez[,selected_cols]
      filtered_data <- as.data.frame(filtered_data)
      colnames(filtered_data) <- colnames(filtered_parpi_genez)[selected_cols]
      rownames(filtered_data) <- NULL
      selected_arm=NULL
      if(colnames(filtered_data) == "gene_key"){
        selected_arm = "Gene symbol:Entrez gene ID"
      }
      if(colnames(filtered_data) == "gene_name"){
        selected_arm = "Gene name"
      }
      if(colnames(filtered_data) == "gene"){
        selected_arm = "Gene symbol"
      }
      if(colnames(filtered_data) == "DDR"){
        #selected_arm = "Process"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>DNA Damage Response annotation for genes</label>"))
          ),
          easyClose = TRUE
        ))
      }
      if(colnames(filtered_data) == "Process"){
        #selected_arm = "Process"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>Process for DNA Damage Response</label>"))
          ),
          easyClose = TRUE
        ))
      }
      if(colnames(filtered_data) == "CGC"){
        #selected_arm = "CGC"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>Cancer Gene Census annotation</label>"))
          ),
          easyClose = TRUE
        ))
      }
      if(colnames(filtered_data) == "CEG"){
        #selected_arm = "CEG"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>Core Essential Gene annotation</label>"))
          ),
          easyClose = TRUE
        ))
      }
      if(str_contains(colnames(filtered_data), "_neg_hits")){
        #selected_arm = "CEG"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>Total number of sensitising hits identified for a drug target (cutoff < -2; scoring methods include zscore-avg, log fold change and normZ)</label>"))
          ),
          easyClose = TRUE
        ))
      }
      if(str_contains(colnames(filtered_data), "_pos_hits")){
        #selected_arm = "CEG"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>Total number of resistance causing hits identified for a drug target (cutoff > 2; scoring methods include zscore-avg, log fold change and normZ)</label>"))
          ),
          easyClose = TRUE
        ))
      }
      if(str_contains(colnames(filtered_data), "MCF7_Olaparib_1.5uM_sf80:human_shRNA_lib:shRNA_hits")){
        #selected_arm = "CEG"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>This column tells whether a gene was a (positive/negative) hit or not in the shRNA screen. Screen details for hits included in this column: Cell line=MCF7;library=OpenBiosystems GIPZ human shRNA library; Drug=Olaparib; dose=1.5uM; SF=80; hit threshold=Drug effect Z score +/- 1.96 with 2 or more shRNA/gene; Publication=PMID: 24240700</label>"))
          ),
          easyClose = TRUE
        ))
      }
      if(str_contains(colnames(filtered_data), "CAL51_Olaparib_1.0uM_sf80:siKINOME_smartpool_lib:siRNA_hits")){
        #selected_arm = "CEG"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>This column tells whether a gene was a (positive/negative) hit or not in the siRNA screen. Screen details for hits included in this column: Cell line=CAL51;library=Dharmacon human siKINOME Smartpool; Drug=Olaparib; dose=1.0uM; SF=80; hit threshold=Drug effect Z score +/- 3; Publication=PMID:18388863</label>"))
          ),
          easyClose = TRUE
        ))
      }
      if(str_contains(colnames(filtered_data), "CAL51_Olaparib_1.0uM_sf80:DNA_damage_response_lib:siRNA_hits")){
        #selected_arm = "CEG"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>This column tells whether a gene was a (positive/negative) hit or not in the siRNA screen. Screen details for hits included in this column: Cell line=CAL51;library=Qiagen human DNA Repair siRNA Set V1.0 siRNA library; Drug=Olaparib; dose=1.0uM; SF=80; hit threshold=Drug effect Z score +/- 3 with >1 siRNA; Publication=PMID: 18832051</label>"))
          ),
          easyClose = TRUE
        ))
      }
      if(str_contains(colnames(filtered_data), "CAL51_Olaparib_1.0uM_sf80:metabolome_smartpool_lib:siRNA_hits")){
        #selected_arm = "CEG"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>This column tells whether a gene was a (positive/negative) hit or not in the siRNA screen. Screen details for hits included in this column: Cell line=CAL51;library=Dharmacon human siDNAmetabolome (custom) Smartpool; Drug=Olaparib; dose=1.0uM; SF=80; hit threshold=Drug effect Z score +/- 3; Publication=PMID: 22933245</label>"))
          ),
          easyClose = TRUE
        ))
      }
      if(str_contains(colnames(filtered_data), "VE.")){
        # selected_arm = "Viability effect"
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>Viability Effect</label>"))
          ),
          easyClose = TRUE
        ))
      } 
      if(str_contains(colnames(filtered_data), drug_columns, logic = "or")){
        metadata_row <- NULL
        if(str_contains(colnames(filtered_data), "gft-C")){
          short_drug_col <- sub(".*:","",colnames(filtered_data))
          metadata_row <- which(parpi_metadata$arm_ID  == short_drug_col)
        } else {
          short_drug_col <- colnames(filtered_data)
          metadata_row <- which(parpi_metadata$short_arm_key == short_drug_col)
        }
        selected_arm <- parpi_metadata[metadata_row, c("cell_line","cas9_type","drug","drug_conc_nM","SF","library","analysis_pipeline","data")]
        showModal(modalDialog(
          title = "column details:",
          size="s",
          wellPanel(
            fluidRow(HTML("<label>Cell line:</label>"), selected_arm$cell_line),
            fluidRow(HTML("<label>cas9 type:</label>"), selected_arm$cas9_type) ,
            fluidRow(HTML("<label>drug:</label>"), selected_arm$drug) ,
            fluidRow(HTML("<label>dose (nM):</label>"), selected_arm$drug_conc_nM),
            fluidRow(HTML("<label>SF:</label>"), selected_arm$SF) ,
            fluidRow(HTML("<label>library:</label>"),selected_arm$library),
            fluidRow(HTML("<label>analysis pipeline:</label>"), selected_arm$analysis_pipeline),
            fluidRow(HTML("<label>data:</label>"), selected_arm$data)
          ),
          easyClose = TRUE
        ))
      }
    }
  })
  
  #Version2 - prepare data for making dotplots with VE based on user selected genes
  filtered_genez_renamed_finalv2 <- eventReactive(input$addGenesVEv2, {
    #get user input - list of genes to plot
    geneNames <- input$geneIDsVEv2
    geneNames <- gsub("\\s","",geneNames)
    geneNames_split <- as.vector(str_split(geneNames,",",simplify=TRUE))
    sel <- input$parpi_metadata_tbl_rows_all
    drug_to_show <- parpi_metadata[sel,"drug"]
    drug_to_show <- as.character(drug_to_show)
    short_name_to_show <- parpi_metadata[sel,"short_arm_key"]
    short_name_to_show <- as.character(short_name_to_show)
    f_parpi_genez <- NULL
    f_parpi_genez <- select(parpi_genez, matches(c("gene",annotation_columns,paste0("^",drug_to_show,".*",short_name_to_show,".*|^VE.*",short_name_to_show,"$|^",short_name_to_show,"$"))))
    
    #alert popup if gene symbols entered by user don't match those in the genez table
    genes_absent <- geneNames_split[!(geneNames_split %in% f_parpi_genez$gene)]
    if(length(genes_absent) > 0){
      shinyalert(
        title = "Please check gene symbols",
        type = "warning",
        showConfirmButton = TRUE,
        confirmButtonCol = "white",
        timer=2000
      )
    }
    
    #table with selected DE columns 
    genez_de_cols <- select(f_parpi_genez, matches(c("^gene$",paste0("^",drug_to_show,".*",short_name_to_show,".*|^",short_name_to_show,"$"))))
    genez_de_cols <- select(genez_de_cols,matches(c("gene","normz","lfc","mageck","stars","zscore-avg"),ignore.case=TRUE))
    
    #DE table with genes of interest 
    genez_de_cols.g <- genez_de_cols[which(genez_de_cols$gene %in% geneNames_split),]
    rownames(genez_de_cols.g) <- NULL
    genez_de_cols.g <- as.data.frame(genez_de_cols.g)
    
    genez_de_cols.g.long <- melt(genez_de_cols.g,id.var="gene")
    genez_de_cols.g.long.2 <- genez_de_cols.g.long
    genez_de_cols.g.long.2$variable <- sub(":gft-C.*|:pub::.*","",genez_de_cols.g.long.2$variable)
    genez_de_cols.g.long.2$variable <- sub("^.*?\\.","",genez_de_cols.g.long.2$variable)
    
    #table with selected VE columns
    genez_ve_cols <- select(parpi_genez, matches(c("^gene$",paste0("^VE.*",short_name_to_show,"$"))))
    genez_ve_cols <- select(genez_ve_cols,matches(c("gene","normz","lfc","mageck","stars","zscore-avg"),ignore.case=FALSE))
    if(ncol(genez_ve_cols) < 2){
      genez_ve_cols[,2] <- 0
      colnames(genez_ve_cols)[2] <- "VE"
      
    }
    
    #add corresponding VE scores=0 for drug columns from published data selected by user
    pub_cols <- select(genez_de_cols, matches(":pub::"))
    if(ncol(pub_cols) > 0){
      colnames(pub_cols) <- sub(":pub::.*","",colnames(pub_cols))
      pub_cols_ve <- paste0("VE.",colnames(pub_cols))
      pub_cols_ve_number <- ncol(pub_cols)
      make_pub_cols_ve <- matrix(0,nrow=nrow(genez_de_cols),ncol=pub_cols_ve_number)
      colnames(make_pub_cols_ve) <- pub_cols_ve
      make_pub_cols_ve <- data.frame(make_pub_cols_ve)
      genez_ve_cols <- cbind(genez_ve_cols,make_pub_cols_ve)
    }
    
    #VE table with genes of interest
    genez_ve_cols.g <- genez_ve_cols[which(genez_ve_cols$gene %in% geneNames_split),]
    rownames(genez_ve_cols.g) <- NULL
    
    genez_ve_cols.g <- as.data.frame(genez_ve_cols.g)
    cols.ve1 <- sub("^.*?\\.","",colnames(genez_ve_cols.g))
    genez_ve_cols.g.2 <- genez_ve_cols.g
    colnames(genez_ve_cols.g.2) <- cols.ve1
    genez_ve_cols.g.2.long <- melt(genez_ve_cols.g.2,id.var="gene")
    
    genez.de.ve <- left_join(genez_de_cols.g.long.2,genez_ve_cols.g.2.long,by=c("gene"="gene","variable"="variable"))
    
    genez.de.ve$variable <- genez_de_cols.g.long$variable
    genez.de.ve <- mutate(genez.de.ve, value = paste(value.y,value.x, sep = ';'))
    genez.de.ve <- select(genez.de.ve, -c("value.x","value.y"))
    genez.de.ve.wide <- genez.de.ve %>% pivot_wider(names_from = variable, values_from = value)
    
    filtered_parpi_genez.g <- genez.de.ve.wide
    
    filtered_parpi_genez.g <- filtered_parpi_genez.g[,colSums(is.na(filtered_parpi_genez.g))<nrow(filtered_parpi_genez.g)]
    rownames(filtered_parpi_genez.g) <- NULL
    filtered_parpi_genez.g <- as.data.frame(filtered_parpi_genez.g)
    filtered_parpi_genez.g <- column_to_rownames(filtered_parpi_genez.g,var="gene")
    filtered_parpi_genez.g.t <- as.data.frame(t(filtered_parpi_genez.g))
    
    #add arm ID to genez table with g.o.i
    i=NULL
    for(i in 1:nrow(filtered_parpi_genez.g.t)){
      filtered_parpi_genez.g.t$armID[i] <- ifelse(str_contains(rownames(filtered_parpi_genez.g.t)[i], "gft-C"), gsub(".*:","",rownames(filtered_parpi_genez.g.t)[i]),rownames(filtered_parpi_genez.g.t)[i])
    }
    
    #join selected rows of metadata and  g.o.i genez rows with common arm ID
    filtered_metadata <- parpi_metadata[sel,]
    complete_filtered_data <- right_join(filtered_metadata,filtered_parpi_genez.g.t,by=c("arm_ID"="armID")) 
    complete_filtered_data$short_arm_key <- as.character(as.factor(complete_filtered_data$short_arm_key))
    
    #metadata to display in interactive dotplots; metadata columns: 1= arm_ID,8=cell_line,11=cas9_type,17=timepoint,27=analysis_pipeline,28=analysis_method,32=short_arm_key,20=drug_conc_nM,18=drug PLUS columns with scores for g.o.i
    #short_filtered_data <- complete_filtered_data[,c(1,8,11,17,28,29,33,21,18,34:ncol(complete_filtered_data))]
    short_filtered_data <- complete_filtered_data[,c(1,8,11,17,28,29,33,21,18,38:ncol(complete_filtered_data))]
    
    #make long table with only score columns and add metadata columns afterwards
    short_filtered_data.m <- (reshape2::melt(short_filtered_data[,-c(1:8)],id.var="drug"))
    short_filtered_data.m$variable <- as.character(as.factor(short_filtered_data.m$variable))
    columns_to_replicate <- short_filtered_data[,1:8]
    repeat_times <- length(geneNames_split)
    columns_to_add <- columns_to_replicate[rep(seq_len(nrow(columns_to_replicate)), repeat_times), ]
    short_filtered_data.table <- data.frame(short_filtered_data.m,columns_to_add)
    rownames(short_filtered_data.table) <- NULL
    
    #split DE and VE scores into separate columns
    ##filtered_genez_renamed  <- short_filtered_data.table %>% rename(gene=variable, VE_DE=value)
    score <- data.frame(do.call('rbind', strsplit(short_filtered_data.table$value,';',fixed=TRUE)))
    #score <- score %>% rename(VE.score=X1, DE.score=X2)
    score$X1 <- as.numeric(as.character(score$X1))
    score$X2 <- as.numeric(as.character(score$X2))
    filtered_genez_renamed_finalv2  <- cbind(
      short_filtered_data.table$drug,
      short_filtered_data.table$variable,
      score$X1,
      score$X2,
      short_filtered_data.table$arm_ID,
      short_filtered_data.table$cell_line,
      short_filtered_data.table$cas9_type,
      short_filtered_data.table$timepoint,
      short_filtered_data.table$analysis_pipeline,
      short_filtered_data.table$analysis_method,
      short_filtered_data.table$short_arm_key,
      short_filtered_data.table$drug_conc_nM
    )
    colnames(filtered_genez_renamed_finalv2) <- c(
      "drug",
      "gene",
      "VE.score",
      "DE.score",
      "arm_ID",
      "cell_line",
      "cas9_type",
      "timepoint",
      "analysis_pipeline",
      "analysis_method",
      "short_arm_key",
      "drug_conc_nM"
    )
    filtered_genez_renamed_finalv2 <- as.data.frame(filtered_genez_renamed_finalv2,stringsAsFactors=FALSE)
    filtered_genez_renamed_finalv2$VE.score <- as.numeric(as.character(filtered_genez_renamed_finalv2$VE.score))
    filtered_genez_renamed_finalv2$DE.score <- as.numeric(as.character(filtered_genez_renamed_finalv2$DE.score))
    filtered_genez_renamed_finalv2$drug_conc_nM <- as.numeric(as.character(filtered_genez_renamed_finalv2$drug_conc_nM))
    
    #replace NAs with 0 in VE.score column
    filtered_genez_renamed_finalv2$VE.score <- ifelse(is.na(filtered_genez_renamed_finalv2$VE.score), 0,filtered_genez_renamed_finalv2$VE.score)
    
    #remove only those rows that have NA in the last column with scores
    completeFun <- function(data, desiredCols) {
      completeVec <- complete.cases(data[, desiredCols])
      return(data[completeVec, ])
    }
    filtered_genez_renamed_finalv2 <- completeFun(filtered_genez_renamed_finalv2, "DE.score")
    
    #categorize analysis methods as rank-based and non-rank based
    m=NULL
    for(m in 1:nrow(filtered_genez_renamed_finalv2)){
      filtered_genez_renamed_finalv2$method[m] <- ifelse(filtered_genez_renamed_finalv2$analysis_method[m] == "STARS" | filtered_genez_renamed_finalv2$analysis_method[m] == "MAGeCK-RRA", "[scoring method: rank-based scores]","[scoring method: Normz/Zscore/LFC/Beta-score]")
    }  
    
    return(filtered_genez_renamed_finalv2)
  })
  
  output$genePlotVEv2 <- renderPlotly({
    req(nrow(filtered_genez_renamed_finalv2()) > 0)
    filtered_genez_renamed_finalv2 <- filtered_genez_renamed_finalv2()
    
    #subset data for geom_hline
    sub.set1 <- filtered_genez_renamed_finalv2 %>%
      filter(method == "[scoring method: Normz/Zscore/LFC/Beta-score]")
    sub.set1 <- sub.set1 %>% mutate(y_val = rep(-2,times=nrow(sub.set1)))
    rownames(sub.set1) <- NULL
    sub.set2 <- filtered_genez_renamed_finalv2 %>%
      filter(method == "[scoring method: Normz/Zscore/LFC/Beta-score]")
    sub.set2 <- sub.set2 %>% mutate(y_val = rep(2,times=nrow(sub.set2)))
    rownames(sub.set2) <- NULL
    sub.set <- bind_rows(sub.set1,sub.set2)
    
    #point shape
    filtered_genez_renamed_finalv2$point.shape <- ifelse(filtered_genez_renamed_finalv2$VE.score > 3 | filtered_genez_renamed_finalv2$VE.score < -3,19,1)
    
    # adjust plot width according to number of gene to be plotted
    if(length(unique(filtered_genez_renamed_finalv2$gene)) <= 2){
      plot_height <- 500
    }
    if(length(unique(filtered_genez_renamed_finalv2$gene)) > 2 & length(unique(filtered_genez_renamed_finalv2$gene)) <= 4){
      plot_height <- 1000
    }
    if(length(unique(filtered_genez_renamed_finalv2$gene)) > 4 & length(unique(filtered_genez_renamed_finalv2$gene)) <= 6){
      plot_height <- 1500
    }
    if(length(unique(filtered_genez_renamed_finalv2$gene)) > 6 & length(unique(filtered_genez_renamed_finalv2$gene)) <= 8){
      plot_height <- 2000
    }
    if(length(unique(filtered_genez_renamed_finalv2$gene)) > 8){
      plot_height <- 2800
    }
    
    #make dotplots
    g.plot.ve <- ggplot(data=filtered_genez_renamed_finalv2,aes(x=drug,y=DE.score,shape=factor(point.shape),colour=short_arm_key,text=paste("VE.score: ",VE.score,"\ncell line: ",cell_line,"\ntimepoint: ",timepoint,"\ndose (nM): ",drug_conc_nM,"\ncas9 type: ",cas9_type,"\nanalysis pipeline: ",analysis_pipeline,"\nanalysis method: ",analysis_method,"\narm ID: ",arm_ID))) +
      geom_beeswarm(alpha = 0.6,size = 2.5) +
      scale_shape_manual(values=c(19,1)) +
      facet_wrap(~gene+method,ncol=2, scales = "free_x") +
      theme(panel.spacing = unit(2, "lines")) +
      theme(strip.text.x = element_text( margin = margin(b = 8, t = 8))) +
      xlab("") +
      ylab("DE score") +
      coord_flip() +
      theme(legend.position = "none") +
      theme(panel.spacing.y = unit(1.0, "lines")) +
      theme(panel.spacing.x = unit(0.1, "lines"))
    gene.plot.ve <- g.plot.ve + geom_hline(data=sub.set, aes(yintercept=y_val),linetype="dotted",col="darkgrey")
    ggplotly(gene.plot.ve,height=plot_height,toolkit="text")
  })
  
  #download table with genes selected for making DE/VE dotplots 
  output$downloadDotplotTableVEv2 <- downloadHandler(
    filename = function() {
      paste('v2_dotplot_DE.VE.data_table-', Sys.time(), '.csv', sep = '')
    },
    content = function(file){
      req(input$addGenesVEv2 > 0)
      output <- filtered_genez_renamed_finalv2()
      output <- output[,-c(ncol(output)-1,ncol(output))]
      write.csv(output, file,row.names=FALSE)
    })
  
  #enable users to select columns to plot from dropdown column lists
  #user selected column for plotting x-axis
  output$colnames1 <- renderUI({
    filtered_parpi_genez <- filtered_parpi_genez()
    filtered_selectInput_data <- select(filtered_parpi_genez,matches(c(":gft-C",":pub::","^VE."),ignore.case = TRUE))
    filtered_selectInput_data <- select(filtered_parpi_genez,matches(c("normz","LFC","ZSCORE-AVG","MAGECK-MLE"),ignore.case = TRUE))
    filtered_colnames <- colnames(filtered_selectInput_data)
    tagList(
      selectInput(
        "colvalue1",
        selectize=TRUE,
        label = " Screen arm comparison (x-axis):", 
        choices = c("Please select",filtered_colnames)
      )
    )
   });
  #user selected column for plotting y-axis
  output$colnames2 <- renderUI({
    filtered_parpi_genez <- filtered_parpi_genez()
    filtered_selectInput_data <- select(filtered_parpi_genez,matches(c(":gft-C",":pub::","^VE."),ignore.case = TRUE))
    filtered_selectInput_data <- select(filtered_parpi_genez,matches(c("normz","LFC","ZSCORE-AVG","MAGECK-MLE"),ignore.case = TRUE))
    filtered_colnames <- colnames(filtered_selectInput_data)
    tagList(
      selectInput(
        "colvalue2",
        selectize=TRUE,
        label = " Screen arm comparison (y-axis):", 
        choices = c("Please select",filtered_colnames)
      )
    )
  });
  
  toListen1<- reactive({
    input$colvalue1
  })
  toListen2<- reactive({
    input$colvalue2
  })
  
  #prepare data for making scatterplots
  selected_cols <- eventReactive(input$scatter_plot, {
    #retrieve selective column names for plotting
    colval1 <- toListen1()
    colval2 <- toListen2()
    sel <- input$parpi_metadata_tbl_rows_all
    drug_to_show <- parpi_metadata[sel,"drug"]
    drug_to_show <- as.character(drug_to_show)
    short_name_to_show <- parpi_metadata[sel,"short_arm_key"]
    short_name_to_show <- as.character(short_name_to_show)
    f_parpi_genez <- NULL
    f_parpi_genez <- select(parpi_genez, matches(c("gene",annotation_columns,paste0("^",drug_to_show,".*",short_name_to_show,".*|^VE.*",short_name_to_show,"$|^",short_name_to_show,"$"))))
    selected_deve_cols <- select(f_parpi_genez,matches(c("^gene$","normz","lfc","mageck-mle","zscore-avg"),ignore.case=TRUE))
    selected_deve_cols <- selected_deve_cols[,c("gene",colval1,colval2)]
    rownames(selected_deve_cols) <- NULL
    selected_deve_cols <- as.data.frame(selected_deve_cols)
    short_deve_colnames <- sub("^.*?:","",colnames(selected_deve_cols)[2:3])
    short_deve_colnames <- sub(":.*$","",short_deve_colnames)
    
    selected_ve_cols <- select(parpi_genez, matches(paste0("^VE.*",short_deve_colnames,"$")))
    selected_ve_cols <- select(selected_ve_cols,matches(c("normz","lfc","mageck","stars","zscore-avg"),ignore.case=FALSE))
    if(ncol(selected_ve_cols) < 1){
      selected_ve_cols$VE <- 0
      colnames(selected_ve_cols) <- "VE"
    } else if (ncol(selected_ve_cols) == 1) {
      colnames(selected_ve_cols) <- "VE1"
    } else {
      colnames(selected_ve_cols) <- c("VE1","VE2")
    }
    selected_cols <- NULL
    selected_cols <- cbind(selected_deve_cols,selected_ve_cols)
    
    cols <- colnames(selected_cols)
    cols <- gsub(":","_",cols)
    cols <- gsub("-","\\.",cols)
    colnames(selected_cols) <- cols
    
    return(selected_cols)
  })
  
  #make scatterplots
  output$scatterPlot <- renderPlotly({
    selected_cols <- selected_cols()
    s.plot <- ggplot(selected_cols, aes_string(x=colnames(selected_cols)[2], y=colnames(selected_cols)[3],Gene=colnames(selected_cols)[1],VE=colnames(selected_cols)[4])) +
      geom_vline(xintercept = 0, linetype = 2, col = "grey50") +
      geom_hline(yintercept = 0, linetype = 2, col = "grey50") +
      geom_point(color="skyblue3", alpha = 0.5, pch = 1) +
      geom_smooth(inherit.aes=FALSE,aes_string(x=colnames(selected_cols)[2], y=colnames(selected_cols)[3]),method="lm",se=FALSE,show.legend = FALSE,size=0.2,col="maroon")
    ggplotly(s.plot, height = 600, tooltip=c("Gene",colnames(selected_cols)[1],colnames(selected_cols)[2],colnames(selected_cols)[3],"VE"))
  })
}
shinyApp(ui, server)



##END##