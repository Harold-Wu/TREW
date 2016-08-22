
library(shiny)
library(sqldf)
library(DT)


db <<- dbConnect(SQLite(), dbname="TREW_1.2.db")
Sites_Info <<- dbReadTable(db, "Sites_Info")
Unique_gene <<- unique(Sites_Info$Gene_ID)
# Define UI for application that draws a histogram
ui <- shinyUI(navbarPage("TREW", inverse = TRUE, collapsible = FALSE,
                         tabPanel("Home", icon=icon("home"), 
                                  navlistPanel("Introduction", widths = c(2, 8),
                                               tabPanel("TREW",
                                                        mainPanel(
                                                          h2("Epitranscriptomic targets of 
                                                             RNA modification readers, erasers 
                                                             and writers"),
                                                          h4("abcdefghijklmnopqrstuvwxyz")
                                                          )),
                                               "Background",
                                               tabPanel("Readers",
                                                        mainPanel(
                                                          h2("Readers"),
                                                          h4("abcdefghijklmnopqrstuvwxyz")
                                                        )),
                                               tabPanel("Erasers",
                                                        mainPanel(
                                                          h2("Erasers"),
                                                          h4("abcdefghijklmnopqrstuvwxyz")
                                                        )),
                                               tabPanel("Writers",
                                                        mainPanel(
                                                          h2("Writers"),
                                                          h4("abcdefghijklmnopqrstuvwxyz")
                                                        )))
                                  ),
                         tabPanel("Table & Figure",
                                  sidebarLayout(
                                    sidebarPanel(width = 3,
                                                 selectInput("Gene_choice", "Gene name", Unique_gene, width = 200),
                                                 selectInput("Species_choice", "Species", 
                                                             c("Drosophila melanogaster","Homo sapiens",  "Mus musculus"), width = 200),
                                                 selectInput("Regulation_choice", "Regulation type", 
                                                             c("Regulator", "Regulatee", "Both of above"), width = 200),
                                                 selectInput("Modification_choice", "Modification type", 
                                                             c("m6A","m1A", "m5C", "Psi", "All of above"), width = 200),
                                                 selectInput("Content_choice", "Output table content", c("Default", "Complete"), 
                                                             width = 200),
                                                 helpText("The default content consists of: chromosome, strand, range start site, 
                                                          range width, consistency, genome assembly type, modification type, technique used,
                                                          target and cell line"),
                                                 submitButton(text = "Search")
                                                 ),
                                    mainPanel(DT::dataTableOutput("x3"),
                                              #textOutput("No_gene"),
                                              #downloadLink('Download_Table', 'Download')
                                              verbatimTextOutput('x4')
                                    )
                                    )),
                         tabPanel("Help")
                         
                         ))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  #db <- dbConnect(SQLite(), dbname="TREW_1.2.db")
  Genome_Location <- dbReadTable(db, "Genome_Location")
  #Sites_Info <- dbReadTable(db, "Sites_Info")
  Source_Info <- dbReadTable(db, "Source_Info")
  #Unique_gene <- unique(Sites_Info$Gene_ID)
  # #=====================================================================================================================
  #1.5 Generate the indexes that can be re-used
  idx_2to1 = table(Genome_Location$Meth_Site_ID)
  idx_3to2 = table(Sites_Info$Source_ID)[Source_Info$DataSet_ID]
  idx_3to1 = table(rep(Sites_Info$Source_ID,idx_2to1))[Source_Info$DataSet_ID]
  
  # #=====================================================================================================================
  #2. Function Input: Several constraints
  #Function Output: Detail or concise.
  GenerateTable <- function(
    Species = NULL,
    Target = NULL,
    Modification = NULL,
    Gene_ID = NULL,
    Return_Detail = FALSE
  )
  {
    # #=====================================================================================================================
    #3. Realize the filter from small table to larger table
    idx2 = !vector(length = dim(Sites_Info)[1])
    idx3 = !vector(length = dim(Source_Info)[1])
    if (!is.null(Species)) {
      idx3 <- Source_Info$Species == Species
    }else{}
    
    if (Target == "Regulator") {
      idx3 <- idx3 & (Source_Info$Target != "YTHDF1" & Source_Info$Target != "YTHDF2")
    }else{
      if (Target == "Regulatee") {
        idx3 <- idx3 & (Source_Info$Target == "YTHDF1" | Source_Info$Target == "YTHDF2")
      }else{}
    }
    
    if (!is.null(Modification) & Modification != "All of above") {
      idx3 <- idx3 & (Source_Info$Modification == Modification)
    }else{}
    
    if (!is.null(Gene_ID)) {
      idx2 <- (Sites_Info$Gene_ID == Gene_ID) & rep(idx3,idx_3to2)
    }else{}
    idx2[which(is.na(idx2))] <- FALSE
    
    idx1 <- rep(idx2,idx_2to1)
    
    table2_tmp <- Sites_Info[rep(which(idx2),idx_2to1[which(idx2)]),]
    
    Full_table <- cbind(Genome_Location[which(idx1),],
                        Sites_Info[rep(which(idx2),idx_2to1[which(idx2)]),],
                        Source_Info[rep(1:dim(Source_Info)[1],idx_3to1)[which(idx1)],])
    if (Return_Detail == "Default") {
      Full_table <- Full_table[, c("Range_start", "Range_width", "Strand", "Chromosome",
                                   "Consistency", "Genome_assembly", "Modification",
                                   "Technique", "Target", "Cell_line")]
    }else{}
    
    DT::datatable(Full_table, style = "bootstrap", extensions = 'Buttons', 
                  options = list(dom = 'Bfrtip', buttons = list('copy', 'csv', 'excel', 'pdf'), 
                                 scrollX = TRUE
                  ))
  }
  
  output$x3 <- DT::renderDataTable(
    GenerateTable(Species = input$Species_choice,
                  Target = input$Regulation_choice,
                  Modification = input$Modification_choice,
                  Gene_ID = input$Gene_choice,
                  Return_Detail = input$Content_choice), server = TRUE
  )
  
  output$x4 <- renderPrint({
    input$x3_rows_selected
  }
  
  )
  
  # output$Download_Table <- downloadHandler(
  #   filename = function() {
  #     paste('data-', Sys.Date(), '.csv', sep='')
  #   },
  #   content = function(file) {
  #     write.csv(Final_Table, file)
  #   }
  # )
  # 
  # 
  # output$No_gene <- renderText({
  #   Empty <- paste("No item found.")
  #   Empty
  # })
  
  
})

# Run the application 
shinyApp(ui = ui, server = server)

