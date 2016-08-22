library(shiny)

# Define UI for application that draws a histogram
ui <- shinyUI(navbarPage("TREW", inverse = TRUE, collapsable = FALSE,
                         tabPanel("Introduction",
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
                                                 selectInput("species_choice", "Species", 
                                                             c("Homo sapiens", "Drosophila melanogaster", 
                                                               "Mus musculus"), width = 200),
                                                 helpText("hahahaha"),
                                                 selectInput("target_choice", "Target", 
                                                             c("dTet", "ALKBH", "Fto", "KIAA1429", 
                                                               "METTL14", "METTL3", "WTAP", "YTHDF2",
                                                               "hPUS1"),width = 200),
                                                 checkboxGroupInput("item_location", "Genome location items choose", 
                                                                    c("Range_start", "Range_width", "Strand", 
                                                                      "Chromosome")),
                                                 checkboxGroupInput("item_sites", "Sites items choose", 
                                                                    c("Diff_p_value", "Diff_fdr", "Doff_log2FoldChange", 
                                                                      "Gene_ID", "Source_ID", "Consistency", "Overlap_UTR5", 
                                                                      "Overlap_CDS", "Overlap_UTR3", "Distance_ConsensusMotif", 
                                                                      "Distance_StartCodon", "Distance_StopCodon")),
                                                 checkboxGroupInput("item_source", "Source items choose", 
                                                                    c("Genome_assembly", "Technique", "Perturbation", 
                                                                      "Date_of_process", "Paper", "Cell_line", "Treatment", 
                                                                      "LiftOver", "Computation_pepline")),
                                                 checkboxGroupInput("item_raw", "Raw data items choose", 
                                                                    c("GEO_RUN", "IP_Input", "Genotype", "Replicate")),
                                                 submitButton(text = "Run")
                                    ),
                                    mainPanel(tableOutput("Final_Table"))
                                  )),
                         tabPanel("Help")
                         
                                  ))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  library(sqldf)
  #Load the database.
  con <- dbConnect(RSQLite::SQLite(), dbname="TREW_1.1.db")
  d3_choice <- sqldf("SElECT * from Source_Info where Species = 'input$species_choice' and
                     Target = 'input$target_choice'",connection = con)
  d1 <- dbReadTable(con, "Genome_Location")
  d2 <- dbReadTable(con, "Sites_Info")
  d3 <- sqldf("SElECT * from Source_Info where Species = 'input$species_choice' and
              Target = 'input$target_choice'",connection = con)
  d4 <- dbReadTable(con, "Raw_data_records")
  
  output$Final_Table <- renderTable({
    #Extract colomu.
    extract1 <- d1[,input$item_location]
    extract1[, length(input$item_location) + 1] <- d1[, "Meth_Site_ID"]
    names(extract1)[length(input$item_location) + 1] <- "Meth_Site_ID"
    
    extract2 <- d2[,input$item_sites]
    extract2[, length(input$item_sites) + 1] <- d2[, "Methylation_ID"]
    names(extract2)[length(input$item_sites) + 1] <- "Methylation_ID"
    extract2[, length(input$item_sites) + 2] <- d2[, "Source_ID"]
    names(extract2)[length(input$item_sites) + 2] <- "Source_ID"
    
    extract3 <- d3[,input$item_source]
    extract3[, length(input$item_source) + 1] <- d3[, "DataSet_ID"]
    names(extract3)[length(input$item_source) + 1] <- "DataSet_ID"
    

    
    #Merge the extract table
    m_e1_e2 <- merge(extract1, extract2, by.x = "Meth_Site_ID", by.y = "Methylation_ID")
    m_e3_e4 <- extract3
    m_e_all <- merge(m_e1_e2, m_e3_e4, by.x = "Source_ID", by.y = "DataSet_ID")
    m_e_all[, c("Meth_Site_ID", "Methylation_ID", "Source_ID", "DataSet_ID", "Data_ID")] <- NULL
    m_e_all
  })
})

# Run the application 
shinyApp(ui = ui, server = server)

