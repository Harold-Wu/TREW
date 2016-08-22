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
                                                 textInput("Gene_choice", "Gene name", width = 200, placeholder = "e.g CDK1"),
                                                 selectInput("Species_choice", "Species", 
                                                             c("Homo sapiens", "Drosophila melanogaster", "Mus musculus"), width = 200),
                                                 selectInput("Regulation_choice", "Regulation type", 
                                                             c("Regulator", "Regulatee", "Both of above"), width = 200),
                                                 selectInput("Modification_choice", "Modification type", 
                                                             c("m1A", "m5C", "m6A", "Psi", "All of above"), width = 200),
                                                 selectInput("Content_choice", "Output table content", c("Default", "Complete"), 
                                                             width = 200),
                                                 submitButton(text = "Search")
                                    ),
                                    mainPanel(DT::dataTableOutput("Final_Table"),
                                              textOutput("No_item"))
                                  )),
                         tabPanel("Help")
                         
                                  ))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  library(sqldf)
  library(DT)
  #Load the database.
  con <- dbConnect(RSQLite::SQLite(), dbname="TREW_1.2.db")
  
  
  output$Final_Table <- DT::renderDataTable({
    #Extract colomu.
    d1 <- dbReadTable(con, "Genome_Location")
    d2 <- dbReadTable(con, "Sites_Info")
    d2 <- d2[which(d2$Gene_ID == input$Gene_choice),]
    d3 <- dbReadTable(con, "Source_Info")
    d3 <- d3[which(d3$Species == input$Species_choice),]
    # d2 <- sqldf("SELECT * from Sites_Info where Gene_ID = 'input$Gene_choice'", connection = con)
    # d3 <- sqldf("SElECT * from Source_Info where Species = 'input$Species_choice'" ,connection = con)
    # d4 <- dbReadTable(con, "Raw_data_records")
    
    if (input$Regulation_choice == "Regulator") {
      d3 <- d3[which(d3$Target != "YTHDF1" & d3$Target != "YTHDF2"),]
    } else {
      if (input$Regulation_choice == "Regulatee") {
        d3 <- d3[which(d3$Target == "YTHDF1" | d3$Target == "YTHDF2"),]
      } else {
        d3 <- d3
      }
    }
    
    if (input$Modification_choice == "All of above") {
      d3 <- d3
    } else {
      d3 <- d3[which(d3$Modification == input$Modification_choice),]
    }
    
    if (input$Content_choice == "Default") {
      d1 <- d1[, c("Meth_Site_ID", "Range_start", "Range_width", "Strand", "Chromosome")]
      d2 <- d2[, c("Methylation_ID", "Source_ID", "Consistency")]
      d3 <- d3[, c("DataSet_ID", "Genome_assembly", "Technique", "Target", "Cell_line")]
      #Merge the extract table
      m1 <- merge(d1, d2, by.x = "Meth_Site_ID", by.y = "Methylation_ID")
      m2 <- merge(m1, d3, by.x = "Source_ID", by.y = "DataSet_ID")
      m2[, c("Meth_Site_ID", "Methylation_ID", "Source_ID", "DataSet_ID")] <- NULL
    } else {
      m1 <- merge(d1, d2, by.x = "Meth_Site_ID", by.y = "Methylation_ID")
      m2 <- merge(m1, d3, by.x = "Source_ID", by.y = "DataSet_ID")
      #m2[, c("Meth_Site_ID", "Methylation_ID", "Source_ID", "DataSet_ID")] <- NULL
    }
    if (length(row(m2)) > 0) {
      DT::datatable(m2, style = "bootstrap", selection = "single", options = 
                      list(scrollX = TRUE, pagingType = "full_numbers"))
    } else {
      output$No_item <- renderText({
        Empty <- paste("No item searched.")
        Empty
      })
    }
  })
})

# Run the application 
shinyApp(ui = ui, server = server)

