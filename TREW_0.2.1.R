################################################################################
# TREW_0.2.1_2016.7.26
# Database version: TREW_1.2
# fluidPage
# Input_Gene name, species, target, modification, content type.
# Output_Table with select and download function.
################################################################################

library(shiny)
library(sqldf)
library(DT)

#### Load the database.
db <- dbConnect(SQLite(), dbname= "TREW_1.2.db")
Genome_Location <- dbReadTable(db, "Genome_Location")
Sites_Info <- dbReadTable(db, "Sites_Info")
Source_Info <- dbReadTable(db, "Source_Info")
Gene_C <- unique(tolower(Sites_Info$Gene_ID))
#### Generate the indexes.
idx_2to1 = table(Genome_Location$Meth_Site_ID)
idx_3to2 = table(Sites_Info$Source_ID)[Source_Info$DataSet_ID]
idx_3to1 = table(rep(Sites_Info$Source_ID,idx_2to1))[Source_Info$DataSet_ID]
#### Load the functions.
#1.Search and generate full table (With the Record_quantity column).
GenerateFull <- function(Species = NULL, 
                         Target = NULL, 
                         Modification = NULL, 
                         Gene_ID = "CDK1",
                         Liftover = NULL)
{ idx2 = !vector(length = dim(Sites_Info)[1])
idx3 = !vector(length = dim(Source_Info)[1])
if (!is.null(Species) & Species != "All of three") {
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
  idx2 <- (tolower(Sites_Info$Gene_ID) == tolower(Gene_ID)) & rep(idx3,idx_3to2)
}else{}

idx2[which(is.na(idx2))] <- FALSE

idx1 <- rep(idx2,idx_2to1)

Full_table <- cbind(Genome_Location[which(idx1),],
                    Sites_Info[rep(which(idx2),idx_2to1[which(idx2)]),],
                    Source_Info[rep(1:dim(Source_Info)[1],idx_3to1)[which(idx1)],])

if (Liftover == "No") {
  Full_table <- Full_table[which(is.na(Full_table[,"LiftOver"])),]
}else{}

Positivity <- vector(mode = "character", length = length(Full_table[,"Target"]))

for (i in 1:length(Full_table[,"Target"])) {
  if (Full_table[i, "Diff_fdr"] < 0.05) {
    Positivity[i] <- "+"
  } else {
    Positivity[i] <- "-"
  }
}

Full_table <- cbind(Full_table, Positivity)

Record_quantity <- vector(mode = "integer", length = length(Full_table[,"Target"]))

for (i in 1:length(Full_table[,"Target"])) {
  Record_quantity[i] <- length(which(Full_table[,"Target"] == Full_table[,"Target"][i] & 
                                     Full_table[,"Modification"] == Full_table[,"Modification"][i]))
}

Positive_Q <- vector(mode = "integer", length = length(Full_table[,"Target"]))

for (i in 1:length(Full_table[,"Target"])) {
  Positive_Q[i] <- length(which(Full_table[,"Target"] == Full_table[,"Target"][i] & 
                                Full_table[,"Modification"] == Full_table[,"Modification"][i] & 
                                Full_table[,"Positivity"] == Full_table[,"Positivity"][i]))
}

Full_table <- cbind(Full_table, Positive_Q, Record_quantity)

Positivity_Percent <- vector(mode = "numeric", length = length(Full_table[,"Target"]))

for (i in 1:length(Full_table[,"Target"])) {
  if (Full_table[i,"Positivity"] == "+") {
    Positivity_P <- Full_table[i,"Positive_Q"] / Full_table[i,"Record_quantity"]
  } else {
    Positivity_P <- 1 - (Full_table[i,"Positive_Q"] / Full_table[i,"Record_quantity"])
  }
  Positivity_Percent[i] <- paste(round(Positivity_P*100, digits = 2), "%", sep = '')
}

Full_table <- cbind(Full_table, Positivity_Percent)

Full_table[, c("Note_t1", "Note_t2", "Note_t3")] <- NULL

rownames(Full_table) <- 1:length(Full_table[, "Target"])
Full_table
}

#1.Verify whether Gene name exists in the Database.
# VerifyGene <- function(Gene_ID = NULL)
# {if (sum(tolower(Sites_Info$Gene_ID) == tolower(Gene_ID)) == 0) {
#   No_Gene <- "The gene name could not be found, please check it."
# } else {
#   No_Gene <- NULL
# }
# 
# }

#2.Generate the general DT table with Record_quantity marked.
GenerateGeneral <- function(Species2 = NULL, 
                            Target2 = NULL, 
                            Modification2 = NULL, 
                            Gene_ID2 = "CDK1",
                            Liftover2 = NULL)
{ FullTable <- GenerateFull(Species = Species2, 
                            Target = Target2, 
                            Modification = Modification2, 
                            Gene_ID = Gene_ID2,
                            Liftover = Liftover2)

  FullGeneral <- FullTable[, c("Target", "Modification", "Positivity_Percent", "Record_quantity")]
  UniqueGeneral <- unique(FullGeneral)
  #DT::datatable(UniqueGeneral,options = list(lengthMenu = list(c(5, 50, -1), c('5', '10', 'All'))))
}

#3.Generate the specific tables (Default & Completed) which were selected.
GenerateSpecific <- function(Species2 = NULL, 
                            Target2 = NULL, 
                            Modification2 = NULL, 
                            Gene_ID2 = "CDK1",
                            Liftover2 = NULL,
                            Content = "Default",
                            Consistency_Positivity = NULL,
                            GeneralTable = NULL,
                            Select_Number = NULL)
{ FullTable <- GenerateFull(Species = Species2, 
                            Target = Target2, 
                            Modification = Modification2, 
                            Gene_ID = Gene_ID2,
                            Liftover = Liftover2)

  Select_Target <- GeneralTable[Select_Number, "Target"]
  Select_Modification <- GeneralTable[Select_Number, "Modification"]
  
  FullSpecific <- FullTable[which(FullTable[,"Target"] %in% Select_Target & 
                                  FullTable[,"Modification"] %in% Select_Modification),]
  
  if (Consistency_Positivity == "Only show positive and consistent results") {
    FullSpecific <- FullSpecific[which(FullSpecific[, "Consistency"] == 1 & 
                                       FullSpecific[, "Positivity"] == "+"),]
  } else {}

  if (Content == "Default") {
    FullSpecific <- FullSpecific[, c("Positivity" ,"Meth_Site_ID", "Source_ID", 
                                     "Consistency", "Genome_assembly", 
                                     "Modification", "Technique", "Target", "Cell_line")]
    DT::datatable(FullSpecific, rownames= FALSE,
                  extensions = list("ColReorder" = NULL,
                                    "Buttons" = NULL,
                                    "FixedHeader" = NULL,
                                    #"Responsive" = NULL,
                                    "Scroller" = NULL),
                  options = list(
                    scrollX = TRUE,
                    deferRender = TRUE,
                    scrollY = 200,
                    scroller = TRUE,
                    dom = 'BRrlftpi',
                    autoWidth=TRUE,
                    fixedHeader = TRUE,
                    lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                    ColReorder = TRUE,
                    buttons =
                      list(
                        'copy',
                        'print',
                        list(
                          extend = 'collection',
                          buttons = c('csv', 'excel'),
                          text = 'Download'
                        ),
                        I('colvis')
                      )
                    ))
  } else {
  DT::datatable(FullSpecific, rownames= FALSE,
                extensions = list("ColReorder" = NULL,
                                  "Buttons" = NULL,
                                  "FixedHeader" = NULL,
                                  #"Responsive" = NULL,
                                  "Scroller" = NULL),
                options = list(
                  scrollX = TRUE,
                  deferRender = TRUE,
                  scrollY = 200,
                  scroller = TRUE,
                  dom = 'BRrlftpi',
                  autoWidth=TRUE,
                  fixedHeader = TRUE,
                  lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                  ColReorder = TRUE,
                  buttons =
                    list(
                      'copy',
                      'print',
                      list(
                        extend = 'collection',
                        buttons = c('csv', 'excel'),
                        text = 'Download'
                      ),
                      I('colvis')
                    ),
                  columnDefs = list(list(
                    targets = 26,
                    render = JS(
                      "function(data, type, row, meta) {",
                      "return type === 'display' && data.length > 15 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 15) + '...</span>' : data;",
                      "}")
                  ))))
  }
}


#############################################################################################
#### UI
ui <- shinyUI(fluidPage(
  theme = "lowdown.css",
  titlePanel("TREW"),
  navlistPanel("Home", widths = c(2, 10),
               tabPanel("Introduction",
                        h2("Welcome to TREW"),
                        h3("epitranscriptomic targets of RNA modification readers, erasers and writers"),
                        br(),
                        hr(),
                        "In eukaryotic cells, the control of mRNA translation and degradation 
                        is critical for managing the quantity and duration of gene expression. 
                        Global translation regulation is typically achieved by modulating both 
                        the activity of translation initiation factors and the availability of 
                        ribosomes. Several process, including cell growth, division, and differentiation, 
                        are mediated by RNA-binding proteins and small complementary RNAs (microRNAs 
                        and short interfering RNAs).",
                        br(),
                        br(),
                        "However, as indicated by most recent evidence, dynamic modifications of mRNA, 
                        including 5-methylcytidine, pseudouridine, and N6-methyladenosinde, have emerged 
                        as potential new mechanisms of post-transcriptional gene regulation. Three types 
                        of proteins, methyltransferases, demethylases and specific binding proteins, which 
                        are also called writers, erasers and readers, are involved in these mechanisms.",
                        br(),
                        br(),
                        "This database, TREW, focus on the epitranscriptomic targets of RNA modification 
                        readers, erasers and writers and it consists of specific details of the targets of 
                        RNA modification. Users could easily obtain the detailed targets information by 
                        using the Search function and get an intuitive impression with the Visualization function."),
               tabPanel("Instruction",
                        h3("Options"),
                        "nvwijnl"),
               "Function",
               tabPanel("Search",
                        fluidRow(
                          column(width = 4, 
                                 helpText("☟Please indicate the gene name, species, regulation type, 
                                          and modification type."),
                                 selectInput("Gene_choice", "Gene name", Gene_C),
                                 selectInput("Species_choice", "Species", c("All of three", "Homo sapiens","Drosophila melanogaster", 
                                                                            "Mus musculus")),
                                 selectInput("Regulation_choice", "Regulation type", 
                                             c("Regulator", "Regulatee", "Both of above")),
                                 selectInput("Modification_choice", "Modification type", 
                                             c("m6A","m1A", "m5C", "Psi", "All of above")),
                                 selectInput("Liftover_choice", "Whether present liftover", 
                                             c("Yes", "No"))
                                 ),
                          column(width = 8,
                                 helpText("☟Choose the row(s) you are interested to see the detail."),
                                 DT::dataTableOutput("G_Table")
                                 )
                          ),
                        hr(),
                        fluidRow(
                          column(width = 4,
                                 helpText("☞Please indicate the content of the specific table.")),
                          column(width = 4,
                                 selectInput("Consistency_Positivity", "Consistency & Positivity", 
                                             c("Only show positive and consistent results", "Show all"))),
                          column(width = 4,
                                 selectInput("Content_choice", "Output table content", 
                                             c("Default", "Completed")))
                        ),
                        hr(),
                        fluidRow(
                          column(width = 12,
                                 DT::dataTableOutput("S_Table"))
                                 
                        )
                        ),
               tabPanel("Visualization"),
               "Help",
               tabPanel("Help")
              )
))

###########################################################################################
#### Server
server <- shinyServer(function(input, output) {
#Generate and output the general table.
  General_Table <<- reactive({GenerateGeneral(Species = input$Species_choice,
                                             Target = input$Regulation_choice,
                                             Modification = input$Modification_choice,
                                             Gene_ID = input$Gene_choice,
                                             Liftover = input$Liftover_choice)
  })
  output$G_Table <- DT::renderDataTable(General_Table(), server = TRUE)
#Select the Records and generate the specific tables (Default & Competed).
  Select_number <<- reactive({as.numeric(input$G_Table_rows_selected)})
  
  Specific_Table <- reactive({GenerateSpecific(Species = input$Species_choice,
                                              Target = input$Regulation_choice,
                                              Modification = input$Modification_choice,
                                              Gene_ID = input$Gene_choice,
                                              Liftover = input$Liftover_choice,
                                              Content = input$Content_choice,
                                              Consistency_Positivity = input$Consistency_Positivity,
                                              GeneralTable = General_Table(),
                                              Select_Number = Select_number())})
  
  output$S_Table <- DT::renderDataTable(Specific_Table(), server = TRUE)
})


#### Run the application 
shinyApp(ui = ui, server = server)