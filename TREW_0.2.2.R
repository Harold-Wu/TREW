################################################################################
# TREW_0.2.2_2016.8.1
# Database version: TREW_0.1.3
# fluidPage
# Input_Gene name, species, target, modification, content type.
# Output_Table with select and download function.
################################################################################

library(shiny)
library(sqldf)
library(DT)

#### Load the database.
db <- dbConnect(SQLite(), dbname= "TREW_0.1.3.db")
Genome_Location <- dbReadTable(db, "Genome_Location")
Sites_Info <- dbReadTable(db, "Sites_Info")
Source_Info <- dbReadTable(db, "Source_Info")
Gene_C <- unique(tolower(Sites_Info$Gene_ID))
#### Generate the indexes.
idx_2to1 = table(Genome_Location$Meth_Site_ID)
idx_3to2 = table(Sites_Info$Source_ID)[Source_Info$DataSet_ID]
idx_3to1 = table(rep(Sites_Info$Source_ID,idx_2to1))[Source_Info$DataSet_ID]
#### Load the functions.
#==#1.Search and generate full table (With the Record_quantity column).
GenerateFull <- function(Species = NULL, 
                         Target = NULL, 
                         Modification = NULL, 
                         Gene_ID = "CDK1",
                         Liftover = NULL)
{# Join the tables through indexes.
idx2 = !vector(length = dim(Sites_Info)[1])
idx3 = !vector(length = dim(Source_Info)[1])
# Option-Species.
if (!is.null(Species) & Species != "All") {
  idx3 <- Source_Info$Species == Species
}else{}
# Option-Regulation type.
if (Target == "Regulator") {
  idx3 <- idx3 & (Source_Info$Target != "YTHDF1" & Source_Info$Target != "YTHDF2")
}else{
  if (Target == "Regulatee") {
    idx3 <- idx3 & (Source_Info$Target == "YTHDF1" | Source_Info$Target == "YTHDF2")
  }else{}
}
# Option-Modification type.
if (!is.null(Modification) & Modification != "All") {
  idx3 <- idx3 & (Source_Info$Modification == Modification)
}else{}
# Option-Gene_ID.
if (!is.null(Gene_ID)) {
  idx2 <- (tolower(Sites_Info$Gene_ID) == tolower(Gene_ID)) & rep(idx3,idx_3to2)
}else{}

idx2[which(is.na(idx2))] <- FALSE

idx1 <- rep(idx2,idx_2to1)

Full_table <- cbind(Genome_Location[which(idx1),],
                    Sites_Info[rep(which(idx2),idx_2to1[which(idx2)]),],
                    Source_Info[rep(1:dim(Source_Info)[1],idx_3to1)[which(idx1)],])
# Option-Liftover.
if (Liftover == "No") {
  Full_table <- Full_table[which(is.na(Full_table[,"LiftOver"])),]
}else{}
# Add column-Positivity.
Positivity <- vector(mode = "character", length = length(Full_table[,"Target"]))

for (i in 1:length(Full_table[,"Target"])) {
  if (Full_table[i, "Diff_fdr"] < 0.05) {
    Positivity[i] <- "Positive"
  } else {
    Positivity[i] <- "Negative"
  }
}
# Add column-Record_Q.
Full_table <- cbind(Full_table, Positivity)

Record_Q <- vector(mode = "integer", length = length(Full_table[,"Target"]))

for (i in 1:length(Full_table[,"Target"])) {
  Record_Q[i] <- length(which(Full_table[,"Target"] == Full_table[,"Target"][i] & 
                                       Full_table[,"Modification"] == Full_table[,"Modification"][i]))
}
# Add column-Consistency_Q.
Consistent_Q <- vector(mode = "integer", length = length(Full_table[,"Target"]))

for (i in 1:length(Full_table[,"Target"])) {
  Consistent_Q[i] <- length(which(Full_table[,"Target"] == Full_table[,"Target"][i] & 
                                           Full_table[,"Modification"] == Full_table[,"Modification"][i] & 
                                           Full_table[,"Consistency"] == 1))
}
# Add column-Positive_Q.
Positive_Q <- vector(mode = "integer", length = length(Full_table[,"Target"]))

for (i in 1:length(Full_table[,"Target"])) {
  Positive_Q[i] <- length(which(Full_table[,"Target"] == Full_table[,"Target"][i] & 
                                  Full_table[,"Modification"] == Full_table[,"Modification"][i] & 
                                  Full_table[,"Positivity"] == "Positive"))
}

Full_table <- cbind(Full_table, Positive_Q, Consistent_Q, Record_Q)
# Add column-Negative_Q.
Negative_Q <- vector(mode = "integer", length = length(Full_table[,"Target"]))

for (i in 1:length(Full_table[,"Target"])) {
  Negative_Q[i] <- Full_table[i,"Record_Q"] - Full_table[i,"Positive_Q"]
}

Full_table <- cbind(Full_table, Negative_Q)
# Add column-Positivity_percent.
Positivity_percent <- vector(mode = "numeric", length = length(Full_table[,"Target"]))

for (i in 1:length(Full_table[,"Target"])) {
  
  Positivity_P <- Full_table[i,"Positive_Q"] / Full_table[i,"Record_Q"]
  
  Positivity_percent[i] <- paste(round(Positivity_P*100, digits = 2), "%", sep = '')
}

Full_table <- cbind(Full_table, Positivity_percent)
# Remove useless columns.
Full_table[, c("Note_t1", "Note_t2", "Note_t3")] <- NULL
# Assign row names.
rownames(Full_table) <- 1:length(Full_table[, "Target"])
Full_table
}

#==#2.Generate the general DT table with Record_quantity marked.
GenerateGeneral <- function(Species2 = NULL, 
                            Target2 = NULL, 
                            Modification2 = NULL, 
                            Gene_ID2 = "CDK1",
                            Liftover2 = NULL)
{# Generate full table.
FullTable <- GenerateFull(Species = Species2, 
                          Target = Target2, 
                          Modification = Modification2, 
                          Gene_ID = Gene_ID2,
                          Liftover = Liftover2)
# Extract needed columns.
FullGeneral <- FullTable[, c("Target", "Target_type", "Modification", "Record_Q", "Positivity_percent", 
                             "Positive_Q", "Negative_Q", "Consistent_Q")]
# Generate unique general table.
UniqueGeneral <- unique(FullGeneral)
}

#==#3.Generate the specific tables (Default & Completed) which were selected.
GenerateSpecific <- function(Species2 = NULL, 
                             Target2 = NULL, 
                             Modification2 = NULL, 
                             Gene_ID2 = "CDK1",
                             Liftover2 = NULL,
                             Content = "Default",
                             #Consistency_Positivity = NULL,
                             GeneralTable = NULL,
                             Select_Number = NULL)
{# Generate full table.
FullTable <- GenerateFull(Species = Species2, 
                          Target = Target2, 
                          Modification = Modification2, 
                          Gene_ID = Gene_ID2,
                          Liftover = Liftover2)

Select_Target <- GeneralTable[Select_Number, "Target"]
Select_Modification <- GeneralTable[Select_Number, "Modification"]

FullSpecific <- FullTable[which(FullTable[,"Target"] %in% Select_Target & 
                                  FullTable[,"Modification"] %in% Select_Modification),]

# if (Consistency_Positivity == "Only show positive and consistent results") {
#   FullSpecific <- FullSpecific[which(FullSpecific[, "Consistency"] == 1 & 
#                                        FullSpecific[, "Positivity"] == "+"),]
# } else {}

if (Content == "Default") {
  FullSpecific <- FullSpecific[, c("Source_ID", "Positivity", "Consistency", "Genome_assembly", 
                                   "Modification", "Technique", "Target", "Cell_line")]
  DT::datatable(FullSpecific, rownames= FALSE, filter = "top",
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
  DT::datatable(FullSpecific, rownames= FALSE, filter = "top",
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
  #theme = "lowdown.css",
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
                          column(width = 6, offset = 3,
                                 h2("TREW Database Search")),
                          column(width = 6, offset = 3,
                                 selectInput("Gene_choice", "", Gene_C, width = 600))
                        ),
                        br(),
                        br(),
                        br(),
                        fluidRow(
                          column(width = 8,
                                 DT::dataTableOutput("G_Table")),
                          column(width = 2,
                                 radioButtons("Species_choice", "Species", c("All", "Homo sapiens","Drosophila melanogaster", 
                                                                             "Mus musculus")),
                                 radioButtons("Regulation_choice", "Regulation type", c("Both", "Regulator", "Regulatee"))),
                          column(width = 2,
                                 radioButtons("Modification_choice", "Modification type", c("All", "m6A","m1A", "m5C", "Psi")),
                                 radioButtons("Liftover_choice", "Whether present liftover", c("Yes", "No"))),
                          br(),
                          br(),
                          column(width = 4,
                                 helpText("☟Please indicate the content of the specific table.")),
                          column(width = 4,
                                 selectInput("Content_choice", "Output table content", c("Default", "Completed")))
                        ),
                        hr(),
                        br(),
                        br(),
                        fluidRow(
                          column(width = 12,
                                 DT::dataTableOutput("S_Table"))

                        )),
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
  DT_General_Table <<- reactive({DT::datatable(General_Table(), rownames= FALSE, 
                                               extensions = list("ColReorder" = NULL,
                                                                 "Buttons" = NULL,
                                                                 "FixedHeader" = NULL,
                                                                 "Scroller" = NULL),
                                               options = list(
                                                 scrollX = TRUE,
                                                 deferRender = TRUE,
                                                 scrollY = 340,
                                                 scroller = TRUE,
                                                 dom = 'BRrlftpi',
                                                 autoWidth=TRUE,
                                                 fixedHeader = TRUE,
                                                 ColReorder = TRUE,
                                                 buttons =
                                                   list(
                                                     I('colvis')
                                                   )
                                               ))})
  output$G_Table <- DT::renderDataTable(DT_General_Table(), server = TRUE)
  #Select the Records and generate the specific tables (Default & Competed).
  Select_number <<- reactive({as.numeric(input$G_Table_rows_selected)})
  
  Specific_Table <- reactive({GenerateSpecific(Species = input$Species_choice,
                                               Target = input$Regulation_choice,
                                               Modification = input$Modification_choice,
                                               Gene_ID = input$Gene_choice,
                                               Liftover = input$Liftover_choice,
                                               Content = input$Content_choice,
                                               #Consistency_Positivity = input$Consistency_Positivity,
                                               GeneralTable = General_Table(),
                                               Select_Number = Select_number())})
  
  output$S_Table <- DT::renderDataTable(Specific_Table(), server = TRUE)
})


#### Run the application 
shinyApp(ui = ui, server = server)