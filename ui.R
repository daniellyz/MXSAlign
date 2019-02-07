library(shiny)
library(V8)
library(shinyjs)
library(d3heatmap)
require(DT, quietly = TRUE) 


shinyUI(navbarPage("MXSAlign: Data alignment from MaXis",
       tabPanel("A) Start a run",
          shinyjs::useShinyjs(),
          shinyjs::extendShinyjs(text = "shinyjs.refresh = function() { location.reload(); }"),         
          
          column(7,         
                 br(),
                 h3("Upload molecular feature files exported from Data Analysis 4.3"),
                 em('File header: #, Max. MW, RT [min], Area, Max. m/z'),
                 fileInput('file1',label='',multiple = TRUE,accept = c('.csv')),
                 
                 h3("Scan mode:"),
                 selectInput('Scan', label='',c('Positive ion mode','Negative ion mode')),

                 h3("Choose the mass variable you want to align:"),
                 selectInput('Mass', label='',c('Max. m/z','Molecular weight')),

                 h3("The start of LC run (min):"),
                 numericInput('Calibrant', label='',0.3,min=0,max=1,step=0.1),
                 h3("The end of LC run (min):"),
                 numericInput('End', label='',10,min=5,max=40,step=0.5),
                 h3("Choose the alignment order:"),
                 selectInput('Order', label='',c('Mass first','Retention time first')),
                 br(),
                 em('Messages from the server:'),
                 br(),
                 textOutput("blank_message1")
                ),
          
         column(5,
                h3("Choose the alignment mode:"),
                selectInput('Mode', label='',c('Fixed tolerance window','Parameter Tuning (Time-consuming)')),
                uiOutput("Controls"),
                checkboxInput("Fixed", label = "Second tolerance window for the RT range:", value = FALSE),
                uiOutput("Control_fix"),
                br(),
                tags$head(
                    tags$style(HTML('#goButton{background-color:lightgreen}'))
                  ),
                actionButton("goButton", "Start the run",style='padding:6px; font-size:150%'),
                br(),
                        
                br(),
                    tags$head(
                     tags$style(HTML('#killButton{background-color:orange}'))),
                actionButton("killButton", "Restart everything",style='padding:6px; font-size:150%')
                
                )),
       
       tabPanel("B) Tuning output",
           br(),
           h3("Number of molecular features by tuning mass tolerance window (row) and retention time window (column)"),
           br(),
           d3heatmapOutput("heat",width = "80%", height = "500px"),
           br(),
           uiOutput("Tune_choose"),
           br(),
           plotOutput("Tune_plot",height=1000, width = 1200)
       ),
      
       tabPanel("C) Alignment & Quality",
          br(),
          textOutput("blank_message2"),
          br(),
          downloadButton("downloadMatrix", "Download aligned matrix",style='padding:6px; font-size:150%'),
          br(),      
          br(),
          navbarPage("",
            tabPanel("General evaluation",
              column(6,   
              h3('Mass deviation distribution:'),
              em('Maximal deviation from the mass average'),
              plotOutput("Mdev",height =450, width = 600),
              h3('Retention time deviation distribution:'),
              em('Maximal deviation from the RT average'),
              plotOutput("Rdev",height =450, width = 600),
              h3('Missing value distribution:'),
              plotOutput("Zeros",height =450, width = 600)),
          
              column(4,
              br(),
              br(),
              dataTableOutput("table1"))),
            
            tabPanel("Mass deviation",
              br(),
              uiOutput("Sample_start1"),
              br(),
              plotOutput("plot_mass_sample",height = 500, width = 1500)
            ),
            
            tabPanel("Retention time deviation",
              br(),
              uiOutput("Sample_start2"),
              br(),
              plotOutput("plot_RT_sample",height = 500, width = 1500)))
            
            ),
       
       tabPanel("D) About ",
         includeMarkdown("About.Rmd"))
  
))
                   

