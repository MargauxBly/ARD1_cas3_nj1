library(shiny)
library(shinythemes)
library(plotly)
library(gridExtra)


shinyUI(
    navbarPage(title = "MODELE ARD 1",
               #theme = shinytheme("cerulean"), #
               # Application title
               tabsetPanel( 
                   tabPanel( title = "Wiener périodique",
                             
                             
                             sidebarLayout(
                                 sidebarPanel(
                                     
                                     selectInput("but", "Veuillez choisir le modele", 
                                                 choices = c("CAS 1", "CAS 2","CAS 3","CAS 4")),
                                     
                                     
                                     withMathJax(),
                                     fluidPage(fluidRow(column(4,
                                     numericInput("mu_x", label = helpText(style="color:black","$$\\mu \\:$$"),
                                                  value = 2)),
                                     
                                    
                                     column(4,
                                     withMathJax(),
                                     numericInput("var1", label = helpText(style="color:black","$$\\sigma^2 \\:$$"),
                                                  value = 5, min=0.0000001)),
                                     column(4,
                                     withMathJax(),
                                     numericInput("rho", label = helpText(style="color:black","$$\\rho $$"),
                                                  value = 0.5, min =0, max = 1)))),
                                     
                                     conditionalPanel(
                                         condition = "input.but != 'Qualité des estimateurs'"
                                         ,tags$hr()
                                         ,numericInput("n", label = "Taille de l'échantillon :"
                                                       ,value = 50, min=1)
                                         
                                     ),
                                     
                           
                                     
                                     numericInput("p", label = "Période:",
                                                  value = 1, min =1),
                                     numericInput("pm_periode", label = "Période de la maintenance :",
                                                  value = 12, min =1),
                          
                                     
                                     numericInput("nbtest",label="Nombre de trajectoires",
                                                  value=10),
                                     numericInput("num",label="Numero de la trajectoire noire",
                                                  value=5),
                                     
                                   
                                     actionButton("button", "Resimuler avec les mêmes paramètres"),
                                     
                                     fluidRow(column(3, verbatimTextOutput("")))
                                     
                                 ),
                                 
                                 
                                 
                                 mainPanel(
                                     tabsetPanel(type = "tabs",
                                                 tabPanel("Simulations"
                                                          ,textOutput (outputId = "Text1")
                                                          , plotlyOutput("mytable")
                                                          ,tableOutput("tabes")
                                                         ,plotOutput("bxARD1")
                                                          
                                                        
                                                          

                                                          
                                                 )
                                                 
                                                 
                                                 ,tabPanel("Performance des estimateurs"
                                                 ,plotOutput("perf")
                                     )
                                                 
                                                 
                                                 
                                     )
                                 )
                             )
                   ),
                   
                   
                   
                   tabPanel( title = "Wiener non forcement periodique",
                             
                             
                             sidebarLayout(
                                 sidebarPanel(
                                     
                                     selectInput("but2", "Veuillez choisir le modele", 
                                                 choices = c("CAS 1", "CAS 2","CAS 3","CAS 4")),
                                     
                                     
                                     fluidPage(fluidRow(column(4,
                                     withMathJax(),
                                     numericInput("mu_x2", label =helpText(style="color:black", "$$\\mu \\:$$"),
                                                  value = 2)),
                                     
                                     
                                     
                                     column(4,
                                     withMathJax(),
                                     numericInput("var2", label = helpText(style="color:black","$$\\sigma^2 \\:$$"),
                                                  value = 5, min=0.0000001)),
                                     column(4,
                                     withMathJax(),
                                     numericInput("rho2", label = helpText(style="color:black","$$\\rho \\:$$"),
                                                  value = 0.5, min =0, max = 1))))
                                  
                                     ,numericInput("n2", label = "Taille de l'échantillon :"
                                                   ,value = 50, min=1
                                     ),
                                     numericInput("pm_periode2", label = "Période de la maintenance :",
                                                  value = 12, min =10),
                                    
                                     
                                     numericInput("nbtest2",label="Nombre de trajectoires",
                                                  value=10, min=1),
                                     numericInput("num2",label="Numero de la trajectoire noire",
                                                  value=10),
                                     
                                     
                                     
                                     
                                  
                                     actionButton("bu2", "$$\\text{Resimuler avec les mêmes } \\Delta \\ \\text{T} $$"),
                                     
                                    
                                     actionButton("button2", "$$\\text{Resimuler avec des } \\Delta \\text{T differents}$$"),
                                     
                                     fluidRow(column(3, verbatimTextOutput("")))
                                     
                                 ),
                                 
                                 
                                 
                                 mainPanel(
                                     tabsetPanel(type = "tabs",
                                                 tabPanel("Simulation"
                                                          ,textOutput (outputId = "Text3")
                                                          , plotlyOutput("mytable2")
                                                          ,tableOutput("tabes2")
                                                          ,plotOutput("bxARD2")
                                                          ,plotOutput("beq2"))
                                                          
                                                 ,tabPanel("Performance des estimateurs"
                                                           ,plotOutput("perf2")
                                                 )
                                                          )
                                     )
                                 )
                             )
                  
                   
                   #  ,tabPanel(title="Comparaison des cas",
                   #           mainPanel(
                   #             tabsetPanel(type="tabs",
                   #                         tabPanel("Biais",
                   #             textOutput("texim"),
                   #             fluidPage(fluidRow(column(4,
                   #                                       
                   #             
                   #            imageOutput("image1")),column(4,imageOutput("image2"))
                   #             ,column(4,imageOutput("image3"))),
                   #             fluidRow(column(4,imageOutput("image4")),
                   #                      column(4,imageOutput("image5")),
                   #                       column(4,imageOutput("image6"))),
                   #             fluidRow(column(4,imageOutput("image7")),
                   #                     column(4,imageOutput("image8")),
                   #                      column(4,imageOutput("image9"))))
                   # 
                   #           ),
                   # tabPanel("EQM",
                   #          textOutput("texim2"),
                   #          fluidPage(fluidRow(column(4,imageOutput("image10")),
                   #                   column(4,imageOutput("image11")),
                   #                   column(4,imageOutput("image12"))),
                   #          fluidRow(column(4,imageOutput("image13")),
                   #                   column(4,imageOutput("image14")),
                   #                   column(4,imageOutput("image15"))),
                   #          fluidRow(column(4,imageOutput("image16")),
                   #                   column(4,imageOutput("image17")),
                   #                   column(4,imageOutput("image18"))))
                   #    )))
                   #  )
 
                                     
                   )
               ))