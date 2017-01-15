# Shiny UI for split test analysis
# See documentation here:

numericInput2 <- function(inputId, label, value = "", pixel_width = 50, ...){
  # This allows us to more specifically modify the input box size (i.e., make it small)
  div(style="display:inline-block",
      tags$label(label, `for` = inputId),
      tags$input(id = inputId, type = "number", value = value, style=paste("width:", pixel_width, "px", sep = ""),...)
  )
 }

shinyUI(fluidPage(
  # Application title
  titlePanel("Multi-arm split test analysis"),
  h6("Questions? Check out the ", a("documentation.", href = "https://github.com/lumoslabs/doctor-fishman/blob/master/README.md")),

  sidebarLayout(

    # Sidebar for data input
    sidebarPanel(

      radioButtons("data_input_type", label = "What type of input?",
                  choices = c("Upload a CSV of conversion data" ="upload",
                              "Type in conversion data" = "typein",
                              "Upload a CSV of cash data" ="upload_cash",
                              "Upload a CSV of survey data (BETA)" = "upload_survey")),

      #upload conversion data file
      conditionalPanel(
        condition = 'input.data_input_type == "upload"',
        fileInput("testfile", label = h4("Upload Conversion Data:") , accept = 'Content-type: text/csv'),
        helpText("This csv should have one row per user or group of users and 3 columns: 1) name of the test variant for this user or group 2) 'trials', the total users in this group 3) total successes in this group"),
        selectInput("split_col", label="Column of split variants", choices = NULL, selected = NA, multiple = F),
        selectInput("trials_col", label="Trials column", choices = NULL, selected = NA),
        selectInput("successes_col", label="Successes column", choices = NULL, selected = NA)
        ),
      #manually type in conversion data
      conditionalPanel(
        condition = 'input.data_input_type == "typein"',
        numericInput2("num_groups", label = h4("How Many Groups In Your Split Test?"), value = 2, pixel_width = 100),
        helpText("Type in data separated by commas, like this: control, 130005, 1086 "),
        uiOutput("dynamic_input_fields")
        ),
      #upload cash file
      conditionalPanel(
        condition = 'input.data_input_type == "upload_cash"',
        fileInput("testfile_cash", label = h4("Upload Cash Data") , accept = 'Content-type: text/csv'),
        helpText("This csv should have one row for every user in the test with 3 columns: 1) unique identifier (e.g., user_id), 2) variant 3) cash collected"),
        br(),
        selectInput("split_col_cash", label="Column of split variants", choices = NULL, selected = NA, multiple = F),
        selectInput("user_id_col", label="Unique Identifier column (e.g., user_id or client_id)", choices = NULL, selected = NA),
        selectInput("cash_col", label="Cash Collected column", choices = NULL, selected = NA)
        ),
    #upload survey file
    conditionalPanel(
      condition = 'input.data_input_type == "upload_survey"',
      fileInput("testfile_survey", label = h4("Upload Survey Data") , accept = 'Content-type: text/csv'),
      helpText("This csv should have one row for every user in the test with 3 columns: 1) unique identifier (e.g., user_id), 2) variant 3) numerical survey response"),
      br(),
      selectInput("split_col_survey", label="Column of split variants", choices = NULL, selected = NA, multiple = F),
      selectInput("user_id_col_survey", label="Unique Identifier column (e.g., user_id or client_id)", choices = NULL, selected = NA),
      selectInput("survey_col", label="Survey Response column", choices = NULL, selected = NA)
      ),
    br(),
    actionButton("run_data_button", label = h4("Analyze Data"),
                 style = "background-color: #0066CC;
                 font-family: 'Lobster';
                 color: white  ")

      ),


    # results are rendered in the main panel
    mainPanel(
      #TODO refactor instructions, make universal object and pull it out so it is BEFORE conditionalPanel, make it reset when you change across types

      # for conversion rate analysis only (upload or typed in)
      conditionalPanel(
        condition = 'input.data_input_type == "typein" || input.data_input_type == "upload"',
        h3(textOutput("header2")), # "Subscription Rate"
        br(),
        tableOutput("analysis_table"),
        tableOutput("ties_table"),
        h3(textOutput("instructions")),
        column(width=6, class = "well",
               h3(textOutput("header3")), # "Probability Density"
               plotOutput(outputId = "group_plot",height = "300px",width="400px")
               ),
        column(width = 6, class = "well",
               h3(textOutput("header4")), # "Differences From Control"
               plotOutput(outputId = "differences_plot",height = "400px",width="400px"),
               uiOutput("diff_plot_note"), # % Diff is defined as  100*(variant - control)/control
               selectInput("baseline_variant1_input", label="Baseline variant to which others are compared", choices = c("control"), selected = "control", multiple = F)
               ),
        br(),
        numericInput2("danger_criteria1", label = "Adjust your worst-case scenario tolerance (%)", value = 5, class = "input-small"),
        numericInput2("min_effect_size1", label = "Adjust your min effect size of interest (%)", value = 0, class = "input-small"),
        selectInput("candidate_winner", label="Adjust your candidate winner", choices = NULL, selected = NULL, multiple = F),
        uiOutput("summary_text")
      ),

      # CASH DATA ANALYSIS
      conditionalPanel(
        condition = 'input.data_input_type == "upload_cash"',
        h3(textOutput("cash_header2")), # "Conversion Rate Data"
        tableOutput("analysis_table_cash"),
        h4(textOutput("cash_header3")), # "Probability Density"
        h3(textOutput("instructions_cash")), # "Instructions"
        plotOutput(outputId = "group_plot_cash",height = "400px",width="600px"),
        br(),
        h3(textOutput("cash_header4")), # "Cash Per Subscription"
        tableOutput("cash_per_sub_table"),
        plotOutput(outputId = "cash_per_sub_plot", height = "500px", width = "500px"),
        br(),
        h3(textOutput("cash_header5")), # "Cash Per Event"
        tableOutput("cash_per_event_table"),
        plotOutput(outputId = "cash_per_event_plot", height = "500px", width = "500px")
      ),

      # SURVEY ANALYSIS
      conditionalPanel(
        condition = 'input.data_input_type == "upload_survey"',
        h3(textOutput("survey_header2")), # "Cash Per Event"
        tableOutput("survey_per_event_table"),
        plotOutput(outputId = "survey_per_event_plot", height = "600px", width = "600px")
      )
    )

  )
)
)
