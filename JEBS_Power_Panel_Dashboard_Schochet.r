# Version 2 - September 2023

suppressWarnings(suppressPackageStartupMessages(library('shiny')))
suppressWarnings(suppressPackageStartupMessages(library('ggplot2')))
suppressWarnings(suppressPackageStartupMessages(library('ggtext')))
suppressWarnings(suppressPackageStartupMessages(library('shinydashboard')))
suppressWarnings(suppressPackageStartupMessages(library('kableExtra')))
suppressWarnings(suppressPackageStartupMessages(library('scales')))

#
# DEFINE THE UI
#

ui <- shinyUI(
  dashboardPage(skin='green',
    dashboardHeader(titleWidth = 329, title=HTML("<em><b><font size='3'>Power_Panel</font></b></em>")),
    dashboardSidebar(width=350,
      tags$style(HTML(".sidebar { height: 90vh; overflow-y: auto; }" )),
      tags$style(HTML(type='text/css', ".irs-min, .irs-max, .irs-grid-text, .irs-grid-pol { visibility: hidden !important;}")),
      sidebarMenu(

        fluidPage(
          tags$head(tags$style(HTML(".textInput {margin-left:5px}"))),    # This adds te left margin
        helpText(
          div(strong(HTML("</p>This dashboard conducts power calculations for </br> commonly used panel data
                          designs using the </br> methods developed in Schochet (Journal of </br>
                          Educational and Behavioral Statistics, 2022) </p> After entering the
                          inputs below, you can press </br> Submit at the bottom to view and save the results"))),
            style = "color: #27CF82;" )),  #27CF82;CAFF70

        #fileInput("file1", label = "Choose the RCT-YES .csv file to upload", accept = c(".csv")),

        fluidPage(
          tags$head(tags$style(HTML(".textInput {margin-left:5px}"))),
        helpText(
          div(strong(HTML("Key design information"))),
          style = "color:#33C5FF"
        )),

        # CREATE THE SIDEBAR PANEL INPUTS

        radioButtons(
          inputId = "did_cits",
          label   = "Design type",
          choices = c(
            "Difference-in-differences (DID)",
            "Comparative interrupted time series (CITS)",
            "Interrupted time series (ITS)"),
          inline = FALSE
        ),

        conditionalPanel(
          'input.did_cits === "Comparative interrupted time series (CITS)"
          | input.did_cits === "Interrupted time series (ITS)"',
          radioButtons(
            inputId = "type_cits",
            label   = "Linear model specification",
            choices = c(
              "Different pre- and post-period slopes",
              "Common pre- and post-period slopes",
              "Discrete post-period indicators (like DID)"),
            inline = FALSE
          )
        ),

        radioButtons(
          inputId = "samp_size",
          label   = "Type of power analysis",
          choices = c(
            "Calculate the minimum detectable effect size (MDE) for a given sample size",
            "Calculate the number of clusters required to achieve a target MDE value"),
          inline = FALSE
        ),

        conditionalPanel(
          'input.samp_size === "Calculate the number of clusters required to achieve a target MDE value"',
          numericInput(
            inputId = "mde",
            label = "Target MDE value",
            value = .20,
            min = .01,
            step = .01
          )
        ),

        radioButtons(
          inputId = "cross_sec",
          label   = "Type of panel data",
          choices = c(
            "Cross-sectional (different people over time)",
            "Longitudinal (same people over time)"),
          inline = FALSE
        ),

        fluidPage(
          tags$head(tags$style(HTML(".textInput {margin-left:5px}"))),
        helpText(
          div(strong(HTML("Time periods and intervals"))),
          style = "color:#33C5FF"
        )
        ),

        numericInput(
          inputId = "ntimep",
          label = "Number of time periods (integer)",
          min = 2,
          step = 1,
          value = 10
        ),

        radioButtons(
          inputId = "even_time",
          label   = "Spacing of time intervals",
          choices = c(
            "Evenly spaced",
            "Not evenly spaced"),
          inline = TRUE
        ),

        #uiOutput("ntboxes"),

        conditionalPanel(
          'input.even_time === "Not evenly spaced"',
          textInput(
            inputId = "time_intervals",
            label   = "Time interval for each time period, from earliest to most recent,
            separated by commas (integers)",
            value   = "1,5,8,9,12,15,19,20,22,23",
          )
        ),

        fluidPage(
          tags$head(tags$style(HTML(".textInput {margin-left:5px}"))),
          helpText(
            div(strong(HTML("Treatment timing groups"))),
            style = "color:#33C5FF"
          )
        ),

        numericInput(
          inputId = "ntimeg",
          label = "Number of treatment timing groups (integer)",
          min = 1,
          step = 1,
          value = 2
        ),

        textInput(
          inputId = "start_time",
          label   = "Start period for each timing group, separated by commas (integers from
          2 to the number of time periods)",
          value   = "4,6"
        ),

        radioButtons(
          inputId = "avg_point",
          label   = "Desired post-treatment period for the power analysis",
          choices = c(
            "Average across all post-periods",
            "Single post-period time point"),
          inline = FALSE
        ),

        conditionalPanel(
          'input.avg_point === "Single post-period time point"',
          radioButtons(
            inputId = "cal_expos",
            label   = "Single post-period measurement unit",
            choices = c(
              "At a specific time period",
              "At a specific treatment exposure point"),
            inline = FALSE
          )
        ),

        conditionalPanel(
          'input.avg_point === "Single post-period time point"
          & input.cal_expos === "At a specific time period"',
          numericInput(
            inputId = "q_time",
            label = "Specific time point (integer, e.g., where 6 is the 6th time point
            starting from the 1st pre-period, even with uneven intervals)",
            min = 1,
            step = 1,
            value = 6
          )
        ),

        conditionalPanel(
          'input.avg_point === "Single post-period time point"
          & input.cal_expos === "At a specific treatment exposure point"',
          numericInput(
            inputId = "l_time",
            label = "Specific exposure point (integer, e.g., where 6 is 6 time periods
            from the treatment start date, even with uneven intervals)",
            min = 1,
            step = 1,
            value = 6
          )
        ),

        fluidPage(
          tags$head(tags$style(HTML(".textInput {margin-left:5px}"))),
          helpText(
            div(strong(HTML("Sample sizes: Clusters are the unit of treatment </br> assignment,
                          such as a state, school, or hospital </br>
                          (see the Power_Panel Documentation.pdf file)"))),
            style = "color:#33C5FF"
          )),

        conditionalPanel(
          'input.samp_size == "Calculate the minimum detectable effect size (MDE) for a given sample size"',
          numericInput(
            inputId = "mt",
            label = "Total number of treatment clusters (integer)",
            min = 2,
            step = 1,
            value = 20),
          textInput(
            inputId = "mtk",
            label   = "Number of treatment clusters in each timing group, separated by commas (integers)",
            value   = "10,10")
        ),

        conditionalPanel(
          'input.samp_size == "Calculate the minimum detectable effect size (MDE) for a given sample size"
          & input.did_cits != "Interrupted time series (ITS)"',
          numericInput(
            inputId = "mc",
            label = "Total number of comparison clusters (integer)",
            min = 2,
            step = 1,
            value = 20),
          textInput(
            inputId = "mck",
            label   = "Number of matched comparison clusters in each timing group, separated by commas (integers)",
            value   = "10,10")
        ),

        conditionalPanel(
          'input.samp_size == "Calculate the number of clusters required to achieve a target MDE value"
          & input.did_cits != "Interrupted time series (ITS)"',
          numericInput(
            inputId = "rt",
            label = "Proportion of all clusters that are treatment clusters (0 to 1)",
            min = .01,
            max = .99,
            step = .01,
            value = .50)
        ),

        conditionalPanel(
          'input.samp_size == "Calculate the number of clusters required to achieve a target MDE value"',
          textInput(
            inputId = "rtk",
            label   = "Proportion of treatment clusters in each timing group,
            separated by commas (0 to 1; must sum to 1)",
            value   = ".50,.50")
        ),

        conditionalPanel(
          'input.samp_size == "Calculate the number of clusters required to achieve a target MDE value"
          & input.did_cits != "Interrupted time series (ITS)"',
          textInput(
            inputId = "rck",
            label   = "Proportion of matched comparison clusters in each timing group,
            separated by commas (0 to 1; must sum to 1)",
            value   = ".50,.50")
        ),

        numericInput(
        inputId = "nsamp",
        label = "Number of individuals per cluster per time period (integer)",
        min = 1,
        step = 1,
        value = 100
        ),

        fluidPage(
          tags$head(tags$style(HTML(".textInput {margin-left:5px}"))),
        helpText(
          div(strong(HTML("Error structure"))),
          style = "color:#33C5FF"
        )),

        radioButtons(
          inputId = "ar1",
          label   = "Autocorrelation structure of cluster-level errors over time",
          choices = c(
            "AR1",
            "Constant",
            "None"),
          inline = TRUE
        ),

        conditionalPanel(
          'input.ar1 != "None"',
          numericInput(
            inputId = "rho",
            label = "Autocorrelation parameter (0 to 1)",
            min = -.99,
            max = .99,
            step = .01,
            value = .40
          )
        ),

        conditionalPanel(
          'input.ar1 != "None" & input.cross_sec==="Longitudinal (same people over time)"',
          numericInput(
            inputId = "phi",
            label = "Autocorrelation parameter of individual-level errors (0 to 1)",
            min = -.99,
            max = .99,
            step = .01,
            value = .40
          )
        ),

        numericInput(
          inputId = "icc",
          label = "Clustering effect: Intraclass correlation coeff. (ICC) within
          cluster-time cells (0 to 1)",
          min = 0,
          max = 1,
          step = .01,
          value = .05
        ),

        fluidPage(
          tags$head(tags$style(HTML(".textInput {margin-left:5px}"))),
          helpText(
            div(strong(HTML("Weights"))),
            style = "color:#33C5FF"
          )
        ),

#        radioButtons(
#          inputId = "wtg",
#          label   = "Weights for the treatment timing groups for pooling to obtain aggregate effects",
#          choices = c(
#            "Equal weighting",
#            "By treatment exposure",
#            "By size",
#            "By exposure*size"
#          ),
#          selected = "By treatment exposure",
#          inline = FALSE),

        numericInput(
          inputId = "deff_wgt",
          label = "Design effect (variance inflation factor) due to various forms of weighting (at least 1)",
          min = 1,
          max = 20,
          step = .1,
          value = 1
        ),

        fluidPage(
          tags$head(tags$style(HTML(".textInput {margin-left:5px}"))),
        helpText(
          div(strong(HTML("Precision gains from model covariates"))),
          style = "color:#33C5FF"
        )),

        numericInput(
          inputId = "r2yx",
          label = "Regression R2 value from model covariates that improves precision (0 to 1)",
          min = 0,
          max = .99,
          step = .01,
          value = 0
        ),

        numericInput(
          inputId = "r2tx",
          label = "Correlation between model treatment indicators and covariates that reduces precision
          (0 to 1)",
          min = 0,
          max = .99,
          step = .01,
          value = 0
        ),

        fluidPage(
          tags$head(tags$style(HTML(".textInput {margin-left:5px}"))),
        helpText(
          div(strong(HTML("Hypothesis testing parameters"))),
          style = "color:#33C5FF"
        )),

        numericInput(
          inputId = "alpha",
          label = "Significance (alpha) level (0 to 1)",
          min = 0,
          max = .99,
          step = .01,
          value = 0.05
        ),

        radioButtons(
          inputId = "two_tailed",
          label   = "Two- or one-tailed test",
          choices = c(
            "Two-tailed",
            "One-tailed"),
          inline = TRUE
        ),

          sliderInput(
            inputId = "power",
            label = "Power level range: Calculations are conducted between min and max values",
            min = .50,
            max = .99,
            #step = 0.5,
            value = c(.60,.90)
          ),

        radioButtons(
          inputId = "power_by",
          label   = "Power level interval",
          choices = c(
            ".01",
            ".05",
            ".10"),
          selected = ".05",
          inline = TRUE
        ),

          div(actionButton("submit_button", "Submit"), class="form-group shiny-input-container")
       )
    ),

    # DEFINE THE MAIN PANEL ELEMENTS. A SplitLayout IS USED SO THAT THERE ARE
    # TWO PANELS FOR PRINTING THE TABLES AND GRAPHS

    dashboardBody(
     splitLayout(
       htmlOutput("tableset"),     # These are critical for printing - need html for kableExtra
       plotOutput("plotset"),      # This is for Panel 2 for the plots
       cellArgs = list(style='white-space: normal;')

       #tableOutput("table"),
       #plotOutput("plot"),
       #textOutput("vv"),
     ),

     hr(),
     textInput(
       inputId = "file_name",
       label   = "Name of file to save the graph (exclude the file path and extension)",
       value = "Power_Panel graph"),
     radioButtons(
       inputId = "png_pdf",
       label   = "Save graph as a png or pdf file",
       choices = c("png" = "png", "pdf" = "pdf"),
       inline = FALSE),
     downloadButton('down', 'Download Graph')
     )
    )
  )

server <- shinyServer(
  function(input, output, session) {

  # The following exits into R when exiting Shiny- this is nice!
    session$onSessionEnded(function() {
      stopApp(returnValue = invisible())
      })

    #
    # THIS FUNCTION CALCULATES THE RHO AND TIME TERMS FOR THE VARIANCES AND COVARIANCES
    #

    calc_rho <- function(corrp,prepost,pres,pref,posts,postf,tbarb1,tbarb2,ar,tpf) {

      if (prepost==1) {                 # PREPOST == 1

        rhog_prek <- 0
        btemp <- pref - pres + 1

        if (btemp>1) {

          vdf1 <- as.numeric(seq(1:btemp)) + pres  - 1
          vdf  <- expand.grid(vdf1,vdf1)
          vv   <- vdf[which(vdf[,2]>vdf[,1]),]         # Pull off observations where vdf2>vdf1

          if (ar==1) {
            sv <- (corrp)^(tpf[vv[,2]]-tpf[vv[,1]])  # Calculates sum of rhos
            rhog_prek <- sum(sv)
          } else if (ar==0)
          {
            rhog_prek <- corrp*nrow(vv)
          }
          btemp1 <- btemp-1
          rhog_prek <- 2*rhog_prek/(btemp*btemp1)
        }

        rho_ret <- rhog_prek
      } else if (prepost==2)                # PREPOST == 2
    {
      rhog_postk <- 0
      atemp <- postf - posts + 1
      if (atemp>1) {

        vdf1 <- as.numeric(seq(1:atemp)) + posts - 1
        vdf  <- expand.grid(vdf1,vdf1)
        vv   <- vdf[which(vdf[,2]>vdf[,1]),]

        if (ar==1) {
          sv <- (corrp)^(tpf[vv[,2]]-tpf[vv[,1]])  # Calculates sum of rhos
          rhog_postk <- sum(sv)
        } else if (ar==0)
        {
          rhog_postk <- corrp*nrow(vv)
        }
        atemp1 <- atemp-1
        rhog_postk <- 2*rhog_postk/(atemp*atemp1)
      }

      rho_ret <- rhog_postk

    } else if (prepost==3)                 # PREPOST == 3
    {
      rhog_ppk <- 0
      btemp <- pref - pres + 1
      atemp <- postf - posts + 1

      vdf1 <- as.numeric(seq(1:btemp)) + pres  - 1
      vdf2 <- as.numeric(seq(1:atemp)) + posts - 1
      vv   <- expand.grid(vdf1,vdf2)

      if (ar==1) {
        sv <- (corrp)^(tpf[vv[,2]]-tpf[vv[,1]])  # Calculates sum of rhos
        rhog_ppk <- sum(as.numeric(sv))
      } else if (ar==0)
      {
        rhog_ppk <- corrp*nrow(vv)
      }

      rho_ret <- rhog_ppk/(atemp*btemp)

    } else if (prepost==4)                  # PREPOST==4
    {
      post_diffk <- 0
      atemp = postf - posts + 1
      if (atemp>1) {

        vv <- as.numeric(seq(1:atemp)) + posts - 1
        sv <- (tpf[vv] - as.numeric(tbarb1))  # Calculates sum of time_star
        post_diffk <- sum(as.numeric(sv))
      }

      rho_ret <- post_diffk

    } else if (prepost==5)                 # PREPOST==5
    {
      post_covk <- 0
      atemp <- postf - posts + 1
      if (atemp>=1) {

        vdf1 <- as.numeric(seq(1:atemp)) + posts - 1
        vv <- expand.grid(vdf1,vdf1)
        sv <- (tpf[vv[,1]] - as.numeric(tbarb1))*(tpf[vv[,2]] - as.numeric(tbarb1))
        post_covk <- sum(as.numeric(sv))
      }

      rho_ret <- post_covk

    } else if (prepost==6)                 # PREPOST==6
    {
      rhog_pre1k <- 0
      btemp <- pref - pres + 1
      if (btemp>1) {

        vdf1 <- as.numeric(seq(1:btemp)) + pres - 1
        vdf  <- expand.grid(vdf1,vdf1)
        vv   <- vdf[which(vdf[,2]>vdf[,1]),]

        time_term <- (tpf[vv[,1]] - as.numeric(tbarb1))*(tpf[vv[,2]] - as.numeric(tbarb1))
        if (ar==1) {
          sv <- time_term*((corrp)^(tpf[vv[,2]]-tpf[vv[,1]]))
          rhog_pre1k <- sum(sv)
        } else if (ar==0)
        {
          sv <- time_term*(corrp)
          rhog_pre1k <- sum(sv)
        }
        btemp1 <- btemp-1
        rhog_pre1k <- 2*rhog_pre1k/(btemp*btemp1)
      }

      rho_ret <- rhog_pre1k

    } else if (prepost==7)                 # PREPOST==7
    {
      rhog_pre2k <- 0
      btemp <- pref - pres + 1

      if (btemp>1) {
        vdf1 <- as.numeric(seq(1:btemp)) + pres - 1
        vdf  <- expand.grid(vdf1,vdf1)
        vv   <- vdf #vdf[which(vdf[,2]!=vdf[,1]),]

        time_term <- (tpf[vv[,1]] - as.numeric(tbarb1))
        if (ar==1) {
          sv <- time_term*((corrp)^abs((tpf[vv[,2]]-tpf[vv[,1]])))
          rhog_pre2k <- sum(sv)
        } else if (ar==0)
        {
          sv <- time_term*(corrp)
          rhog_pre2k <- sum(sv)
        }
        btemp1 <- btemp-1
        rhog_pre2k <- rhog_pre2k/(btemp*btemp) #(btemp*btemp1)
      }

      rho_ret <- rhog_pre2k

    } else if (prepost==8)                 # PREPOST==8
    {
      rhog_pp1k <- 0
      btemp <- pref - pres + 1
      atemp <- postf - posts + 1

      vdf1 <- as.numeric(seq(1:btemp)) + pres - 1
      vdf2 <- as.numeric(seq(1:atemp)) + posts - 1
      vv   <- expand.grid(vdf1,vdf2)

      time_term <- (tpf[vv[,1]] - as.numeric(tbarb1))

      if (ar==1) {
        sv <- time_term*((corrp)^(tpf[vv[,2]]-tpf[vv[,1]]))
        rhog_pp1k <- sum(as.numeric(sv))
      } else if (ar==0)
      {
        sv <- time_term*(corrp)
        rhog_pp1k <- sum(as.numeric(sv))
      }

      rho_ret = rhog_pp1k/(btemp*atemp)

    } else if (prepost==9)                 # PREPOST==9
    {
      rhog_post1k <- 0
      atemp <- postf - posts + 1
      if (atemp>1) {

        vdf1 <- as.numeric(seq(1:atemp)) + posts - 1
        vdf  <- expand.grid(vdf1,vdf1)
        vv   <- vdf[which(vdf[,2]>vdf[,1]),]

        time_term <- (tpf[vv[,1]] - as.numeric(tbarb1))*(tpf[vv[,2]] - as.numeric(tbarb1))
        if (ar==1) {
          sv <- time_term*((corrp)^(tpf[vv[,2]]-tpf[vv[,1]]))
          rhog_post1k <- sum(sv)
        } else if (ar==0)
        {
          sv <- time_term*(corrp)
          rhog_post1k <- sum(sv)
        }
        atemp1 <- atemp-1
        rhog_post1k <- 2*rhog_post1k/(atemp*atemp1)
      }

      rho_ret <- rhog_post1k

    } else if (prepost==10)                 # PREPOST==10
    {
      rhog_post2k <- 0
      atemp <- postf - posts + 1

      if (atemp>1) {
        vdf1 <- as.numeric(seq(1:atemp)) + posts - 1
        vdf  <- expand.grid(vdf1,vdf1)
        vv   <- vdf #vdf[which(vdf[,2]!=vdf[,1]),]

        time_term <- (tpf[vv[,1]] - as.numeric(tbarb1))
        if (ar==1) {
          sv <- time_term*((corrp)^abs((tpf[vv[,2]]-tpf[vv[,1]])))
          rhog_post2k <- sum(sv)
        } else if (ar==0)
        {
          sv <- time_term*(corrp)
          rhog_post2k <- sum(sv)
        }
        atemp1 <- atemp-1
        rhog_post2k <- rhog_post2k/(atemp*atemp) #(atemp*atemp1)
      }

      rho_ret <- rhog_post2k

    } else if (prepost==11)                 # PREPOST==11
    {
      rhog_pp2k <- 0
      btemp <- pref - pres + 1
      atemp <- postf - posts + 1

      vdf1 <- as.numeric(seq(1:btemp)) + pres - 1
      vdf2 <- as.numeric(seq(1:atemp)) + posts - 1
      vv   <- expand.grid(vdf1,vdf2)

      time_term <- (tpf[vv[,2]] - as.numeric(tbarb1))

      if (ar==1) {
        sv <- time_term*((corrp)^(tpf[vv[,2]]-tpf[vv[,1]]))
        rhog_pp2k <- sum(as.numeric(sv))
      } else if (ar==0)
      {
        sv <- time_term*(corrp)
        rhog_pp2k <- sum(as.numeric(sv))
      }

      rho_ret = rhog_pp2k/(btemp*atemp)

    } else if (prepost==12)                 # PREPOST==12
    {
      rhog_pp4k <- 0
      btemp <- pref - pres + 1
      atemp <- postf - posts + 1

      vdf1 <- as.numeric(seq(1:btemp)) + pres - 1
      vdf2 <- as.numeric(seq(1:atemp)) + posts - 1
      vv   <- expand.grid(vdf1,vdf2)

      time_term <- (tpf[vv[,1]] - as.numeric(tbarb1))*(tpf[vv[,2]] - as.numeric(tbarb2))

      if (ar==1) {
        sv <- time_term*((corrp)^(tpf[vv[,2]]-tpf[vv[,1]]))
        rhog_pp4k <- sum(as.numeric(sv))
      } else if (ar==0)
      {
        sv <- time_term*(corrp)
        rhog_pp4k <- sum(as.numeric(sv))
      }

      rho_ret = rhog_pp4k/(btemp*atemp)

    } else if (prepost==90)                 # PREPOST==90
    {
      pre_diffk <- 0
      btemp <- pref - pres + 1

      if (btemp>1) {

        vv <- as.numeric(seq(1:btemp)) + pres - 1
        sv <- (tpf[vv] - as.numeric(tbarb1))*(tpf[vv] - as.numeric(tbarb2))
        pre_diffk <- sum(as.numeric(sv))
      }

      rho_ret <- pre_diffk

    } else if (prepost==100) {              # PREPOST==100

      rhog_pre3k <- 0
      btemp <- pref - pres + 1

      if (btemp>1) {
        vdf1 <- as.numeric(seq(1:btemp)) + pres - 1
        vdf  <- expand.grid(vdf1,vdf1)
        vv   <- vdf[which(vdf[,2]!=vdf[,1]),]

        time_term <- (tpf[vv[,1]] - as.numeric(tbarb1))*(tpf[vv[,2]] - as.numeric(tbarb2))
        if (ar==1) {
          sv <- time_term*((corrp)^abs((tpf[vv[,2]]-tpf[vv[,1]])))
          rhog_pre3k <- sum(sv)
        } else if (ar==0)
        {
          sv <- time_term*(corrp)
          rhog_pre3k <- sum(sv)
        }
        btemp1 <- btemp-1
        rhog_pre3k <- rhog_pre3k/(btemp*btemp1)
      }

      rho_ret <- rhog_pre3k

    } else if (prepost==110)                # PREPOST==110
    {
      rhog_pp2k <- 0
      btemp <- pref - pres + 1
      atemp <- postf - posts + 1

      vdf1 <- as.numeric(seq(1:btemp)) + pres - 1
      vdf2 <- as.numeric(seq(1:atemp)) + posts - 1
      vv   <- expand.grid(vdf1,vdf2)

      time_term <- (tpf[vv[,1]] - as.numeric(tbarb1))*(tpf[vv[,2]] - as.numeric(tbarb2))

      if (ar==1) {
        sv <- time_term*((corrp)^(tpf[vv[,2]]-tpf[vv[,1]]))
        rhog_pp2k <- sum(as.numeric(sv))
      } else if (ar==0)
      {
        sv <- time_term*(corrp)
        rhog_pp2k <- sum(as.numeric(sv))
      }

      rho_ret = rhog_pp2k/(btemp*atemp)
    }

      # End of function
      return(rho_ret)
    }

    #
    # THIS IS THE MAIN PROGRAM THAT DOES ALL THE CALCULATIONS AND RETURNS THE RESULTS
    # INTO POW_RES. IT IS REACTIVE SO UPDATES WHEN SIDE BAR INPUTS ARE CHANGED
    #

    pow_res <- reactive({

      # Read in inputs and save them to variables which makes it easier to refer to them

      if (input$did_cits == "Difference-in-differences (DID)") {
        did_cits <- 1
      } else if (input$did_cits == "Comparative interrupted time series (CITS)") {
        did_cits <- 2
      } else if (input$did_cits == "Interrupted time series (ITS)") {
        did_cits <- 3
      }

      if (input$type_cits == "Different pre- and post-period slopes") {
        type_cits <- 1
      } else if (input$type_cits == "Common pre- and post-period slopes") {
        type_cits <- 2
      } else if (input$type_cits == "Discrete post-period indicators (like DID)") {
        type_cits <- 3
      }

      if (input$samp_size == "Calculate the minimum detectable effect size (MDE) for a given sample size") {
        samp_size <- 1
      } else if (input$samp_size == "Calculate the number of clusters required to achieve a target MDE value") {
        samp_size <- 2
      }

      mde <- input$mde

      if (input$cross_sec=="Cross-sectional (different people over time)") {
        cross_sec <- 1
      } else if (input$cross_sec=="Longitudinal (same people over time)") {
        cross_sec <- 2
      }

      ntimep <- input$ntimep

      if (input$even_time=="Evenly spaced") {
        even_time <- 1
      } else if (input$even_time=="Not evenly spaced") {
        even_time <- 2
      }

      tpa <- input$time_intervals

      ntimeg <- input$ntimeg

      ska <- input$start_time

      if (input$avg_point == "Average across all post-periods") {
        avg_point <- 1
      } else if (input$avg_point == "Single post-period time point") {
        avg_point <- 2
      }

      if (input$cal_expos == "At a specific time period") {
        cal_expos <- 1
      } else if (input$cal_expos == "At a specific treatment exposure point") {
        cal_expos <- 2
      }

      q_time <- input$q_time
      l_time <- input$l_time

      wtg <- 2

      nsamp <- input$nsamp
      mt    <- input$mt
      mtka  <- input$mtk
      rt    <- input$rt
      rtka  <- input$rtk
      mc    <- input$mc
      mcka  <- input$mck
      rcka  <- input$rck

      if (input$ar1=="AR1") {
        ar1 <- 1
      } else if (input$ar1=="Constant") {
        ar1 <- 2
      } else if (input$ar1 == "None") {
        ar1 <- 3
      }

      rhoa <- input$rho
      phia <- input$phi
      icca <- input$icc
      deff_wgta <- input$deff_wgt

      alphaa <- input$alpha

      if (input$two_tailed=="Two-tailed") {
        two_tailed <- 1
      } else if (input$two_tailed=="One-tailed") {
        two_tailed <- 0
      }

      power <- input$power
      power_s <- power[1]
      power_f <- power[2]

      if (input$power_by == ".05") {
        power_by <- .05
      } else if (input$power_by == ".01") {
        power_by <- .01
      } else if (input$power_by == ".10") {
        power_by <- .10
      }

      r2yxa <- input$r2yx
      r2txa <- input$r2tx

      # Remove commas and other weird symbols from the TextInput numeric inputs

      tpb  <- gsub("\\D"," ",tpa)
      skb  <- gsub("\\D"," ",ska)
      rtkb <- gsub("[^0123456789.]"," ",rtka)
      mtkb <- gsub("\\D"," ",mtka)
      rckb <- gsub("[^0123456789.]"," ",rcka)
      mckb <- gsub("\\D"," ",mcka)

      alphab  <- gsub("[^0123456789.]"," ",alphaa)
      iccb    <- gsub("[^0123456789.]"," ",icca)
      rhob    <- gsub("[^0123456789.]"," ",rhoa)
      phib    <- gsub("[^0123456789.]"," ",phia)
      r2yxb   <- gsub("[^0123456789.]"," ",r2yxa)
      r2txb   <- gsub("[^0123456789.]"," ",r2txa)
      deff_wgtb  <- gsub("[^0123456789.]"," ",deff_wgta)

      # These steps are needed to make the text variables numeric

      tp  <- as.numeric(unlist(strsplit(tpb, ' ')))
      sk  <- as.numeric(unlist(strsplit(skb, ' ')))
      rtk <- as.numeric(unlist(strsplit(rtkb, ' ')))
      mtk <- as.numeric(unlist(strsplit(mtkb, ' ')))
      rck <- as.numeric(unlist(strsplit(rckb, ' ')))
      mck <- as.numeric(unlist(strsplit(mckb, ' ')))

      alpha <- as.numeric(unlist(strsplit(alphab, ' ')))
      icc   <- as.numeric(unlist(strsplit(iccb, ' ')))
      rho   <- as.numeric(unlist(strsplit(rhob, ' ')))
      phi   <- as.numeric(unlist(strsplit(phib, ' ')))
      r2yx  <- as.numeric(unlist(strsplit(r2yxb, ' ')))
      r2tx  <- as.numeric(unlist(strsplit(r2txb, ' ')))
      deff_wgt <- as.numeric(unlist(strsplit(deff_wgtb, ' ')))

      # CHECK FOR ERRORS

      # Set RHO and PHI to 0 for AR1 = 3 and change AR1 to 1 and 0 to fit older coding

      ar1_orig <- ar1
      if (ar1==3) {
        rho <- 0
        phi <- 0
      }
      if (ar1>1) {
        ar1 <- 0
      }

      # Write error messages and crash variable for key inputs
      # Test for missing values, decimal values for categories, and out of ranges

      crash <- 0
      nerr  <- 0
      err_mess <- matrix(ncol=100,nrow=1)

      # NTIMEP

      bad_per <- 0
      if ((!exists("ntimep")) | (is.na(ntimep)) | (ntimep<=0) | (ntimep == "")
          | (ntimep != round(ntimep))) {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid number of time periods")
        bad_per <- 1
      } else if ((did_cits>1) & (type_cits==1) & (ntimep<6))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("The number of time periods must be at least 6 for the fully-interacted
                              CITS and ITS designs")
        bad_per <- 1
      } else if ((did_cits>1) & (type_cits==2) & (ntimep<6))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("The number of time periods must be at least 6 for CITS and ITS
                              designs with common slopes")
        bad_per <- 1
      } else if ((did_cits>1) & (type_cits==3) & (ntimep<4))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("The number of time periods must be at least 4 for the
                              CITS and ITS designs with discrete post-period indicators")
        bad_per <- 1
      }

      # NTIMEG and S1-SK

      bad_time <- 0
      if ((!exists("ntimeg")) | (is.na(ntimeg)) | (ntimeg<1)
          | (ntimeg != round(ntimeg))) {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid number of timing groups")
        bad_time <- 1
      } else if ((ntimeg>0) & (ntimeg != length(sk)))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- sprintf("The number of elements in the list showing the
        start times for each treatment timing group, %d, does not match the number of
        timing groups, %d.",length(sk),ntimeg)

        bad_time <- 1
      } else if (ntimeg>0)
      {
        bads <- 0
        bad_aftp <- 0
        for (i in 1:ntimeg) {
          if ((is.na(sk[i])) | (sk[i]<=1) | ((sk[i]>ntimep) & (!is.na(ntimep)))) {
            crash <- 1
            nerr  <- nerr + 1
            err_mess[1,nerr] <- sprintf("Invalid start time for timing group %d",i)
            bad_time <- 1
          } else if (sk[i]>1)
          {
            if (bad_per==0) {
              aftp <- ntimep - sk[i] + 1
              if ((did_cits>1) & (type_cits<3) & (aftp<3)) {
                bad_aftp <- 1
              }
            }
            if (sk[i] <= bads) {
              crash <- 1
              nerr  <- nerr + 1
              err_mess[1,nerr] <- sprintf("Invalid start time for timing group %d:
                                          it does not occur after the previous one",i)
              bad_time <- 1
            }
            bads <- sk[i]
          }
        }
        if ((did_cits>1) & (type_cits>=1) & (sk[1]<4)) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- c("The start times must be at least 4 for the CITS and ITS designs")
        }
        if (bad_aftp==1) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- c("The start times must allow for at least 3 post-periods
                                for the CITS and ITS designs with post-period trendlines")
        }
      }

      # Q_TIME and L_TIME

      bad_q <- 0
      if ((avg_point==2) & (cal_expos==1)) {
        if ((!exists("q_time")) | (is.na(q_time)) | (q_time<=0) | (q_time == "")
            | (q_time != round(q_time)) | ((q_time>ntimep) & (!is.na(ntimep)))) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- c("Invalid post-period time point for the
                                power analysis")
          bad_q <- 1
        }
        if ((bad_q==0) & (bad_time==0) & (q_time<sk[1])) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- c("Invalid post-period time point for the
                                power analysis: it must be a post-period for at
                                least one timing group")
        }
      }

      bad_l <- 0
      if ((avg_point==2) & (cal_expos==2)) {
        if ((!exists("l_time")) | (is.na(l_time)) | (l_time<=0) | (l_time == "")
            | (l_time != round(l_time))) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- c("Invalid exposure time point for the
                                power analysis")
          bad_l <- 1
        }
        if ((bad_l==0) & (bad_time==0)) {
          l_cal <- l_time + sk[1] - 1
          if ((l_cal < sk[1]) | ((l_cal > ntimep) & (!is.na(ntimep)))) {
            crash <- 1
            nerr  <- nerr + 1
            err_mess[1,nerr] <- c("Invalid exposure time point for the
                                power analysis: it must correspond to an observed
                                post-period for at least one timing group")
          }
        }
      }

      # NSAMP

      if ((!exists("nsamp")) | (is.na(nsamp)) | (nsamp<1)
          | (nsamp != round(nsamp))) {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid number of individuals")
      }

      # MT, MTk, MC, and MCk

      if (did_cits < 3) {
        mclus <- mt + mc
      } else if (did_cits == 3)
      {
        mclus <- mt
      }

      if (samp_size==1) {
        if ((!exists("mt")) | (is.na(mt)) | (mt<=1)) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- c("Invalid number of treatment clusters")
        }
        if ((!exists("mc")) | (is.na(mc)) | (mc<=1)) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- c("Invalid number of comparison clusters")
        }
        if (mt != round(mt)) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- c("Invalid number of treatment clusters")
        }
        if (mc != round(mc)) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- c("Invalid number of comparison clusters")
        }

        if ((!is.na(ntimeg)) & (ntimeg != length(mtk))) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- sprintf("The number of elements in the list showing the number of
          treatment clusters in each timing group, %d, does not match
          the number of timing groups, %d.",length(mtk),ntimeg)

          bad_time <- 1
        } else if ((!is.na(mt)) & (mt>0) & (!is.na(ntimeg)) & (ntimeg>0) & (bad_time==0))
        {
          sum_mt = 0;
          for (i in 1:ntimeg) {
            if ((is.na(mtk[i])) | (mtk[i]<=1) | (mtk[i]>mt) | (mtk[i]!=round(mtk[i]))) {
              crash <- 1
              nerr  <- nerr + 1
              err_mess[1,nerr] <- sprintf("Invalid number of treatment clusters in timing group %d",i)
            } else
            {
              sum_mt = sum_mt + mtk[i]
            }
          }

          if (sum_mt != mt) {
            crash <- 1
            nerr  <- nerr + 1
            err_mess[1,nerr] <- sprintf("The number of treatment clusters across the timing groups
            does not sum to the total number of treatment clusters, %d.",mt)
          }
        }

        if (did_cits < 3) {
          if ((!is.na(ntimeg)) & (ntimeg != length(mck))) {
            crash <- 1
            nerr  <- nerr + 1
            err_mess[1,nerr] <- sprintf("The number of elements in the list showing the number of
            comparison clusters in each timing group, %d, does not match
            the number of timing groups, %d.",length(mck),ntimeg)

            bad_time <- 1
          } else if ((!is.na(mc)) & (mc>0) & (!is.na(ntimeg)) & (ntimeg>0) & (bad_time==0))
          {
            sum_mc = 0;
            for (i in 1:ntimeg) {
              if ((is.na(mck[i])) | (mck[i]<=1) | (mck[i]>mc) | (mck[i]!=round(mck[i]))) {
                crash <- 1
                nerr  <- nerr + 1
                err_mess[1,nerr] <- sprintf("Invalid number of comparison clusters in timing group %d",i)
              } else
              {
                sum_mc = sum_mc + mck[i]
              }
            }

            if (sum_mc != mc) {
              crash <- 1
              nerr  <- nerr + 1
              err_mess[1,nerr] <- sprintf("The number of comparison clusters across the timing groups
                                          does not sum to the total number of comparison clusters, %d.",mc)
            }
          }
        }
      }

      # MDE, RT, RTk, and RCk

      if (samp_size == 2) {

        if ((!exists("mde")) | (is.na(mde)) | (mde <= 0)) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- c("Invalid MDE")
        }

        if ((did_cits==1) | (did_cits==2)) {
          if ((!exists("rt")) | (is.na(rt)) | (rt<=0) | (rt>=1)) {
            crash <- 1
            nerr  <- nerr + 1
            err_mess[1,nerr] <- c("Invalid proportion of all clusters that are treatment clusters")
          }
        } else if (did_cits == 3)
        {
          rt <- .5
        }

        if ((!is.na(rt)) & (rt>0) & (rt<1)) {
          if ((!is.na(ntimeg)) & (ntimeg != length(rtk))) {
            crash <- 1
            nerr  <- nerr + 1
            err_mess[1,nerr] <- sprintf("The number of elements in the list showing the proportions of
            treatment clusters in each timing group, %d, does not match the
            number of timing groups, %d.",length(rtk),ntimeg)
            bad_time <- 1
          } else if ((!is.na(ntimeg)) & (ntimeg == length(rtk)))
          {
            sum_rt <- 0
            for (i in 1:ntimeg) {
              if ((is.na(rtk[i])) | (rtk[i]<=0) | (rtk[i]>1)) {
                crash <- 1
                nerr  <- nerr + 1
                err_mess[1,nerr] <- sprintf("Invalid proportion of treatment clusters in timing group %d",i)
              } else
              {
                sum_rt <- sum_rt + rtk[i]
              }
            }

            if (sum_rt!=1) {
              crash <- 1
              nerr  <- nerr + 1
              err_mess[1,nerr] <- c("The proportions of treatment clusters in the timing groups do not sum to 1")
            }
          }

          if (did_cits < 3) {
            if ((!is.na(ntimeg)) & (ntimeg != length(rck))) {
              crash <- 1
              nerr  <- nerr + 1
              err_mess[1,nerr] <- sprintf("The number of elements in the list showing the proportions of
              comparison clusters in each timing group, %d, does not match the
                                          number of timing groups, %d.",length(rck),ntimeg)
              bad_time <- 1
            } else if ((!is.na(ntimeg)) & (ntimeg == length(rck)))
            {
              sum_rc <- 0
              for (i in 1:ntimeg) {
                if ((is.na(rck[i])) | (rck[i]<=0) | (rck[i]>1)) {
                  crash <- 1
                  nerr  <- nerr + 1
                  err_mess[1,nerr] <- sprintf("Invalid proportion of comparison clusters in timing group %d",i)
                } else
                {
                  sum_rc <- sum_rc + rck[i]
                }
              }

              if (sum_rc!=1) {
                crash <- 1
                nerr  <- nerr + 1
                err_mess[1,nerr] <- c("The proportions of comparison clusters in the timing groups do not sum to 1")
              }
            }
          }
        }
      }

      # EVEN_TIME = 2 with TP

      if ((even_time==2) & (!is.na(ntimep)) & (ntimep>0) & (bad_per==0)) {
        badp = 0

        if (ntimep != length(tp)) {
          crash <- 1
          nerr  <- nerr + 1
          err_mess[1,nerr] <- sprintf("The number of elements in the list showing the time intervals, %d,
            does not match the number of time periods, %d.",length(tp),ntimep)
        } else if (ntimep == length(tp))
        {
        for (i in 1:ntimep) {

          if ((is.na(tp[i])) | (tp[i]<=0) | (tp[i] != round(tp[i]))) {
            crash <- 1
            nerr  <- nerr + 1
            err_mess[1,nerr] <- sprintf("Invalid time interval for period %d: must be a positive integer",i)
          } else if (tp[i]>1)
          {
            if (tp[i] <= badp) {
              crash <- 1
              nerr  <- nerr + 1
              err_mess[1,nerr] <- sprintf("Invalid time interval for period %d: Occurs before the last one",i)
            }
            badp = tp[i]
          }
        }
        }
      }

      # ALPHA

      if ((!exists("alpha")) | (length(alpha) != 1)) {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid alpha level")
      } else if ((is.na(alpha)) | (alpha>=1) | (alpha<=0))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid alpha level")
      }

      # ICC

      if ((!exists("icc")) | (length(icc) != 1)) {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid intraclass correlation")
      } else if ((is.na(icc)) | (icc > 1) | (icc < 0))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid intraclass correlation")
      }

      # RHO

      if ((!exists("rho")) | (length(rho) != 1)) {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid cluster-level autocorrelation")
      } else if ((is.na(rho)) | (rho >= 1) | (rho <= -1))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid cluster-level autocorrelation")
      }

      # PHI

      if ((!exists("phi")) | (length(phi) != 1)) {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid individual-level autocorrelation")
      } else if ((is.na(phi)) | (phi >= 1) | (phi <= -1))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid individual-level autocorrelation")
      }


      # DEFF_WGT

      if ((!exists("deff_wgt")) | (length(deff_wgt) != 1)) {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid design effect due to weighting")
      } else if ((is.na(deff_wgt)) | (deff_wgt < 1))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid design effect due to weighting")
      }


      # R2YX

      if ((!exists("r2yx")) | (length(r2yx) != 1)) {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid regression R2 value from model covariates")
      } else if ((is.na(r2yx)) | (r2yx >= 1) | (r2yx < 0))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid regression R2 value from model covariates")
      }

      # R2TX

      if ((!exists("r2tx")) | (length(r2tx) != 1)) {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid correlation between model treatment indicators and covariates")
      } else if ((is.na(r2tx)) | (r2tx >= 1) | (r2tx < 0))
      {
        crash <- 1
        nerr  <- nerr + 1
        err_mess[1,nerr] <- c("Invalid correlation between model treatment indicators and covariates")
      }

      # NOW PERFORM THE VARIANCE-COVARIANCE CALCULATIONS IF crash==0

      if (crash==0) {

        # SET MT, MTK, MC, and MCK variables if SAMP_SIZE = 2 and
        # RT, RTK, RC, and RCK Variables if SAMP_SIZE = 1

        if (samp_size==2) {
          mclus <- 1
          mt    <- rt
          mtk_its <- rtk
          mtk     <- rt*rtk
          mc      <- 1 - rt
          mck     <- (1-rt)*rck
        } else if (samp_size==1)
        {
          mde <- .2
          rt  <- .5
          rc  <- 1 - rt
          mtk_its <- mtk
        }

        # Calculate the number of controls

        # mc <- mclus-mt

        # Calculate the t cutoff alpha depending on a two- or one-tailed test

        if (two_tailed==1) {
          alpha2 <- 1 - (alpha/2)
        } else if (two_tailed==0)
        {
          alpha2 <- 1 - alpha
        }

        # Create time variables if evenly spaced

        if (even_time==2) {
          tpp <- tp
        } else if (even_time==1) {
          tpp <- matrix(0,ntimep,1)
          for (i in 1:ntimep)
          tpp[i] <- i
        }

        #print(paste('tp',tp))
        #print(paste('tpp',tpp))

        # Create variable that has variation in treatment timing- that is, >1 timing group

        if (ntimeg>1) {
          has_var_tt <- 1
        } else if (ntimeg==1)
        {
          has_var_tt <- 0
        }

        # For CITS design - get tbar for the baseline period for each timing group
        # as well as sums of squares in baseline period and post-period sums from baseline means

        tbar  <- sk
        ssqt  <- sk
        tbara <- sk
        ssqta <- sk
        pp_diff <- sk

        for (k in 1:ntimeg) {

          stk   <- sk[k]
          timep <- tpp[1:(stk-1)]                # Pull off baseline period
          tbar[k]   <- mean(as.numeric(timep))    # Take mean of time variable in baseline period
          pre_diff2 <- (timep - tbar[k])^2
          ssqt[k]   <- sum(as.numeric(pre_diff2)) # Sum of squares of time variable in baseline period

          timepa <- tpp[stk:ntimep]                # Pull off post-period
          tbara[k]   <- mean(as.numeric(timepa))    # Take mean of time variable in post-period
          post_diff2 <- (timepa - tbara[k])^2
          ssqta[k]   <- sum(as.numeric(post_diff2)) # Sum of squares of time variable in post-period

          pp_diff[k] <- tbara[k]-tbar[k]

          #print(paste('postd',postd))

        } # end k loop

        tbar_full   <- sum(as.numeric(tpp))/ntimep  # Full period time mean
        diff2_full  <- (tpp - tbar_full)^2
        ssqt_full   <- sum(as.numeric(diff2_full)) # Sum of squares of time variable in full-period

        #print(paste('tpp',tpp))
        #print(paste('tbar_full',tbar_full))
        #print(paste('ssqt_full',ssqt_full))

        # Calculate net R2 value that also includes the design effects due to weighting

        if ((r2yx >= 0) & (r2yx < 1) & (r2tx >= 0) & (r2tx < 1)) {
          r2 <- (1-r2yx)/(1-r2tx)
        } else
        {
          r2 <- 1
        }

        # Define BK, AK, and the weight vector and get sums of weights for later

        bk <- sk-1
        ak <- ntimep-sk+1
        sumak <- sum(as.numeric(ak))

        wk <- sk

        if (avg_point==1) {
          wk <- ak
        } else if (avg_point==2) {
          wk <- sk - sk + 1
        }

        #print(paste("wk",wk))

        sumwk <- sum(as.numeric(ak))  # I fixed this since weights are ak

        #
        # VARIANCE CALCULATIONS
        #

        term1      <- 0
        term1_long <- 0
        term1_noc  <- 0
        term1_long_noc <- 0

        term1q      <- 0
        term1q_long <- 0
        term1q_noc  <- 0
        term1q_long_noc <- 0

        sumiq  <- 0
        mclusq <- 0
        aclusq <- 0

        cits_term1 <- 0
        its_term1  <- 0
        cits_term1_long <- 0
        its_term1_long  <- 0

        rho_pre  <- sk
        rho_post <- sk
        rho_pp   <- sk
        phi_pre  <- sk
        phi_post <- sk
        phi_pp   <- sk

        for (k in 1:ntimeg) {

          # For DID, cross-sectional design

          #calc_rho(corrp=&rho,ar1=0,prepost=1,pres=1,pref=&&b&k,posts=0,postf=0,tbarb1=0,tbarb2=0)

          k <- as.numeric(k)
          bkk <- bk[k]
          skk <- sk[k]
          akk <- ak[k]
          wkk <- wk[k]
          mtkk <- mtk[k]
          mckk <- mck[k]
          mtkk_its <- mtk_its[k]
          tbarkk   <- tbar[k]
          ssqtkk   <- ssqt[k]
          tbarakk   <- tbara[k]
          ssqtakk   <- ssqta[k]
          pp_diffkk <- pp_diff[k]

          iqkk <- 1

          # DID for average post-period

          if (avg_point==1) {

            rho_pre  <- calc_rho(rho,1,1,bkk,0,0,0,0,ar1,tpp)
            rho_post <- calc_rho(rho,2,0,0,skk,ntimep,0,0,ar1,tpp)
            rho_pp   <- calc_rho(rho,3,1,bkk,skk,ntimep,0,0,ar1,tpp)

            theta_term1 <- (1 + rho_post*(akk-1))/akk
            theta_term2 <- (1 + rho_pre*(bkk-1))/bkk
            theta_term3 <- rho_pp
            theta_term  <- icc*(theta_term1 + theta_term2 - 2*theta_term3)

            # e_term is for the cross-section design
            e_term <- ((1-icc)/nsamp)*((1/akk) + (1/bkk))

            term1 <- term1 + (akk^2)*((1/mtkk)+(1/mckk))*r2*(theta_term + e_term)

            # For the ITS design
            term1_noc <- term1_noc + (akk^2)*((1/mtkk_its))*r2*(theta_term + e_term)

            if (cross_sec==2) {

              phi_pre  <- calc_rho(phi,1,1,bkk,0,0,0,0,ar1,tpp)
              phi_post <- calc_rho(phi,2,0,0,skk,ntimep,0,0,ar1,tpp)
              phi_pp   <- calc_rho(phi,3,1,bkk,skk,ntimep,0,0,ar1,tpp)

              e_term1  <- (1 + phi_post*(akk-1))/akk
              e_term2  <- (1 + phi_pre*(bkk-1))/bkk
              e_term3  <- phi_pp
              e_term_long <- ((1-icc)/nsamp)*(e_term1 + e_term2 - 2*e_term3);

              term1_long <- term1_long + (akk^2)*((1/mtkk)+(1/mckk))*r2*(theta_term + e_term_long)
              term1_long_noc <- term1_long_noc + (akk^2)*((1/mtkk_its))*r2*(theta_term + e_term_long)
            }
          }

          # DID for a specific post-period point

          if (avg_point==2) {

            if (cal_expos==1) {
              qkk <- q_time
            } else if (cal_expos==2) {
              qkk <- l_time + skk - 1
            }

            if ((qkk>=skk) & (qkk<=ntimep)) {
              iqkk <- 1
              sumiq <- sumiq + 1
              mclusq <- mclusq + (mtkk+mckk)
              aclusq <- aclusq + akk
            } else {
              iqkk <- 0
            }

            if (iqkk==1) {

              rho_preq  <- calc_rho(rho,1,1,bkk,0,0,0,0,ar1,tpp)
              rho_postq <- calc_rho(rho,2,0,0,skk,ntimep,0,0,ar1,tpp)
              rho_ppq   <- calc_rho(rho,3,1,bkk,qkk,qkk,0,0,ar1,tpp)

              theta_term1q <- 1
              theta_term2q <- (1 + rho_preq*(bkk-1))/bkk
              theta_term3q <- rho_ppq
              theta_termq  <- icc*(theta_term1q + theta_term2q - 2*theta_term3q)

              e_termq <- ((1-icc)/nsamp)*(1 + (1/bkk))

              term1q <- term1q + ((1/mtkk)+(1/mckk))*r2*(theta_termq + e_termq)

              #print(paste('DID theta_termq',theta_termq))
              #print(paste('DID e_termq',e_termq))

              # For the ITS design
              term1q_noc <- term1q_noc + ((1/mtkk_its))*r2*(theta_termq + e_termq)

              # This repeats this for the DID variance needed for the CITS design

              rho_pre  <- calc_rho(rho,1,1,bkk,0,0,0,0,ar1,tpp)
              rho_post <- calc_rho(rho,2,0,0,skk,ntimep,0,0,ar1,tpp)
              rho_pp   <- calc_rho(rho,3,1,bkk,skk,ntimep,0,0,ar1,tpp)

              theta_term1 <- (1 + rho_post*(akk-1))/akk
              theta_term2 <- (1 + rho_pre*(bkk-1))/bkk
              theta_term3 <- rho_pp
              theta_term  <- icc*(theta_term1 + theta_term2 - 2*theta_term3)

              e_term <- ((1-icc)/nsamp)*((1/akk) + (1/bkk))

              term1 <- term1 + ((1/mtkk)+(1/mckk))*r2*(theta_term + e_term)
              term1_noc <- term1_noc + ((1/mtkk_its))*r2*(theta_term + e_term)

              #print(paste('theta_term',theta_term))

              if (cross_sec==2) {

                phi_pre  <- calc_rho(phi,1,1,bkk,0,0,0,0,ar1,tpp)
                phi_post <- calc_rho(phi,2,0,0,skk,ntimep,0,0,ar1,tpp)
                phi_pp   <- calc_rho(phi,3,1,bkk,qkk,qkk,0,0,ar1,tpp)

                e_term1q  <- 1
                e_term2q  <- (1 + phi_pre*(bkk-1))/bkk
                e_term3q  <- phi_pp
                e_term_longq <- ((1-icc)/nsamp)*(e_term1q + e_term2q - 2*e_term3q)

                term1q_long <- term1q_long + ((1/mtkk)+(1/mckk))*r2*(theta_termq + e_term_longq)
                term1q_long_noc <- term1q_long_noc + ((1/mtkk_its))*r2*(theta_termq + e_term_longq)

                e_term1 <- (1 + phi_post*(akk-1))/akk
                e_term2 <- (1 + phi_pre*(bkk-1))/bkk
                e_term3 <- phi_pp
                e_term_long <- ((1-icc)/nsamp)*(e_term1 + e_term2 - 2*e_term3)

                term1_long <- term1_long + ((1/mtkk)+(1/mckk))*r2*(theta_term + e_term_long)
                term1_long_noc <- term1_long_noc + ((1/mtkk_its))*r2*(theta_term + e_term_long)

              }
            } # iqkk==1
          } # avg_point==2

          # CITS - fully interacted model

          if ((did_cits>1) & (type_cits==1) & (avg_point==1)) {

            rho_pre1  <- calc_rho(rho,6,1,bkk,0,0,tbarkk,0,ar1,tpp)
            rho_pre2  <- calc_rho(rho,7,1,bkk,0,0,tbarkk,0,ar1,tpp)
            rho_pp1   <- calc_rho(rho,8,1,bkk,skk,ntimep,tbarkk,0,ar1,tpp)

            # CITS, Cross-section

            cits_theta_term1 <- (pp_diffkk^2)*((1/ssqtkk)+((bkk-1)*bkk*rho_pre1/(ssqtkk^2)))
            cits_theta_term2 <- 2*pp_diffkk*bkk*rho_pre2 / ssqtkk
            cits_theta_term3 <- 2*pp_diffkk*bkk*rho_pp1 / ssqtkk
            cits_theta_term  <- icc*(cits_theta_term1 + cits_theta_term2 - cits_theta_term3)

            #print(paste('cits_theta_term',cits_theta_term))

            # CITS, e_term for the cross-section design

            cits_e_term4 <- (pp_diffkk^2) / ssqtkk
            cits_e_term  <- ((1-icc)/nsamp)*(cits_e_term4)

            #print(paste('cits_e_term4',cits_e_term4))

            # CITS and ITS AGGREGATE

            cits_term1 <- cits_term1 + (akk^2)*((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term)
            its_term1  <- its_term1  + (akk^2)*((1/mtkk_its))*r2*(cits_theta_term + cits_e_term)

            # CITS e_term_long is for the longitudinal design

            if (cross_sec==2) {

              phi_pre1  <- calc_rho(phi,6,1,bkk,0,0,tbarkk,0,ar1,tpp)
              phi_pre2  <- calc_rho(phi,7,1,bkk,0,0,tbarkk,0,ar1,tpp)
              phi_pp1   <- calc_rho(phi,8,1,bkk,skk,ntimep,tbarkk,0,ar1,tpp)

              cits_e_term1 <- (pp_diffkk^2)*((1/ssqtkk)+((bkk-1)*bkk*phi_pre1/(ssqtkk^2)))
              cits_e_term2 <- 2*pp_diffkk*bkk*phi_pre2 / ssqtkk
              cits_e_term3 <- 2*pp_diffkk*bkk*phi_pp1 / ssqtkk
              cits_e_term_long  <- ((1-icc)/nsamp)*(cits_e_term1 + cits_e_term2 - cits_e_term3)

              #print(paste('phi_pre1',phi_pre1))
              #print(paste('rho_pre1',rho_pre1))
              #print(paste('ssqtkk',ssqtkk))
              #print(paste('cits_e_term',cits_e_term))
              #print(paste('cits_e_term1',cits_e_term1))
              #print(paste('cits_e_term2',cits_e_term2))
              #print(paste('cits_e_term3',cits_e_term3))
              #print(paste('cits_e_term_long',cits_e_term_long))

              cits_term1_long <- cits_term1_long + (akk^2)*((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term_long)
              its_term1_long  <- its_term1_long  + (akk^2)*((1/mtkk_its))*r2*(cits_theta_term + cits_e_term_long)

            }   # end cross_sec==2
          }   # end (did_cits>1) & (type_cits==1) & (avg_point==1)

          if ((did_cits>1) & (type_cits==1) & (avg_point==2) & (iqkk==1)) {

            rho_pre1   <- calc_rho(rho,6,1,bkk,0,0,tbarkk,0,ar1,tpp)
            rho_pre2   <- calc_rho(rho,7,1,bkk,0,0,tbarkk,0,ar1,tpp)
            rho_post1  <- calc_rho(rho,9,0,0,skk,ntimep,tbarakk,0,ar1,tpp)
            rho_post2  <- calc_rho(rho,10,0,0,skk,ntimep,tbarakk,0,ar1,tpp)
            rho_pp2    <- calc_rho(rho,11,1,bkk,skk,ntimep,tbarakk,0,ar1,tpp)
            rho_pp3    <- calc_rho(rho,8,1,bkk,skk,ntimep,tbarkk,0,ar1,tpp)
            rho_pp4    <- calc_rho(rho,12,1,bkk,skk,ntimep,tbarkk,tbarakk,ar1,tpp)

            diffa_qkk <- tpp[qkk]-tbarakk
            diff_qkk  <- tpp[qkk]-tbarkk
            cits_theta_term1 <- (diffa_qkk^2)*((1/ssqtakk)+((akk-1)*akk*rho_post1/(ssqtakk^2)))
            cits_theta_term2 <- (diff_qkk^2)*((1/ssqtkk)+((bkk-1)*bkk*rho_pre1/(ssqtkk^2)))
            cits_theta_term3 <- 2*diffa_qkk*akk*rho_post2/ ssqtakk
            cits_theta_term4 <- 2*diff_qkk*bkk*rho_pre2/ ssqtkk
            cits_theta_term5 <- 2*diffa_qkk*akk*rho_pp2/ ssqtakk
            cits_theta_term6 <- 2*diff_qkk*bkk*rho_pp3/ ssqtkk
            cits_theta_term7 <- 2*diffa_qkk*diff_qkk*akk*bkk*rho_pp4/ (ssqtakk*ssqtkk)

            cits_theta_term  <- icc*(cits_theta_term1 + cits_theta_term2 + cits_theta_term3 + cits_theta_term4
                                     - cits_theta_term5 - cits_theta_term6 - cits_theta_term7)

            #print(paste('cits_theta_term',cits_theta_term))

            cits_e_term8 <- (diffa_qkk^2)*(1/ssqtakk)
            cits_e_term9 <- (diff_qkk^2)*(1/ssqtkk)
            cits_e_term  <- ((1-icc)/nsamp)*(cits_e_term8 + cits_e_term9)

            #print(paste('cits_e_term',cits_e_term))

            cits_term1 <- cits_term1 + ((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term)
            its_term1  <- its_term1  + ((1/mtkk_its))*r2*(cits_theta_term + cits_e_term)

            if (cross_sec==2) {

              phi_pre1   <- calc_rho(phi,6,1,bkk,0,0,tbarkk,0,ar1,tpp)
              phi_pre2   <- calc_rho(phi,7,1,bkk,0,0,tbarkk,0,ar1,tpp)
              phi_post1  <- calc_rho(phi,9,0,0,skk,ntimep,tbarakk,0,ar1,tpp)
              phi_post2  <- calc_rho(phi,10,0,0,skk,ntimep,tbarakk,0,ar1,tpp)
              phi_pp2    <- calc_rho(phi,11,1,bkk,skk,ntimep,tbarakk,0,ar1,tpp)
              phi_pp3    <- calc_rho(phi,8,1,bkk,skk,ntimep,tbarkk,0,ar1,tpp)
              phi_pp4    <- calc_rho(phi,12,1,bkk,skk,ntimep,tbarkk,tbarakk,ar1,tpp)

              diffa_qkk <- tpp[qkk]-tbarakk
              diff_qkk  <- tpp[qkk]-tbarkk
              cits_e_term1 <- (diffa_qkk^2)*((1/ssqtakk)+((akk-1)*akk*phi_post1/(ssqtakk^2)))
              cits_e_term2 <- (diff_qkk^2)*((1/ssqtkk)+((bkk-1)*bkk*phi_pre1/(ssqtkk^2)))
              cits_e_term3 <- 2*diffa_qkk*akk*phi_post2/ ssqtakk
              cits_e_term4 <- 2*diff_qkk*bkk*phi_pre2/ ssqtkk
              cits_e_term5 <- 2*diffa_qkk*akk*phi_pp2/ ssqtakk
              cits_e_term6 <- 2*diff_qkk*bkk*phi_pp3/ ssqtkk
              cits_e_term7 <- 2*diffa_qkk*diff_qkk*akk*bkk*phi_pp4/ (ssqtakk*ssqtkk)

              cits_e_term_long  <- ((1-icc)/nsamp)*(cits_e_term1 + cits_e_term2 + cits_e_term3 + cits_e_term4
                                                    - cits_e_term5 - cits_e_term6 - cits_e_term7)

              cits_term1_long <- cits_term1_long + ((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term_long)
              its_term1_long  <- its_term1_long  + ((1/mtkk_its))*r2*(cits_theta_term + cits_e_term_long)

            } # end cross_sec==2

          } # end (did_cits>1) & (type_cits==1) & (avg_point==2) & (iqkk==1)

          # CITS Common Slopes

          if ((did_cits>1) & (type_cits==2)) {

            postd <- tpp               # Create post-period indicator
            postd[1:(skk-1)] <- 0
            postd[skk:ntimep] <- 1

            postmean   <- akk/ntimep #mean(as.numeric(postd))
            #chk_mean   <- akk/ntimep

            #print(paste('postmean',postmean))
            #print(paste('postd',postd))
            #print(paste('chk_mean',chk_mean))

            # Calculate rho_full since involves new parameters and don't want to update all the old functions to add new inputs

            rho_full1 <- 0
            atemp <- ntimep
            if (atemp>1) {

              vdf1 <- as.numeric(seq(1:atemp))
              vdf  <- expand.grid(vdf1,vdf1)
              vv   <- vdf[which(vdf[,2]>vdf[,1]),]

              time_term <- (postd[vv[,1]] - as.numeric(postmean))*(postd[vv[,2]] - as.numeric(postmean))

              if (ar1==1) {
                sv <- time_term*((rho)^(tpp[vv[,2]]-tpp[vv[,1]]))
                rho_full1 <- sum(sv)
              } else if (ar1==0)
              {
                sv <- time_term*(rho)
                rho_full1 <- sum(sv)
              }

            atemp1 <- atemp-1
            rho_full1 <- 2*rho_full1/(atemp*atemp1)

            #print(paste('rho_full1',rho_full1))

            }

            rho_full2 <- 0
            atemp <- ntimep

            if (atemp>1) {
              vdf1 <- as.numeric(seq(1:ntimep))
              vdf  <- expand.grid(vdf1,vdf1)
              vv   <- vdf[which(vdf[,2]!=vdf[,1]),]

              time_term <- (tpp[vv[,1]] - as.numeric(tbar_full))*(postd[vv[,2]] - as.numeric(postmean))
              if (ar1==1) {
                sv <- time_term*((rho)^abs((tpp[vv[,2]]-tpp[vv[,1]])))
                rho_full2 <- sum(sv)
              } else if (ar1==0)
              {
                sv <- time_term*(rho)
                rho_full2 <- sum(sv)
              }
              atemp1 <- atemp-1
              rho_full2 <- rho_full2/(atemp*atemp1)

              #print(paste('rho_full2',rho_full2))

            }

            rho_full3   <- calc_rho(rho,6,1,ntimep,0,0,tbar_full,0,ar1,tpp)

            #print(paste('rho_full3',rho_full3))

            ssqt_term        <- ssqt_full / (ssqtkk+ssqtakk)
            cits_theta_term1 <- ssqt_term
            cits_theta_term2 <- ((1/akk)+(1/bkk))*ntimep*(ntimep-1)*(ssqt_term^2)*rho_full1
            cits_theta_term3 <- 2*ntimep*(ntimep-1)*(ssqt_full/(ssqtkk+ssqtakk))*(pp_diffkk/(ssqtkk+ssqtakk))*rho_full2
            cits_theta_term4 <- akk*bkk*(ntimep-1)*((pp_diffkk/(ssqtkk+ssqtakk))^2)*rho_full3

            cits_theta_term  <- icc*((1/akk)+(1/bkk))*(cits_theta_term1 + cits_theta_term2 - cits_theta_term3 + cits_theta_term4)

            cits_e_term1     <-  ssqt_term
            cits_e_term      <- ((1-icc)/nsamp)*((1/akk)+(1/bkk))*(cits_e_term1)

            if (avg_point==1) {
              cits_term1 <- cits_term1 + (akk^2)*((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term)
              its_term1  <- its_term1  + (akk^2)*((1/mtkk_its))*r2*(cits_theta_term + cits_e_term)

              #print(paste('cits_theta_term1',cits_theta_term1))
              #print(paste('cits_theta_term2',cits_theta_term2))
              #print(paste('cits_theta_term3',cits_theta_term3))
              #print(paste('cits_theta_term4',cits_theta_term4))
              #print(paste('rho_full1',rho_full1))
              #print(paste('rho_full2',rho_full2))
              #print(paste('rho_full3',rho_full3))
              #print(paste('ssqt_term',ssqt_term))

              #print(paste('cits_term1',cits_term1))
              #print(paste('its_term1',its_term1))

            } else if ((avg_point==2) & (iqkk==1)) {
              cits_term1 <- cits_term1 + ((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term)
              its_term1  <- its_term1  + ((1/mtkk_its))*r2*(cits_theta_term + cits_e_term)
            }

            # Longitudinal

            if (cross_sec==2) {

              phi_full1 <- 0
              atemp <- ntimep
              if (atemp>1) {

                vdf1 <- as.numeric(seq(1:atemp))
                vdf  <- expand.grid(vdf1,vdf1)
                vv   <- vdf[which(vdf[,2]>vdf[,1]),]

                time_term <- (postd[vv[,1]] - as.numeric(postmean))*(postd[vv[,2]] - as.numeric(postmean))

                if (ar1==1) {
                  sv <- time_term*((phi)^(tpp[vv[,2]]-tpp[vv[,1]]))
                  phi_full1 <- sum(sv)
                } else if (ar1==0)
                {
                  sv <- time_term*(phi)
                  phi_full1 <- sum(sv)
                }

                atemp1 <- atemp-1
                phi_full1 <- 2*phi_full1/(atemp*atemp1)

              }

              phi_full2 <- 0
              atemp <- ntimep

              if (atemp>1) {
                vdf1 <- as.numeric(seq(1:ntimep))
                vdf  <- expand.grid(vdf1,vdf1)
                vv   <- vdf[which(vdf[,2]!=vdf[,1]),]

                time_term <- (tpp[vv[,1]] - as.numeric(tbar_full))*(postd[vv[,2]] - as.numeric(postmean))
                if (ar1==1) {
                  sv <- time_term*((phi)^abs((tpp[vv[,2]]-tpp[vv[,1]])))
                  phi_full2 <- sum(sv)
                } else if (ar1==0)
                {
                  sv <- time_term*(phi)
                  phi_full2 <- sum(sv)
                }
                atemp1 <- atemp-1
                phi_full2 <- phi_full2/(atemp*atemp1)

                #print(paste('rho_full2',rho_full2))

              }

              phi_full3   <- calc_rho(phi,6,1,ntimep,0,0,tbar_full,0,ar1,tpp)

              #print(paste('rho_full3',rho_full3))

              #ssqt_term        <- ssqt_full / (ssqtkk+ssqtakk)
              cits_e_term1 <- ssqt_term
              cits_e_term2 <- ((1/akk)+(1/bkk))*ntimep*(ntimep-1)*(ssqt_term^2)*phi_full1
              cits_e_term3 <- 2*ntimep*(ntimep-1)*(ssqt_full/(ssqtkk+ssqtakk))*(pp_diffkk/(ssqtkk+ssqtakk))*phi_full2
              cits_e_term4 <- akk*bkk*(ntimep-1)*((pp_diffkk/(ssqtkk+ssqtakk))^2)*phi_full3

              cits_e_term_long  <- ((1-icc)/nsamp)*((1/akk)+(1/bkk))*(cits_e_term1 + cits_e_term2 - cits_e_term3 + cits_e_term4)

              if (avg_point==1) {

                cits_term1_long <- cits_term1_long + (akk^2)*((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term_long)
                its_term1_long  <- its_term1_long  + (akk^2)*((1/mtkk_its))*r2*(cits_theta_term + cits_e_term_long)

              } else if ((avg_point==2) & (iqkk==1)) {

                cits_term1_long <- cits_term1_long + ((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term_long)
                its_term1_long  <- its_term1_long  + ((1/mtkk_its))*r2*(cits_theta_term + cits_e_term_long)

              }

            }
          } # end type_cits==2

          # CITS - Discrete model

          if ((did_cits>1) & (type_cits==3) & (avg_point==1)) {

            rho_pre1  <- calc_rho(rho,6,1,bkk,0,0,tbarkk,0,ar1,tpp)
            rho_pre2  <- calc_rho(rho,7,1,bkk,0,0,tbarkk,0,ar1,tpp)
            rho_pp1   <- calc_rho(rho,8,1,bkk,skk,ntimep,tbarkk,0,ar1,tpp)

            # CITS, Cross-section

            cits_theta_term1 <- (pp_diffkk^2)*((1/ssqtkk)+((bkk-1)*bkk*rho_pre1/(ssqtkk^2)))
            cits_theta_term2 <- 2*pp_diffkk*bkk*rho_pre2 / ssqtkk
            cits_theta_term3 <- 2*pp_diffkk*bkk*rho_pp1 / ssqtkk
            cits_theta_term  <- icc*(cits_theta_term1 + cits_theta_term2 - cits_theta_term3)

            cits_e_term4 <- (pp_diffkk^2) / ssqtkk
            cits_e_term  <- ((1-icc)/nsamp)*(cits_e_term4)

            cits_term1 <- cits_term1 + (akk^2)*((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term)
            its_term1  <- its_term1  + (akk^2)*((1/mtkk_its))*r2*(cits_theta_term + cits_e_term)

            # CITS longitudinal design

            if (cross_sec==2) {

              phi_pre1  <- calc_rho(phi,6,1,bkk,0,0,tbarkk,0,ar1,tpp)
              phi_pre2  <- calc_rho(phi,7,1,bkk,0,0,tbarkk,0,ar1,tpp)
              phi_pp1   <- calc_rho(phi,8,1,bkk,skk,ntimep,tbarkk,0,ar1,tpp)

              cits_e_term1 <- (pp_diffkk^2)*((1/ssqtkk)+((bkk-1)*bkk*phi_pre1/(ssqtkk^2)))
              cits_e_term2 <- 2*pp_diffkk*bkk*phi_pre2 / ssqtkk
              cits_e_term3 <- 2*pp_diffkk*bkk*phi_pp1 / ssqtkk
              cits_e_term_long  <- ((1-icc)/nsamp)*(cits_e_term1 + cits_e_term2 - cits_e_term3)

              cits_term1_long <- cits_term1_long + (akk^2)*((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term_long)
              its_term1_long  <- its_term1_long  + (akk^2)*((1/mtkk_its))*r2*(cits_theta_term + cits_e_term_long)

            }
          }   # end (type_cits==3) & (avg_point==1)

          if ((did_cits>1) & (type_cits==3) & (avg_point==2) & (iqkk==1)) {

            rho_pre1  <- calc_rho(rho,6,1,bkk,0,0,tbarkk,0,ar1,tpp)
            rho_pre2  <- calc_rho(rho,7,1,bkk,0,0,tbarkk,0,ar1,tpp)
            rho_pp1   <- calc_rho(rho,8,1,bkk,qkk,qkk,tbarkk,0,ar1,tpp)

            # CITS, Cross-section

            diff_qkk  <- tpp[qkk]-tbarkk

            cits_theta_term1 <- (diff_qkk^2)*((1/ssqtkk)+((bkk-1)*bkk*rho_pre1/(ssqtkk^2)))
            cits_theta_term2 <- 2*diff_qkk*bkk*rho_pre2 / ssqtkk
            cits_theta_term3 <- 2*diff_qkk*bkk*rho_pp1 / ssqtkk
            cits_theta_term  <- icc*(cits_theta_term1 + cits_theta_term2 - cits_theta_term3)

            cits_e_term4 <- (diff_qkk^2) / ssqtkk
            cits_e_term  <- ((1-icc)/nsamp)*(cits_e_term4)

            #print(paste('cits_theta_term',cits_theta_term))
            #print(paste('cits_e_term',cits_e_term))

            cits_term1 <- cits_term1 + ((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term)
            its_term1  <- its_term1  + ((1/mtkk_its))*r2*(cits_theta_term + cits_e_term)

            # CITS longitudinal design

            if (cross_sec==2) {

              phi_pre1  <- calc_rho(phi,6,1,bkk,0,0,tbarkk,0,ar1,tpp)
              phi_pre2  <- calc_rho(phi,7,1,bkk,0,0,tbarkk,0,ar1,tpp)
              phi_pp1   <- calc_rho(phi,8,1,bkk,qkk,qkk,tbarkk,0,ar1,tpp)

              cits_e_term1 <- (diff_qkk^2)*((1/ssqtkk)+((bkk-1)*bkk*phi_pre1/(ssqtkk^2)))
              cits_e_term2 <- 2*diff_qkk*bkk*phi_pre2 / ssqtkk
              cits_e_term3 <- 2*diff_qkk*bkk*phi_pp1 / ssqtkk
              cits_e_term_long  <- ((1-icc)/nsamp)*(cits_e_term1 + cits_e_term2 - cits_e_term3)

              cits_term1_long <- cits_term1_long + ((1/mtkk)+(1/mckk))*r2*(cits_theta_term + cits_e_term_long)
              its_term1_long  <- its_term1_long  + ((1/mtkk_its))*r2*(cits_theta_term + cits_e_term_long)

            }
          }   # end (type_cits==3) & (avg_point==2)

        }   # end k

        #
        # CALCULATE TOTAL VARIANCE AND COVARIANCES
        #

        did_tot  <- NA
        cits_tot <- NA
        its_tot  <- NA

        sumak2 <- sumwk^2
        sumiq2 <- sumiq^2

        #if (avg_point==1) {
        #  sumk2 <- sumwk^2
        #} else if (avg_point==2) {
        #  sumk2 <- sumiq^2
        #}

        # Total DID Variances

        if (cross_sec==1) {

          if (avg_point==1) {
            did_tot <- deff_wgt*(1/sumak2)*term1
          } else if (avg_point==2) {
            did_tot <- deff_wgt*(1/sumiq2)*term1q
          }

        } else if (cross_sec==2) {

          if (avg_point==1) {
            did_tot <- deff_wgt*(1/sumak2)*term1_long
          } else if (avg_point==2) {
            did_tot <- deff_wgt*(1/sumiq2)*term1q_long
          }
        }

        #if (did_cits==1) {
          #print(paste('did_tot',did_tot))
        #}

        # Total CITS and ITS Covariances and Variances - Adds in the DID values

        if ((did_cits>1) & (type_cits==1)) {

          if (cross_sec==1) {

            if (avg_point==1) {

              cits_t1 <- term1 + cits_term1
              cits_tot <- deff_wgt*(1/sumak2)*cits_t1

              its_t1 <- term1_noc + its_term1
              its_tot <- deff_wgt*(1/sumak2)*its_t1

              #print(paste('term1',term1))
              #print(paste('cits_term1',cits_term1))
              #print(paste('cits_theta_term1',cits_theta_term1))

            } else if (avg_point==2) {

              cits_t1 <- term1 + cits_term1
              cits_tot <- deff_wgt*(1/sumiq2)*cits_t1

              its_t1 <- term1_noc + its_term1
              its_tot <- deff_wgt*(1/sumiq2)*its_t1

            }

          } else if (cross_sec==2) {

            if (avg_point==1) {

              cits_long_t1 <- term1_long + cits_term1_long
              cits_tot <- deff_wgt*(1/sumak2)*cits_long_t1

              its_long_t1 <- term1_long_noc + its_term1_long
              its_tot <- deff_wgt*(1/sumak2)*its_long_t1

              #print(paste('term1_long',term1_long))
              #print(paste('cits_term1_long',cits_term1_long))
              #print(paste('cits_theta_term1',cits_theta_term1))

            } else if (avg_point==2) {

              cits_long_t1 <- term1_long + cits_term1_long
              cits_tot <- deff_wgt*(1/sumiq2)*cits_long_t1

              its_long_t1 <- term1_long_noc + its_term1_long
              its_tot <- deff_wgt*(1/sumiq2)*its_long_t1

            }
          }

          #print(paste('cits_tot',cits_tot))
          #print(paste('its_tot',its_tot))

        }

        if ((did_cits>1) & (type_cits==2)) {

          if (cross_sec==1) {

            if (avg_point==1) {

              cits_t1 <- cits_term1
              cits_tot <- deff_wgt*(1/sumak2)*cits_t1

              its_t1 <- its_term1
              its_tot <- deff_wgt*(1/sumak2)*its_t1

            } else if (avg_point==2) {

              cits_t1 <- cits_term1
              cits_tot <- deff_wgt*(1/sumiq2)*cits_t1

              its_t1 <- its_term1
              its_tot <- deff_wgt*(1/sumiq2)*its_t1
            }

          } else if (cross_sec==2) {

            if (avg_point==1) {

              cits_long_t1 <- cits_term1_long
              cits_tot <- deff_wgt*(1/sumak2)*cits_long_t1

              its_long_t1 <- its_term1_long
              its_tot <- deff_wgt*(1/sumak2)*its_long_t1

            } else if (avg_point==2) {

              cits_long_t1 <- cits_term1_long
              cits_tot <- deff_wgt*(1/sumiq2)*cits_long_t1

              its_long_t1 <- its_term1_long
              its_tot <- deff_wgt*(1/sumiq2)*its_long_t1
            }
          }

          #print(paste('cits_tot',cits_tot))
          #print(paste('its_tot',its_tot))
        }

          if ((did_cits>1) & (type_cits==3)) {

            if (cross_sec==1) {

              if (avg_point==1) {

                cits_t1 <- term1 + cits_term1
                cits_tot <- deff_wgt*(1/sumak2)*cits_t1

                its_t1 <- term1_noc + its_term1
                its_tot <- deff_wgt*(1/sumak2)*its_t1

              } else if (avg_point==2) {

                #print(paste('term1q',term1q))
                #print(paste('cits_term1',cits_term1))

                cits_t1 <- term1q + cits_term1
                cits_tot <- deff_wgt*(1/sumiq2)*cits_t1

                its_t1 <- term1q_noc + its_term1
                its_tot <- deff_wgt*(1/sumiq2)*its_t1

              }

            } else if (cross_sec==2) {

              if (avg_point==1) {

                cits_long_t1 <- term1_long + cits_term1_long
                cits_tot <- deff_wgt*(1/sumak2)*cits_long_t1

                its_long_t1 <- term1_long_noc + its_term1_long
                its_tot <- deff_wgt*(1/sumak2)*its_long_t1

              } else if (avg_point==2) {

                cits_long_t1 <- term1q_long + cits_term1_long
                cits_tot <- deff_wgt*(1/sumiq2)*cits_long_t1

                its_long_t1 <- term1q_long_noc + its_term1_long
                its_tot <- deff_wgt*(1/sumiq2)*its_long_t1
              }
            }

            #print(paste('cits_tot',cits_tot))
            #print(paste('its_tot',its_tot))

          }

        #
        # CONDUCT SAMPLE SIZE CALCULATIONS FOR EACH POWER LEVEL
        #

        # Function to calculate required M for SAMP_SIZE==2

        calc_m <- function(var_term) {

          # Use Secant Method to Compute M

          if (var_term>0) {

            maxiter <- 25
            converge <- .000001

            # Compute initial function values

            m_opt   <- 30
            m_opt1  <- 20

            if (did_cits==1) {
              if (avg_point==1) {
                df  <- m_opt*ntimep - m_opt - ntimeg*ntimep - sumak
                df1 <- m_opt1*ntimep - m_opt1 - ntimeg*ntimep - sumak
              } else if (avg_point==2) {
                df  <- m_opt*ntimep - m_opt - sumiq*ntimep - sumiq
                df1 <- m_opt1*ntimep - m_opt1 - ntimeg*ntimep - sumiq
              }

            } else if ((did_cits==2) & (type_cits==1)) {
              if (avg_point==1) {
                df  <- m_opt*ntimep - 8*ntimeg
                df1 <- m_opt1*ntimep - 8*ntimeg
              } else if (avg_point==2) {
                df  <- m_opt*ntimep - 8*sumiq
                df1 <- m_opt1*ntimep - 8*sumiq
              }

            } else if ((did_cits==3) & (type_cits==1)) {
              if (avg_point==1) {
                df  <- m_opt*ntimep - 4*ntimeg
                df1 <- m_opt1*ntimep - 4*ntimeg
              } else if (avg_point==2) {
                df  <- m_opt*ntimep - 4*sumiq
                df1 <- m_opt1*ntimep - 4*sumiq
              }

            } else if ((did_cits==2) & (type_cits==2)) {
              if (avg_point==1) {
                df  <- m_opt*ntimep - 6*ntimeg
                df1 <- m_opt1*ntimep - 6*ntimeg
              } else if (avg_point==2) {
                df  <- m_opt*ntimep - 6*sumiq
                df1 <- m_opt1*ntimep - 6*sumiq
              }

            } else if ((did_cits==3) & (type_cits==2)) {
              if (avg_point==1) {
                df  <- m_opt*ntimep - 3*ntimeg
                df1 <- m_opt1*ntimep - 3*ntimeg
              } else if (avg_point==2) {
                df  <- m_opt*ntimep - 3*sumiq
                df1 <- m_opt1*ntimep - 3*sumiq
              }

            } else if ((did_cits==2) & (type_cits==3)) {
              if (avg_point==1) {
                df  <- m_opt*ntimep - 4*ntimeg - sumak
                df1 <- m_opt1*ntimep - 4*ntimeg - sumak
              } else if (avg_point==2) {
                df  <- m_opt*ntimep - 4*sumiq - sumiq
                df1 <- m_opt1*ntimep - 4*sumiq - sumiq
              }

            } else if ((did_cits==3) & (type_cits==3)) {
              if (avg_point==1) {
                df  <- m_opt*ntimep - 2*ntimeg - sumak
                df1 <- m_opt1*ntimep - 2*ntimeg - sumak
              } else if (avg_point==2) {
                df  <- m_opt*ntimep - 2*sumiq - sumiq
                df1 <- m_opt1*ntimep - 2*sumiq - sumiq
              }
            }

            if (df<0)  {df  <- 1}
            if (df1<0) {df1 <- 1}

            inv_alpha2 <- qt(alpha2,df)
            inv_power  <- qt(power,df)

            factor <- inv_alpha2 + inv_power

            func <- m_opt - ((factor^2)*var_term/(mde^2))

            inv_alpha2 <- qt(alpha2,df1)
            inv_power  <- qt(power,df1)

            factor <- inv_alpha2 + inv_power

            func1 <- m_opt1 - ((factor^2)*var_term/(mde^2))

            # Iterate

            iter <- 1
            while ((iter <= maxiter) & (max(abs(func))>converge)) {

              delta <- func*(m_opt - m_opt1)/(func - func1)
              m_opt1 <- m_opt
              m_opt  <- m_opt - delta

              func1 <- func

              if (did_cits==1) {
                if (avg_point==1) {
                  df  <- m_opt*ntimep - m_opt - ntimeg*ntimep - sumak
                } else if (avg_point==2) {
                  df  <- m_opt*ntimep - m_opt - sumiq*ntimep - sumiq
                }

              } else if ((did_cits==2) & (type_cits==1)) {
                if (avg_point==1) {
                  df  <- m_opt*ntimep - 8*ntimeg
                } else if (avg_point==2) {
                  df  <- m_opt*ntimep - 8*sumiq
                }

              } else if ((did_cits==3) & (type_cits==1)) {
                if (avg_point==1) {
                  df  <- m_opt*ntimep - 4*ntimeg
                } else if (avg_point==2) {
                  df  <- m_opt*ntimep - 4*sumiq
                }

              } else if ((did_cits==2) & (type_cits==2)) {
                if (avg_point==1) {
                  df  <- m_opt*ntimep - 6*ntimeg
                } else if (avg_point==2) {
                  df  <- m_opt*ntimep - 6*sumiq
                }

              } else if ((did_cits==3) & (type_cits==2)) {
                if (avg_point==1) {
                  df  <- m_opt*ntimep - 3*ntimeg
                } else if (avg_point==2) {
                  df  <- m_opt*ntimep - 3*sumiq
                }

              } else if ((did_cits==2) & (type_cits==3)) {
                if (avg_point==1) {
                  df  <- m_opt*ntimep - 4*ntimeg - sumak
                } else if (avg_point==2) {
                  df  <- m_opt*ntimep - 4*sumiq - sumiq
                }

              } else if ((did_cits==3) & (type_cits==3)) {
                if (avg_point==1) {
                  df  <- m_opt*ntimep - 2*ntimeg - sumak
                } else if (avg_point==2) {
                  df  <- m_opt*ntimep - 2*sumiq - sumiq
                }
              }

              if (df<0) {
                df <- 1
              }

              inv_alpha2 <- qt(alpha2,df)
              inv_power  <- qt(power,df)

              factor <- inv_alpha2 + inv_power

              func <- m_opt - ((factor^2)*var_term/(mde^2));

              iter <- iter + 1
            }

            if (iter <= maxiter) {
              m_optf <- round(m_opt)

            } else if (iter >= maxiter) {
              m_optf <- c("No Convergence")
            }

          } else if (var_term <= 0) {
            m_optf <- NA
            iter   <- NA
            func   <- NA
          }

          # Output a data frame

          m_optg <- data.frame(m_optf,iter,func)

          return(m_optg)
        }  # END of Optimization function

        #
        # CALCULATE POWER RESULTS FOR SAMPLE_SIZE = 1 or 2
        #

        ps <- round(power_s*100)
        pf <- round(power_f*100)
        pb <- round(power_by*100)

        countp <- 1
        for (powert in seq(ps,pf,pb)) {

          power <- powert/100

          mde_val  <- NA
          m_opt <- NA
          iter  <- NA
          func  <- NA

          # Do first for SAMP_SIZE = 1

          if (samp_size==1) {

            if (did_cits==1) {

              if (avg_point==1) {
                df  <- mclus*ntimep - mclus - ntimeg*ntimep - sumak
              } else if (avg_point==2) {
                df  <- mclusq*ntimep - mclusq - sumiq*ntimep - sumiq
              }

              inv_alpha2 <- qt(alpha2,df)
              inv_power  <- qt(power,df)

              factor <- inv_alpha2 + inv_power

              mde_val <- factor*(did_tot^.5)

            } else if ((did_cits==2) & (type_cits==1)) {

              if (avg_point==1) {
                df  <- mclus*ntimep - 8*ntimeg
              } else if (avg_point==2) {
                df  <- mclusq*ntimep - 8*sumiq
              }

              inv_alpha2 <- qt(alpha2,df)
              inv_power  <- qt(power,df)

              factor <- inv_alpha2 + inv_power

              mde_val <- factor*(cits_tot^.5)

            } else if ((did_cits==3) & (type_cits==1)) {

              if (avg_point==1) {
                df  <- mclus*ntimep - 4*ntimeg
              } else if (avg_point==2) {
                df  <- mclusq*ntimep - 4*sumiq
              }

              inv_alpha2 <- qt(alpha2,df)
              inv_power  <- qt(power,df)

              factor <- inv_alpha2 + inv_power

              mde_val <- factor*(its_tot^.5)

            } else if ((did_cits==2) & (type_cits==2)) {
              if (avg_point==1) {
                df  <- mclus*ntimep - 6*ntimeg
              } else if (avg_point==2) {
                df  <- mclusq*ntimep - 6*sumiq
              }

              inv_alpha2 <- qt(alpha2,df)
              inv_power  <- qt(power,df)

              factor <- inv_alpha2 + inv_power

              mde_val <- factor*(cits_tot^.5)

            } else if ((did_cits==3) & (type_cits==2)) {
              if (avg_point==1) {
                df  <- mclus*ntimep - 3*ntimeg
              } else if (avg_point==2) {
                df  <- mclusq*ntimep - 3*sumiq
              }

              inv_alpha2 <- qt(alpha2,df)
              inv_power  <- qt(power,df)

              factor <- inv_alpha2 + inv_power

              mde_val <- factor*(its_tot^.5)

            } else if ((did_cits==2) & (type_cits==3)) {
              if (avg_point==1) {
                df  <- mclus*ntimep - 4*ntimeg - sumak
              } else if (avg_point==2) {
                df  <- mclusq*ntimep - 4*sumiq - sumiq
              }

              inv_alpha2 <- qt(alpha2,df)
              inv_power  <- qt(power,df)

              factor <- inv_alpha2 + inv_power

              mde_val <- factor*(cits_tot^.5)

            } else if ((did_cits==3) & (type_cits==3)) {
              if (avg_point==1) {
                df  <- mclus*ntimep - 2*ntimeg - sumak
              } else if (avg_point==2) {
                df  <- mclusq*ntimep - 2*sumiq - sumiq
              }

              inv_alpha2 <- qt(alpha2,df)
              inv_power  <- qt(power,df)

              factor <- inv_alpha2 + inv_power

              mde_val <- factor*(its_tot^.5)

            }

          } else if (samp_size==2) {

            if (did_cits==1) {
              m_opt_out <- calc_m(did_tot)               # Calls optimization function

            } else if (did_cits==2) {
              m_opt_out <- calc_m(cits_tot)              # Calls optimization function

            } else if (did_cits==3) {
              m_opt_out <- calc_m(its_tot)               # Calls optimization function
            }

            m_opt <- m_opt_out$m_optf
            iter  <- m_opt_out$iter
            func  <- m_opt_out$func

          }

          # Write results to the data frame resgt and stack across power levels in resg

          if (samp_size==1) {
            power1   <- format(power, digits=2, nsmall=2)
            mde_val1 <- format(mde_val, digits=2, nsmall=2)
            resgt <- data.frame(power,mde_val,power1,mde_val1)
            #colnames(resgt) <- c("Power","MDE")
          }
          else if (samp_size==2) {
            power1 <- format(power, digits=2, nsmall=2)
            m_opt1 <- format(m_opt, digits=1, big.mark=",")
            resgt <- data.frame(power, m_opt, power1, m_opt1)
            #colnames(resgt) <- c("Power","Clusters", )
          }

          if (countp==1) {
            resg <- resgt
          } else if (countp>1) {
            resg <- rbind(resg,resgt)
          }

          countp <- countp + 1

        }  # End of Power Level loop

      } # End of if Crash=0

      #
      # WRITE RESULTS TO POW_RES WHICH IS A DATA FRAME
      #

      if (crash==0) {
        pow_res <- resg

      } else if (crash==1) {
        err_out <- data.frame(err_mess[1:nerr])
        colnames(err_out) <- c("Errors")
        pow_res <- err_out
      }

    }) # end pow_res reactive

    #
    # RETRIEVES THE CRASH VARIABLE NEEDED TO DEFINE THE TABLES AND PLOTS
    #

    crashv <- reactive({
      if (ncol(pow_res())==1) {
        crashv <- 1
      } else if (ncol(pow_res())>1)
      {
        crashv <- 0
      }
    })

    #
    # CREATES THE KABLEEXTRA CODE FOR THE TABLES TO RUN LATER
    # DEPENDING ON THE CRASH VARIABLE
    #

    # First get caption for the tables

    cap <- reactive({

      if (input$did_cits == "Difference-in-differences (DID)") {
        did_citsz <- 1
      } else if (input$did_cits == "Comparative interrupted time series (CITS)") {
        did_citsz <- 2
      } else if (input$did_cits == "Interrupted time series (ITS)") {
        did_citsz <- 3
      }

      if (input$type_cits == "Different pre- and post-period slopes") {
        type_citsz <- 1
      } else if (input$type_cits == "Common pre- and post-period slopes") {
        type_citsz <- 2
      } else if (input$type_cits == "Discrete post-period indicators (like DID)") {
        type_citsz <- 3
      }

      if (input$samp_size == "Calculate the minimum detectable effect size (MDE) for a given sample size") {
        samp_sizez <- 1
      } else if (input$samp_size == "Calculate the number of clusters required to achieve a target MDE value") {
        samp_sizez <- 2
      }

      mdez <- input$mde

      if (input$cross_sec=="Cross-sectional (different people over time)") {
        cross_secz <- 1
      } else if (input$cross_sec=="Longitudinal (same people over time)") {
        cross_secz <- 2
      }

      mtz    <- input$mt
      mcz    <- input$mc
      mclusz <- mtz + mcz

      if ((did_citsz==1) & (cross_secz==1)) {
        dsgn <- c("DID Cross-Sectional Design")
      } else if ((did_citsz==2) & (cross_secz==1) & (type_citsz==1)) {
        dsgn <- c("CITS Fully-Interacted Cross-Sectional Design")
      } else if ((did_citsz==3) & (cross_secz==1) & (type_citsz==1)) {
        dsgn <- c("ITS Fully-Interacted Cross-Sectional Design")
      } else if ((did_citsz==1) & (cross_secz==2)) {
        dsgn <- c("DID Longitudinal Design")
      } else if ((did_citsz==2) & (cross_secz==2) & (type_citsz==1)) {
        dsgn <- c("CITS Fully-Interacted Longitudinal Design")
      } else if ((did_citsz==3) & (cross_secz==2) & (type_citsz==1)) {
        dsgn <- c("ITS Fully-Interacted Longitudinal Design")
      } else if ((did_citsz==2) & (cross_secz==1) & (type_citsz==2)) {
        dsgn <- c("CITS Common-Slopes Cross-Sectional Design")
      } else if ((did_citsz==3) & (cross_secz==1) & (type_citsz==2)) {
        dsgn <- c("ITS Common-Slopes Cross-Sectional Design")
      } else if ((did_citsz==2) & (cross_secz==2) & (type_citsz==2)) {
        dsgn <- c("CITS Common-Slopes Longitudinal Design")
      } else if ((did_citsz==3) & (cross_secz==2) & (type_citsz==2)) {
        dsgn <- c("ITS Common-SLopes Longitudinal Design")
      } else if ((did_citsz==2) & (cross_secz==1) & (type_citsz==3)) {
        dsgn <- c("CITS Discrete Cross-Sectional Design")
      } else if ((did_citsz==3) & (cross_secz==1) & (type_citsz==3)) {
        dsgn <- c("ITS Discrete Cross-Sectional Design")
      } else if ((did_citsz==2) & (cross_secz==2) & (type_citsz==3)) {
        dsgn <- c("CITS Discrete Longitudinal Design")
      } else if ((did_citsz==3) & (cross_secz==2) & (type_citsz==3)) {
        dsgn <- c("ITS Discrete Longitudinal Design")
      }

      if (did_citsz==3) {
        mcz <- 0
      }

      mdez    <- format(mdez,digits=2, nsmall=2)
      mtz     <- format(mtz,digits=1)
      mcz     <- format(mcz,digits=1)

      if (crashv()==0) {
        if ((samp_sizez==1) & (did_citsz<=2)) {

      cap1 <- sprintf("MDE Results for the %s",dsgn)
      cap2 <- sprintf("with %s Treatment and %s Comparison Clusters",mtz,mcz)
      cap  <- data.frame(cap1,cap2)

      } else if ((samp_sizez==1) & (did_citsz==3)) {
      cap1 <- sprintf("MDE Results for the %s",dsgn)
      cap2 <- sprintf("with %s Treatment Clusters",mtz)
      cap  <- data.frame(cap1,cap2)

      } else if (samp_sizez==2) {
      cap1 <- sprintf("Required Sample Sizes for the %s",dsgn)
      cap2 <- sprintf("to Achieve an MDE of %s",mdez)
      cap  <- data.frame(cap1,cap2)
        }
      } else NULL
    })

    cap_tab <- reactive({
      if (crashv()==0) {
        cap_tab <- paste(cap()[1],cap()[2])
      } else NULL
    })

    cap_grph <- reactive({
      if (crashv()==0) {
        cap_grph <- paste(cap()[1],"\n",cap()[2])
      } else NULL
    })

    # Now for column names

    cn <- reactive({

      if (input$samp_size == "Calculate the minimum detectable effect size (MDE) for a given sample size") {
        ssz <- 1
      } else if (input$samp_size == "Calculate the number of clusters required to achieve a target MDE value") {
        ssz <- 2
      }

      if (ssz==1) {
        cn <- c("Power Level", "MDE Value")
      } else if (ssz==2) {
        cn <- c("Power Level", "Required Clusters")
      }
    })

    # Now for the footnote

    fn <- reactive({

      if (input$did_cits == "Difference-in-differences (DID)") {
        did_citsz <- 1
      } else if (input$did_cits == "Comparative interrupted time series (CITS)") {
        did_citsz <- 2
      } else if (input$did_cits == "Interrupted time series (ITS)") {
        did_citsz <- 3
      }

      if (input$type_cits == "Different pre- and post-period slopes") {
        type_citsz <- 1
      } else if (input$type_cits == "Common pre- and post-period slopes") {
        type_citsz <- 2
      } else if (input$type_cits == "Discrete post-period indicators (like DID)") {
        type_citsz <- 3
      }

      if (input$samp_size == "Calculate the minimum detectable effect size (MDE) for a given sample size") {
        samp_sizez <- 1
      } else if (input$samp_size == "Calculate the number of clusters required to achieve a target MDE value") {
        samp_sizez <- 2
      }

      mdez <- input$mde

      if (input$cross_sec=="Cross-sectional (different people over time)") {
        cross_secz <- 1
      } else if (input$cross_sec=="Longitudinal (same people over time)") {
        cross_secz <- 2
      }

      ntimepz <- input$ntimep

      tpz <- input$time_intervals

      ntimegz <- input$ntimeg

      skz <- input$start_time

      if (input$avg_point == "Average across all post-periods") {
        avg_pointz <- 1
      } else if (input$avg_point == "Single post-period time point") {
        avg_pointz <- 2
      }

      if (input$cal_expos == "At a specific time period") {
        cal_exposz <- 1
      } else if (input$cal_expos == "At a specific treatment exposure point") {
        cal_exposz <- 2
      }

      q_timez <- input$q_time
      l_timez <- input$l_time

      wtgz <- c("Treatment exposure")

      nsampz <- input$nsamp
      mtz    <- input$mt
      mtkz   <- input$mtk
      rtz    <- input$rt
      rtkz   <- input$rtk
      mcz    <- input$mc
      mckz   <- input$mck
      rckz   <- input$rck

      mclusz <- mtz + mcz

      ar1z <- input$ar1

      rhoz <- input$rho
      phiz <- input$phi
      iccz <- input$icc
      r2yxz <- input$r2yx
      r2txz <- input$r2tx
      deff_wgtz <- input$deff_wgt

      alphaz <- input$alpha

      if (input$two_tailed=="Two-tailed") {
        twotz <- c("Two-tailed")
      } else if (input$two_tailed=="One-tailed") {
        twotz <- c("One-tailed")
      }

      if ((did_citsz==1) & (cross_secz==1)) {
        dsgn <- c("DID Cross-Sectional Design")
      } else if ((did_citsz==2) & (cross_secz==1) & (type_citsz==1)) {
        dsgn <- c("CITS FUlly-Interacted Cross-Sectional Design")
      } else if ((did_citsz==3) & (cross_secz==1) & (type_citsz==1)) {
        dsgn <- c("ITS Fully-Interacted Cross-Sectional Design")
      } else if ((did_citsz==1) & (cross_secz==2)) {
        dsgn <- c("DID Longitudinal Design")
      } else if ((did_citsz==2) & (cross_secz==2) & (type_citsz==1)) {
        dsgn <- c("CITS Fully-Interacted Longitudinal Design")
      } else if ((did_citsz==3) & (cross_secz==2) & (type_citsz==1)) {
        dsgn <- c("ITS Fully-Interacted Longitudinal Design")
      } else if ((did_citsz==2) & (cross_secz==1) & (type_citsz==2)) {
        dsgn <- c("CITS Common-Slopes Cross-Sectional Design")
      } else if ((did_citsz==3) & (cross_secz==1) & (type_citsz==2)) {
        dsgn <- c("ITS Common-Slopes Cross-Sectional Design")
      } else if ((did_citsz==2) & (cross_secz==2) & (type_citsz==2)) {
        dsgn <- c("CITS Common-Slopes Longitudinal Design")
      } else if ((did_citsz==3) & (cross_secz==2) & (type_citsz==2)) {
        dsgn <- c("ITS Common-SLopes Longitudinal Design")
      } else if ((did_citsz==2) & (cross_secz==1) & (type_citsz==3)) {
        dsgn <- c("CITS Discrete Cross-Sectional Design")
      } else if ((did_citsz==3) & (cross_secz==1) & (type_citsz==3)) {
        dsgn <- c("ITS Discrete Cross-Sectional Design")
      } else if ((did_citsz==2) & (cross_secz==2) & (type_citsz==3)) {
        dsgn <- c("CITS Discrete Longitudinal Design")
      } else if ((did_citsz==3) & (cross_secz==2) & (type_citsz==3)) {
        dsgn <- c("ITS Discrete Longitudinal Design")
      }

      if (did_citsz==3) {
        rtz <- 1
        mcz <- 0
      }

      if (cross_secz==2) {
        phiz <- format(phiz,digits=2)
      } else if (cross_secz==1) {
        phiz <- NA
      }

      ar1z    <- format(ar1z,)
      rhoz    <- format(rhoz,digits=2)
      iccz    <- format(iccz,digits=2)
      deff_wgtz  <- format(deff_wgtz,digits=2)
      r2yxz   <- format(r2yxz,digits=2, nsmall=2)
      r2txz   <- format(r2txz,digits=2, nsmall=2)
      alphaz  <- format(alphaz,digits=2, nsmall=2)
      ntimepz <- format(ntimepz,digits=1)
      ntimegz <- format(ntimegz,digits=1)
      skz     <- format(skz,digits=1)
      mdez    <- format(mdez,digits=2, nsmall=2)
      nsampz  <- format(nsampz,digits=1)
      rtz     <- format(rtz,digits=2, nsmall=2)
      mtz     <- format(mtz,digits=1)
      mtkz    <- format(mtkz,digits=1)
      rtkz    <- format(rtkz,digits=2, nsmall = 2)
      mcz     <- format(mcz,digits=1)
      mckz    <- format(mckz,digits=1)
      rckz    <- format(rckz,digits=2, nsmall = 2)
      q_timez <- format(q_timez,digits=1)
      l_timez <- format(l_timez,digits=1)


      if (avg_pointz==1) {
        avgp <- c("Inputs: Average post-period treatment effect,")
      } else if ((avg_pointz==2) & (cal_exposz==1)) {
        avgp <- sprintf(c("Inputs: Treatment effect in post-period %s,"),q_timez)
      } else if ((avg_pointz==2) & (cal_exposz==2)) {
        avgp <- sprintf(c("Inputs: Treatment effect at exposure point %s,"),l_timez)
      }

      if (crashv()==0) {

        if (samp_sizez==1) {

          #txt1 <- c("Inputs: Periods = %s, Subjects = %s, Number of timing groups (TGs) = %s,")
          txt1 <- c("Periods = %s, Subjects = %s, Number of timing groups (TGs) = %s,")
          txt2 <- c("TG start times = %s, TG treatment samples = %s, TG comparison samples = %s,")
          txt3 <- c("Autocorrelation = %s, Rho = %s, Phi = %s, ICC = %s,")
#          txt4 <- c("Timing group weight = %s, Design effect due to various forms of weighting = %s, R2yx = %s, R2tx = %s, Alpha = %s,")
          txt4 <- c("Design effect due to weighting = %s, R2yx = %s, R2tx = %s, Alpha = %s,")
          txt5 <- c("%s test")
          txtg <- paste(avgp,txt1,txt2,txt3,txt4,txt5)

          fn <- sprintf(txtg,ntimepz,nsampz,ntimegz,paste(skz,collapse=" "),
                        paste(mtkz,collapse=" "), paste(mckz,collapse=" "), ar1z,
                        rhoz,phiz,iccz,deff_wgtz,r2yxz,r2txz,alphaz,twotz)
        } else if (samp_sizez==2) {

          txt1 <- c("Periods = %s, Subjects = %s, Treatment share = %s, Number of timing groups (TGs) = %s,")
          txt2 <- c("TG start times = %s, TG group treatment shares = %s, TG group comparison shares = %s,")
#          txt3 <- c("Autocorrelation = %s, Rho = %s, Phi = %s, ICC = %s, Timing group weight = %s,")
          txt3 <- c("Autocorrelation = %s, Rho = %s, Phi = %s, ICC = %s,")
          txt4 <- c("Design effect due to weighting = %s, R2yx = %s, R2tx = %s, Alpha = %s, %s test")
          txtg <- paste(avgp,txt1,txt2,txt3,txt4)

          fn <- sprintf(txtg,ntimepz,nsampz,rtz,ntimegz,paste(skz,collapse=" "),
                        paste(rtkz,collapse=" "), paste(rckz,collapse=" "), ar1z,
                        rhoz,phiz,iccz,deff_wgtz,r2yxz,r2txz,alphaz,twotz)
        } else NULL
      }
    })

    ktab <- reactive({

      if (crashv() == 0) {
      kbl(pow_res()[3:4],
          caption = paste("<center><strong>", cap_tab(), "<center><strong>"),
          align=c("c","c"),escape=FALSE,col.names=cn()) %>%
          kable_styling(full_width=FALSE,position="left",
          bootstrap_options = "striped") %>%
          column_spec(1, color ="#00aff5",bold=TRUE) %>%          #007bbc
          column_spec(2, color ="#00aff5",bold=TRUE) %>%
          footnote(general = fn(), fixed_small_size = TRUE)

      } else if (crashv() == 1)
      {
      kbl(pow_res()) %>%
          kable_styling(full_width=FALSE,position="left") %>%
          column_spec(1, color = "red", width = 5)
      }
    })

    #
    # CREATES THE GGPLOT2 CODE FOR THE PLOTS TO RUN LATER
    # DEPENDING ON THE CRASH AND SAMP_SIZE VARIABLE
    #

    grph <- reactive({

      if ((crashv() == 0) & (input$samp_size=="Calculate the minimum detectable effect size (MDE) for a given sample size")) {
        #print(pow_res())
        ggplot(
          pow_res()[1:2],                     # Reads the data
          aes(x = power,                      # aes sets the axes
              y = mde_val,
              group = 1)) +
          #geom_line()
          geom_line(color="#00aff5",size=1.2) +
          labs(x = cn()[1], y = cn()[2]) +
          theme(axis.title.x = element_text(face="bold")) +
          theme(axis.title.y = element_text(face="bold"))
      } else if ((crashv() == 0) & (input$samp_size=="Calculate the number of clusters required to achieve a target MDE value"))
      {
        ggplot(
          pow_res()[1:2],
          aes(x = power,
              y = m_opt,
              group = 1)) +
          geom_line()
#          geom_line(color="#00aff5",size=1.2) +
#          scale_y_continuous(label=comma) +
#          labs(x = cn()[1], y = cn()[2]) +
#          theme(axis.title.x = element_text(face="bold")) +
#          theme(axis.title.y = element_text(face="bold"))
      } else NULL

    })

    # This one is for downloading

    grph1 <- reactive({

      if ((crashv() == 0) & (input$samp_size=="Calculate the minimum detectable effect size (MDE) for a given sample size")) {
        ggplot(
          pow_res()[1:2],
          aes(x = power,
              y = mde_val,
              group = 1)) +
          geom_line()
#          geom_line(color="#00aff5",size=1.2) +
#          labs(title = cap_grph(), caption = fn(), x = cn()[1], y = cn()[2]) +
#          # Allows the footnote to wrap and adds space before it
#          theme(plot.caption=element_textbox_simple(hjust=0, size=10, margin=margin(10,0,0,0))) +
#          theme(plot.title = element_text(hjust = 0.5)) +        # This centers the title
#          theme(axis.title.x = element_text(face="bold")) +
#          theme(axis.title.y = element_text(face="bold"))

      } else if ((crashv() == 0) & (input$samp_size=="Calculate the number of clusters required to achieve a target MDE value"))
      {
#        ggplot(
#          pow_res()[1:2],
#          aes(x = power,
#              y = m_opt,
#              group = 1)) #+
#          geom_line(color="#00aff5",size=1.2) +
#          scale_y_continuous(label=comma) +
#          labs(title = cap(), caption = fn(), x = cn()[1], y = cn()[2]) +
#          theme(plot.caption=element_textbox_simple(hjust=0, size=10, margin=margin(10,0,0,0))) +
#          theme(plot.title = element_text(hjust = 0.5)) +
#          theme(axis.title.x = element_text(face="bold")) +
#          theme(axis.title.y = element_text(face="bold"))
      } else NULL

    })

  #
  # CREATE SUBMIT BUTTON AND MAKE TABLES IN PANEL 1 AND PLOTS IN PANEL 2
  #

    # Need to reset the output file name for those who leave it blank or the download produces an error
    filen <- reactive({
      ifelse(input$file_name != "", input$file_name, "Power_Panel graph")
    })

  observe({

    if (input$submit_button > 0) {
      #output$vv <- renderText({ tt() })
      #output$table <- renderTable({ crashv() })

      output$tableset <- renderText({ ktab() })

      if (is.null(grph())) {
        output$plotset <- renderPlot(NULL)
      } else
      {
        output$plotset <- renderPlot({ grph() })

        #Download table to file
        output$down <- downloadHandler(
          filename = function() {
            paste(filen(), input$png_pdf, sep=".")
          },
          content = function(file) {
            plot_function <- match.fun(input$png_pdf)

            plot_function(file)
            print( grph1() )
            #save_kable(ktab())
            dev.off()
         }
        )
      }
    }
  })

})  # End server

suppressMessages(
  suppressWarnings(
    runApp(shinyApp(ui = ui, server = server))
))




