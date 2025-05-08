# Libraries
library(data.table)
library(shiny)
library(stringr)
library(tools)

# Define UI for the application
ui <- fluidPage(
  tags$head(tags$style(HTML('
    hr {
      border-top: 1px solid #999 !important;
    }
    .valid-ok {
      color: darkgreen;
      font-weight: bold;
      font-size: 18px;
    }
    .valid-fail {
      color: red;
      font-weight: bold;
      font-size: 18px;
    }
  '))),

  titlePanel("Jordan Application: Command-line Wrapper"),

  sidebarLayout(
    sidebarPanel(
      textInput("data_path", "Input Data Path", value = "/data/project/chr1.bed"),
      div(style = "margin-top: -10px; margin-bottom: 10px;",
          h6("Provide full path to a VCF, BCF, or PLINK-format file. For PLINK, give the path to the `.bed` file; other required files will be inferred.")
      ),

      tags$hr(style = "margin-top: 5px; margin-bottom: 10px;"),

      textInput("snp_path", "SNP List Path", value = "/data/project/snp_list.txt"),
      div(style = "margin-top: -10px; margin-bottom: 10px;",
          h6("SNP list must include columns: CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE, and BETA (or OR).")
      ),

      tags$hr(style = "margin-top: 5px; margin-bottom: 10px;"),

      textInput("out_path", "Output Path", value = "/data/project/output"),
      div(style = "margin-top: -10px; margin-bottom: 10px;",
          h6("Define the output directory. A new directory will be created if it doesn't exist.")
      ),

      tags$hr(style = "margin-top: 5px; margin-bottom: 10px;"),

        checkboxInput("dosage", "Input is dosage (--dosage)", value = FALSE),
        checkboxInput("freq", "Calculate variant frequencies (--freq)", value = FALSE),
        checkboxInput("exclude", "Calculate PRS with and without APOE (--exclude)", value = FALSE),
        checkboxInput("multiple", "Multiple PLINK files should be used as input (--multiple)", value = FALSE),
        checkboxInput("plot", "Draw density plot (--plot)", value = FALSE),
        checkboxInput("directEffects", "Use variant effect-sizes as they are in the SNP list file instead of converting to risk (--directEffects)", value = FALSE),
        checkboxInput("keepDosage", "Keep dosages as additional file (--keepDosage)", value = FALSE),

      actionButton("run_btn", "Run Jordan")
    ),

    mainPanel(
      h4("Jordan Command Preview"),
      div(style = "border: 1px solid #ccc; background-color: #f8f8f8; padding: 8px;",
        verbatimTextOutput("cmd_preview")
    ),

      conditionalPanel(
        condition = "output.data_type_loaded",
        tags$hr(),
        h4("Input Data Info"),
        div(style = "border: 1px solid #ccc; padding: 8px; background-color: #f8f8f8;",
            htmlOutput("data_type_info")
        )
      ),

      tags$hr(),
      h4("SNP List Info and Preview"),
      div(style = "border: 1px solid #ccc; padding: 8px; background-color: #f8f8f8;",
          htmlOutput("snp_info")
      ),
      div(style = "max-height: 200px; overflow-y: scroll; border: 1px solid #ccc; padding: 5px;",
          tableOutput("snp_preview")
      ),

      tags$hr(),
      h4("Run Log"),
      verbatimTextOutput("log_output"),
      downloadButton("download_vcf", "Download VCF Output"),
      downloadButton("download_log", "Download Log")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  args <- commandArgs()
  script_path <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
  jordan_dir <- str_replace_all(script_path, '/shiny_app', '')

  output_file_path <- reactiveVal()
  log_file_path <- reactiveVal()

  generate_cmd <- reactive({
    req(input$data_path, input$snp_path)

    flags <- c(
      if (input$dosage) "--dosage",
      if (input$freq) "--freq",
      if (input$exclude) "--exclude",
      if (input$multiple) "--multiple"
    )

    out_prefix <- file.path(tempdir(), "jordan_output")

    cmd <- paste(
      file.path(jordan_dir, "jordan"),
      "--vcf", input$data_path,
      "--bed", input$snp_path,
      "--prefix", out_prefix,
      paste(flags, collapse = " ")
    )

    return(list(cmd = cmd, prefix = out_prefix))
  })

  output$cmd_preview <- renderText({
    generate_cmd()$cmd
  })

  observeEvent(input$run_btn, {
    cmd_info <- generate_cmd()
    log_path <- paste0(cmd_info$prefix, ".log")

    result <- system(cmd_info$cmd, intern = TRUE)
    writeLines(result, log_path)

    output$log_output <- renderText(paste(result, collapse = "\n"))

    output_file_path(paste0(cmd_info$prefix, ".vcf"))
    log_file_path(log_path)
  })

  output$download_vcf <- downloadHandler(
    filename = function() "output.vcf",
    content = function(file) file.copy(output_file_path(), file)
  )

  output$download_log <- downloadHandler(
    filename = function() "output.log",
    content = function(file) file.copy(log_file_path(), file)
  )

  output$data_type_info <- renderUI({
    req(input$data_path)
    fpath <- input$data_path
    fname <- basename(fpath)
    exists <- file.exists(fpath)

    if (!exists) {
      showNotification("Input data file does not exist.", type = "error", duration = 5)
      return(HTML(paste(fname, "-- incorrect path -- data type: Unknown", "<span class='valid-fail'>&#10007;</span>")))
    }

    ext <- tolower(file_ext(fname))
    base <- file_path_sans_ext(fname)

    type <- switch(
      ext,
      "vcf" = "VCF",
      "gz" = {
        if (grepl("\\.vcf\\.gz$", fname, ignore.case = TRUE)) "VCF"
        else if (grepl("\\.bcf\\.gz$", fname, ignore.case = TRUE)) "BCF"
        else NA
      },
      "bcf" = "BCF",
      "bed" = "PLINK",
      "bim" = "PLINK",
      "fam" = "PLINK",
      "pgen" = "PLINK2",
      "pvar" = "PLINK2",
      "psam" = "PLINK2",
      NA
    )

    valid <- if (type %in% c("VCF", "BCF", "PLINK", "PLINK2"))
      "<span class='valid-ok'>&#10003;</span>" else "<span class='valid-fail'>&#10007;</span>"

    msg <- if (!is.na(type)) {
      if (type %in% c("PLINK", "PLINK2")) {
        paste0(base, " -- file exists -- data type: ", type)
      } else {
        paste0(fname, " -- file exists -- data type: ", type)
      }
    } else {
      paste0(fname, " -- file exists -- data type: Unknown (consider changing this)")
    }

    HTML(paste(msg, valid))
  })

  output$data_type_loaded <- reactive({
    !is.null(input$data_path) && input$data_path != ""
  })
  outputOptions(output, "data_type_loaded", suspendWhenHidden = FALSE)

  output$snp_info <- renderUI({
    req(input$snp_path)
    fpath <- input$snp_path
    fname <- basename(fpath)
    exists <- file.exists(fpath)

    if (!exists) {
      return(HTML(paste(fname, "-- incorrect path -- file type: SNP list", "<span class='valid-fail'>&#10007;</span>")))
    }

    cols <- tryCatch({
      colnames(fread(fpath, header = TRUE, nrows = 1, check.names = FALSE))
    }, error = function(e) NULL)

    required <- c("CHROM", "POS", "EFFECT_ALLELE", "OTHER_ALLELE")
    ok <- all(toupper(required) %in% toupper(cols)) &&
      any(toupper("BETA") %in% toupper(cols) || toupper("OR") %in% toupper(cols))

    valid <- if (ok) "<span class='valid-ok'>&#10003;</span>" else "<span class='valid-fail'>&#10007;</span>"
    HTML(paste(fname, "-- file exists -- file type: SNP list", valid))
  })

  output$snp_loaded <- reactive({
    file.exists(input$snp_path)
  })
  outputOptions(output, "snp_loaded", suspendWhenHidden = FALSE)

  output$snp_preview <- renderTable({
    req(input$snp_path)
    path <- input$snp_path
    validate(
      need(file.exists(path), "SNP list file does not exist.")
    )
    fread(path, header = TRUE, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)
}

shinyApp(ui, server)