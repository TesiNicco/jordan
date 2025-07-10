# Libraries - suppress startup messages
suppressPackageStartupMessages({
  library(data.table)
  library(shiny)
  library(stringr)
  library(tools)
  library(shinyjs)
  library(processx)
})

# Define UI
ui <- fluidPage(
  clss = "container-fluid",
  useShinyjs(),

  tags$style(HTML("
    html, body {
      height: 100%;
      margin: 0;
      padding: 0;
    }

    .container-fluid {
      min-height: 100vh;
      display: flex;
      flex-direction: column;
      justify-content: space-between;
    }

    .spinner {
      border: 4px solid #f3f3f3;
      border-top: 4px solid #007BFF;
      border-radius: 50%;
      width: 20px;
      height: 20px;
      animation: spin 1s linear infinite;
      display: inline-block;
      margin-right: 8px;
    }

    @keyframes spin {
      0% { transform: rotate(0deg); }
      100% { transform: rotate(360deg); }
    }

    .image-container-foot {
      background-color: #ffffff;
      display: flex;
      justify-content: space-around;
      align-items: center;
      padding: 12px 0;
      border-top: 1px solid #ccc;
      margin-top: auto;
      width: 100%;
    }

    .image-container-foot img {
      height: 80px;
      max-width: 30%;
      object-fit: contain;
    }
  ")),

# Top logo or image banner
  tags$div(
    style = "width: 100%; overflow: hidden; margin-bottom: 20px;",
    tags$img(
      src = "images/dna.png",
      style = "width: 100%; height: auto; display: block;",
      alt = "Top Banner"
    )
  ),

  tags$head(tags$style(HTML('.valid-ok { color: darkgreen; font-weight: bold; font-size: 18px; }
    .valid-fail { color: red; font-weight: bold; font-size: 18px; }'))),

  titlePanel("Jordan: Command-line Wrapper"),

  tags$script(HTML("
    Shiny.addCustomMessageHandler('toggleRunButton', function(enabled) {
      const btn = document.getElementById('run_btn');
      if (enabled) {
        btn.removeAttribute('disabled');
      } else {
        btn.setAttribute('disabled', 'disabled');
      }
    });
  ")),

  sidebarLayout(
    sidebarPanel(
      textInput("data_path", "Input Data Path", value = ""),
      div(style = "margin-top: -10px; margin-bottom: 10px;", h6("Provide full path to a VCF, BCF, or PLINK-format file. For PLINK, give the path to the `.bed` file; other required files will be inferred.")),

      tags$hr(style = "margin-top: 5px; margin-bottom: 10px;"),

      textInput("snp_path", "SNP List Path", value = ""),
      div(style = "margin-top: -10px; margin-bottom: 10px;", h6("SNP list must include columns: CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE, and BETA (or OR).")),

      tags$hr(style = "margin-top: 5px; margin-bottom: 10px;"),

      textInput("out_path", "Output Path", value = ""),
      div(style = "margin-top: -10px; margin-bottom: 10px;", h6("Define the output directory. A new directory will be created if it doesn't exist.")),

      tags$hr(style = "margin-top: 5px; margin-bottom: 10px;"),

      checkboxInput("dosage", "Input is dosage (--dosage)", value = FALSE),
      checkboxInput("freq", "Calculate variant frequencies (--freq)", value = FALSE),
      checkboxInput("exclude", "Calculate PRS with and without APOE (--exclude)", value = FALSE),
      checkboxInput("multiple", "Multiple PLINK files should be used as input (--multiple)", value = FALSE),
      checkboxInput("plot", "Draw plots (--plot)", value = FALSE),
      # add conditional panel for plot
      conditionalPanel(
        condition = "input.plot == true",
        checkboxInput("exclude_na", "Exclude NAs from plots", value = FALSE),
      ),
      checkboxInput("directEffects", "Use variant effect-sizes as they are in the SNP list file instead of converting to risk (--directEffects)", value = FALSE),
      checkboxInput("keepDosage", "Keep dosages as additional file (--keepDosage)", value = FALSE),
      checkboxInput("maf", "Minor Allele Frequency filter (--maf)", value = FALSE),
      textInput("addWeight", "Additional weight for PRS (eg pathway-PRS or cell-PRS)", value = "", placeholder = "Column name in SNP list file to consider for double weighted-PRS"),
      checkboxInput("assoc_analysis", "Association analysis", value = FALSE),

      # add conditional panel for maf: input should be a number
      conditionalPanel(
        condition = "input.maf == true",
        numericInput("maf_value", "MAF Value", value = 0.01, min = 0, max = 1, step = 0.005),
        div(style = "margin-top: -10px; margin-bottom: 10px;", h6("Set the MAF threshold for filtering variants.")),
      ),

      # add conditional panel for association analysis
      conditionalPanel(
        condition = "input.assoc_analysis == true",
        textInput("pheno_path", "Phenotype Data Path", value = ""),
        selectInput("test_type", "Test Type", choices = c("PRS", "SNP", "Both"), selected = "PRS"),
        textInput("pheno_outcomes", "Outcomes (comma-separated)", value = ""),
        textInput("pheno_covariates", "Covariates (comma-separated)", value = ""),
        checkboxInput("sex_strata", "Sex-stratified analyses", value = FALSE),
        textInput("survival_var", "Variables for survival analysis", value = "", placeholder = "Variable name as reported in the phenotype file, eg SURVIVAL"),
        checkboxInput("tiles_analysis", "Tile-based analysis", value = FALSE),
        
        # add conditional panel for the tiles analysis
        conditionalPanel(
          condition = "input.tiles_analysis == true",
          textInput("n_tiles", "Number of tiles", value = "", placeholder = "Number of tiles to create, eg 10"),
          textInput("variable_tiles", "Variable to consider", value = "", placeholder = "Variable name as reported in the phenotype file, eg AD_CAT"),
          textInput("variable_tiles_reference", "Reference group", value = "", placeholder = "Reference group for the variable, eg 0"),
        ),

        checkboxInput("split_analysis", "Split-based analysis", value = FALSE),
        # add conditional panel for the split analysis
        conditionalPanel(
          condition = "input.split_analysis == true",
          textInput("variable_split", "Variable to consider", value = "", placeholder = "Variable name as reported in the phenotype file, eg MMSE"),
          textInput("variable_split_threshold", "Thresholds for the splits", value = "", placeholder = "Thresholds for splitting the variable, eg 1-10-20"),
        ),

        # add conditional panel for the association of tiles and splits
        conditionalPanel(
          condition = "input.tiles_analysis == true || input.split_analysis == true",
          checkboxInput("assoc_tiles_splits", "Association of tiles and splits", value = FALSE),
        ),
      ),
    ),

    mainPanel(
      conditionalPanel(
        condition = "output.data_type_loaded",
        tags$hr(),
        h4("Input Data Info"),
        div(style = "border: 1px solid #ccc; padding: 8px; background-color: #f8f8f8;", htmlOutput("data_type_info"))
      ),

      tags$hr(),
      h4("SNP List Info and Preview"),
      div(style = "border: 1px solid #ccc; padding: 8px; background-color: #f8f8f8;", htmlOutput("snp_info")),
      div(style = "max-height: 200px; overflow-y: scroll; border: 1px solid #ccc; padding: 5px;", tableOutput("snp_preview")),

      tags$hr(),
      h4("Output Path Info"),
      div(style = "border: 1px solid #ccc; padding: 8px; background-color: #f8f8f8;",
          htmlOutput("out_path_info")
      ),

      conditionalPanel(
        condition = "input.assoc_analysis == true",
        tags$hr(),
        h4("Phenotype Data Info and Preview"),
        div(style = "border: 1px solid #ccc; padding: 8px; background-color: #f8f8f8;", htmlOutput("pheno_info")),
        div(style = "max-height: 200px; overflow-y: scroll; border: 1px solid #ccc; padding: 5px;", tableOutput("pheno_preview"))
      ),

      tags$hr(),
      h4("Jordan Command Preview"),
      tags$div(
        style = "border: 1px solid #ccc; background-color: #f8f8f8; padding: 8px;",
        # Output text (with wrapping)
        verbatimTextOutput("cmd_preview"),
        # CSS for output formatting
        tags$style(HTML("
          #cmd_preview {
            white-space: pre-wrap;
            word-break: break-word;
            overflow-x: hidden;
            margin-bottom: 8px;
          }
        ")),
        # Copy button aligned bottom-left
        tags$div(
          style = "text-align: left;",
          tags$button(
            "Copy",
            id = "copy_btn",
            style = "font-size: 12px; padding: 4px 8px; border: 1px solid #ccc; background-color: #fff;"
          )
        ),
        # JavaScript for clipboard copy
        tags$script(HTML("
          document.getElementById('copy_btn').addEventListener('click', function() {
            var text = document.getElementById('cmd_preview').innerText;
            navigator.clipboard.writeText(text).then(function() {
              alert('Copied to clipboard!');
            }, function(err) {
              alert('Failed to copy text: ' + err);
            });
          });
        "))
      ),

      tags$hr(),
      uiOutput("run_button_ui"),
      uiOutput("spinner_container"),
      tags$div(style = "height: 10px;"),
      uiOutput("output_tabs")
      # Add vertical space

    )
  ),

  tags$div(class = "image-container-foot",
  tags$img(src = "images/amstUMC.jpg", alt = "ams"),
  tags$img(src = "images/tudelft1.png", alt = "tud"),
  tags$img(src = "images/github.png", alt = "gith"))
)

# Define server logic
server <- function(input, output, session) {
  # Directory to hold PDF plots (served via www/)
  plot_dir <- file.path("www", "tmp_plots")
  # Clear existing PDFs
  old_pdfs <- list.files(plot_dir, pattern = "\\.pdf$", full.names = TRUE)
  file.remove(old_pdfs)
  
  # Derive the path to the Jordan binary
  get_script_dir <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    match <- grep("--file=", cmdArgs)
    if (length(match) > 0) {
      # Called via Rscript app.R or Rscript /path/to/app.R
      return(normalizePath(dirname(sub("--file=", "", cmdArgs[match]))))
    } else {
      # Likely run inside RStudio or via source() — fallback to working dir
      return(normalizePath(getwd()))
    }
  }
  jordan_path <- str_replace_all(get_script_dir(),  '/bin/shiny_app', '')
  print(jordan_path)

  # Identify data type based on file extension
  data_type <- reactive({
      req(input$data_path)
      fname <- basename(input$data_path)
      ext <- tolower(tools::file_ext(fname))

      switch(
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
    })
  
  # Set default example paths
  jordan_dir <- file.path(jordan_path, "bin", "jordan.R")
  example_data_path <- file.path(jordan_path, "example_data", "example_data_plink.bed")
  example_snp_path  <- file.path(jordan_path, "example_data", "AD_snps.txt")
  example_out_path  <- file.path(jordan_path, "example_data", "test")
  example_pheno_path <- file.path(jordan_path, "example_data", "example_assoc_file.txt")
  example_pheno_outcomes <- "AD,MMSE,AD_CAT"
  example_pheno_covariates <- "SEX,PC1,PC2,PC3"

  # Update the UI fields with the computed values
  observe({
    updateTextInput(session, "data_path", value = example_data_path)
    updateTextInput(session, "snp_path", value = example_snp_path)
    updateTextInput(session, "out_path", value = example_out_path)
    updateTextInput(session, "pheno_path", value = example_pheno_path)
    updateTextInput(session, "pheno_outcomes", value = example_pheno_outcomes)
    updateTextInput(session, "pheno_covariates", value = example_pheno_covariates)
  })

  observe({
    valid <- file.exists(input$data_path) && file.exists(input$snp_path) && nzchar(input$out_path)
    if (isTRUE(input$assoc_analysis)) {
      valid <- valid && file.exists(input$pheno_path)
    }
    if (valid) {
      shinyjs::enable("run_btn")
    } else {
      shinyjs::disable("run_btn")
    }
  })

  output_file_path <- reactiveVal()
  log_file_path <- reactiveVal()
  prs_results_path <- reactiveVal()
  freq_results_path <- reactiveVal()
  assoc_prs_results_path <- reactiveVal()
  assoc_snp_results_path <- reactiveVal()
  log_lines <- character(0)
  proc <- NULL  # to store the running process
  proc_trigger <- reactiveVal(FALSE)  # this will signal polling start
  check_proc = reactiveTimer(500)
  proc_done = reactiveVal(FALSE)
  plot_index <- reactiveVal(1)

  output$run_button_ui <- renderUI({
    fluidRow(
      column(width = 2,
        actionButton("run_btn", "Run Jordan", class = "btn btn-primary")
      ),
      column(width = 1,
        uiOutput("spinner_container")  # spinner appears here
      )
    )
  })

  generate_cmd <- reactive({
    req(input$data_path, input$snp_path)

    flags <- c(
      if (input$dosage) "--dosage",
      if (input$freq) "--freq",
      if (input$exclude) "--exclude",
      if (input$multiple) "--multiple",
      if (input$addWeight != "") paste0("--addWeight ", input$addWeight),
      if (input$plot) "--plot",
      if (input$exclude_na) "exclude_NA",
      if (input$directEffects) "--directEffects",
      if (input$keepDosage) "--keepDosage",
      if (input$assoc_analysis) { paste0("--assoc ", input$test_type, " ", input$pheno_path, " --assoc-var ", input$pheno_outcomes, " --assoc-cov ", input$pheno_covariates)},
      if (input$maf) { paste0("--maf ", str_replace_all(input$maf_value, ",", ".")) },
      if (input$sex_strata) "--sex-strata",
      if (input$tiles_analysis) { paste0("--tiles '", input$n_tiles, ";", input$variable_tiles, ";", input$variable_tiles_reference, "'") },
      if (input$split_analysis) { paste0("--split '", input$variable_split, ";", input$variable_split_threshold, "'") },
      if (input$assoc_tiles_splits) "--assoc-split-tiles",
      if (input$survival_var != "") paste0("--survival-var ", input$survival_var)
      )

      input_jordan <- ifelse(
        data_type() %in% c("PLINK", "PLINK2"),
        sub("\\.(bed|bim|fam|pgen|pvar|psam)$", "", input$data_path),
        input$data_path
      )

    cmd <- paste(
      file.path(jordan_dir),
      "--genotype", input_jordan,
      "--snplist", input$snp_path,
      "--outname", input$out_path,
      paste(flags, collapse = " ")
    )

    return(list(cmd = cmd))
  })

  output$cmd_preview <- renderText({
    generate_cmd()$cmd
  })

  observeEvent(input$run_btn, {
    shinyjs::disable("run_btn")
    log_lines <<- character(0)

    output$spinner_container <- renderUI({
      tags$div(class = "spinner", title = "Running...", style = "margin-top: 6px;")
    })

    output$log_output <- renderText({ paste(tail(log_lines, 100), collapse = "\n") })
    cmd_info <- generate_cmd()
    jordan_cmd <- cmd_info$cmd

    # Start the process and assign to global proc
    proc <<- processx::process$new(
      command = "bash",
      args = c("-c", jordan_cmd),
      stdout = "|", stderr = "|"
    )

    # Reset output paths
    output_file_path(paste0(input$out_path, "/output.vcf"))
    log_file_path(paste0(input$out_path, "/output.log"))
    prs_results_path(file.path(input$out_path, "PRS_table.txt"))
    freq_results_path(file.path(input$out_path, "frequencies.txt"))
    assoc_prs_results_path(file.path(input$out_path, "association_results_PRS.txt"))
    assoc_snp_results_path(file.path(input$out_path, "association_results_single.txt"))

    # Trigger polling observer
    proc_done(FALSE)
    proc_trigger(TRUE)
    plot_index(1)  # Reset slideshow to first plot
  })

  plot_files <- reactive({
    req(proc_done())
    list.files(plot_dir, pattern = "\\.pdf$", full.names = FALSE)
  })

  output$output_tabs <- renderUI({
  tabs <- list(
    tabPanel("Log", verbatimTextOutput("log_output"))
  )

  if (proc_done()) {
    tabs <- append(tabs, list(
      tabPanel("PRS", 
        htmlOutput("results_output"),
        div(
          style = "max-height: 200px; overflow-y: scroll; border: 1px solid #ccc; padding: 5px;",
          tableOutput("prs_table")
        )
      )
    ))

    if (input$freq) {
      tabs <- append(tabs, list(
        tabPanel("Freq", 
          htmlOutput("freq_output"),
          div(
            style = "max-height: 200px; overflow-y: scroll; border: 1px solid #ccc; padding: 5px;",
            tableOutput("freq_table")
          )
        )
      ))
    }

    if (input$plot) {
      tabs <- append(tabs, list(
        tabPanel("Plot",
          uiOutput("plot_controls"),
          uiOutput("pdf_viewer")
        )
      ))
    }

    if (input$assoc_analysis) {
      if (input$test_type == "PRS") {
        tabs <- append(tabs, list(
          tabPanel("PRS-association", 
            htmlOutput("assoc_prs_output"),
            div(
              style = "max-height: 200px; overflow-y: scroll; border: 1px solid #ccc; padding: 5px;",
              tableOutput("assoc_prs_table")
            )
          )
        ))
      } else if (input$test_type == "SNP") {
        tabs <- append(tabs, list(
          tabPanel("SNP-association", 
            htmlOutput("assoc_snp_output"),
            div(
              style = "max-height: 200px; overflow-y: scroll; border: 1px solid #ccc; padding: 5px;",
              tableOutput("assoc_snp_table")
            )
          )
        ))
      } else if (input$test_type == "Both") {
        tabs <- append(tabs, list(
          tabPanel("PRS-association", 
            htmlOutput("assoc_prs_output"),
            div(
              style = "max-height: 200px; overflow-y: scroll; border: 1px solid #ccc; padding: 5px;",
              tableOutput("assoc_prs_table")
            )
          ),
          tabPanel("SNP-association", 
            htmlOutput("assoc_snp_output"),
            div(
              style = "max-height: 200px; overflow-y: scroll; border: 1px solid #ccc; padding: 5px;",
              tableOutput("assoc_snp_table")
            )
          )
        ))
      }
    }
  }

  do.call(tabsetPanel, c(list(id = "results_tabs"), tabs))
})

  observe({
  req(proc_trigger())
  invalidateLater(500)

  if (!is.null(proc) && proc$is_alive()) {
    proc$poll_io(100)

    if (proc$is_incomplete_output()) {
      new_output <- proc$read_output_lines(n = 1)
      log_lines <<- c(log_lines, new_output)
    }

    if (proc$is_incomplete_error()) {
      new_error <- proc$read_error_lines(n = 1)
      log_lines <<- c(log_lines, new_error)
    }

    output$log_output <- renderText({ paste(tail(log_lines, 100), collapse = "\n") })
  } else if (!is.null(proc)) {
  # Final read after process exits
  repeat {
    more_out <- proc$read_output_lines()
    more_err <- proc$read_error_lines()
    if (length(more_out) == 0 && length(more_err) == 0) break
    log_lines <<- c(log_lines, more_out, more_err)
  }

  # ✅ Now show complete final log (not just tail)
  output$log_output <- renderText({ paste(log_lines, collapse = "\n") })

  writeLines(log_lines, log_file_path())
  output$spinner_container <- renderUI({ NULL })
  shinyjs::enable("run_btn")

  # Update plots
  old_pdfs <- list.files(plot_dir, pattern = "\\.pdf$", full.names = TRUE)
  file.remove(old_pdfs)
  new_pdfs <- list.files(input$out_path, pattern = "\\.pdf$", full.names = TRUE)
  file.copy(new_pdfs, plot_dir, overwrite = TRUE)

  plot_index(1)
  proc_trigger(FALSE)
  proc_done(TRUE)
  proc <<- NULL
  }
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
    ext <- tolower(file_ext(fname))
    base <- file_path_sans_ext(fname)

    if (!file.exists(fpath)) {
      return(HTML(paste(fname, "-- Incorrect path. File does NOT exist", "<span class='valid-fail'>&#10007;</span>")))
    }

    dtype <- data_type()

    valid <- if (dtype %in% c("VCF", "BCF", "PLINK", "PLINK2"))
      "<span class='valid-ok'>&#10003;</span>"
    else
      "<span class='valid-fail'>&#10007;</span>"

    msg <- if (!is.na(dtype)) {
      if (dtype %in% c("PLINK", "PLINK2")) {
        paste0(base, " -- file exists -- data type: ", dtype)
      } else {
        paste0(fname, " -- file exists -- data type: ", dtype)
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
      return(HTML(paste(fname, "-- Incorrect path -- file type: SNP list", "<span class='valid-fail'>&#10007;</span>")))
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
      need(file.exists(path), "SNP list file does not exist. Define the path to the SNP list file to view the preview.")
    )
    fread(path, header = TRUE, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)

  output$pheno_info <- renderUI({
    req(input$assoc_analysis, input$pheno_path)
    fpath <- input$pheno_path
    fname <- basename(fpath)
    exists <- file.exists(fpath)

    if (!exists) {
      return(HTML(paste(fname, "-- Incorrect path -- file type: phenotype data", "<span class='valid-fail'>&#10007;</span>")))
    }

    valid <- "<span class='valid-ok'>&#10003;</span>"
    HTML(paste(fname, "-- file exists -- file type: phenotype data", valid))
  })

  output$pheno_preview <- renderTable({
    req(input$assoc_analysis, input$pheno_path)
    path <- input$pheno_path
    validate(
      need(file.exists(path), "Phenotype data file does not exist. Define the path to the phenotype data file to view the preview.")
    )
    fread(path, header = TRUE, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)

  output$out_path_info <- renderUI({
    req(input$out_path)
    out_path <- normalizePath(input$out_path, winslash = "/", mustWork = FALSE)
    parent_dir <- dirname(out_path)

    exists <- dir.exists(out_path)
    parent_exists <- dir.exists(parent_dir)

    msg <- if (exists) {
      paste(out_path, "-- directory exists")
    } else if (parent_exists) {
      paste(out_path, "-- will be created")
    } else {
      paste(out_path, "-- Incorrect path. Parent directory missing")
    }

    icon <- if (exists || parent_exists) {
      "<span class='valid-ok'>&#10003;</span>"
    } else {
      "<span class='valid-fail'>&#10007;</span>"
    }

    HTML(paste(msg, icon))
  })

  output$results_output <- renderUI({
    if (!proc_done()) {
    return(HTML("<i>Jordan is running. Results will be displayed once the run has completed.</i>"))
  }
    prs_path <- prs_results_path()  
    if (!file.exists(prs_path)) {
      return(HTML("<span class='valid-fail'>&#10007;</span> <b>PRS_table.txt</b> not found in the output directory."))
    }
    HTML("<h5>PRS Table</h5>")
  })
  output$prs_table <- renderTable({
    req(proc_done())
    req(file.exists(prs_results_path()))
    fread(prs_results_path())
  }, striped = TRUE, bordered = TRUE)

  output$freq_output <- renderUI({
    invalidateLater(1000, session)
    if (!proc_done()) {
    return(HTML("<i>Jordan is running. Results will be displayed once the run has completed.</i>"))
  }
    freq_path <- freq_results_path()
    if (!file.exists(freq_path)) {
      return(HTML("<span class='valid-fail'>&#10007;</span> <b>frequencies.txt</b> not found in the output directory."))
    }
    HTML("<h5>Frequencies Table</h5>")
  })
  output$freq_table <- renderTable({
    req(proc_done())
    req(file.exists(freq_results_path()))
    fread(freq_results_path())
  }, striped = TRUE, bordered = TRUE)

  output$assoc_prs_output <- renderUI({
    invalidateLater(1000, session)
    if (!proc_done()) {
    return(HTML("<i>Jordan is running. Results will be displayed once the run has completed.</i>"))
  }
    assoc_prs_path <- assoc_prs_results_path()
    if (!file.exists(assoc_prs_path)) {
      return(HTML("<span class='valid-fail'>&#10007;</span> <b>association_results_PRS.txt</b> not found in the output directory."))
    }
    HTML("<h5>Association Results (PRS)</h5>")
  })
  output$assoc_prs_table <- renderTable({
    req(proc_done())
    req(file.exists(assoc_prs_results_path()))
    fread(assoc_prs_results_path())
  }, striped = TRUE, bordered = TRUE)

  output$assoc_snp_output <- renderUI({
    invalidateLater(1000, session)
    if (!proc_done()) {
    return(HTML("<i>Jordan is running. Results will be displayed once the run has completed.</i>"))
  }
    assoc_snp_path <- assoc_snp_results_path()
    if (!file.exists(assoc_snp_path)) {
      return(HTML("<span class='valid-fail'>&#10007;</span> <b>association_results_single.txt</b> not found in the output directory."))
    }
    HTML("<h5>Association Results (SNP)</h5>")
  })
  output$assoc_snp_table <- renderTable({
    req(proc_done())
    req(file.exists(assoc_snp_results_path()))
    fread(assoc_snp_results_path())
  }, striped = TRUE, bordered = TRUE)

  output$plot_controls <- renderUI({
    files <- plot_files()
    if (!proc_done()) {
    return(HTML("<i>Jordan is running. Results will be displayed once the run has completed.</i>"))
  }
    if (length(files) == 0) return("No plots available.")

    tagList(
      actionButton("prev_plot", "Previous"),
      span(paste("Plot", plot_index(), "of", length(files))),
      actionButton("next_plot", "Next")
    )
  })

  output$pdf_viewer <- renderUI({
  if (!proc_done()) {
    return(HTML("<i>Jordan is running. Plots will be displayed once available.</i>"))
  }

  files <- plot_files()
  if (length(files) == 0) return(NULL)
  idx <- plot_index()
  if (idx > length(files)) idx <- length(files)
  file_to_show <- files[[idx]]
  timestamp <- as.numeric(Sys.time())  # cache buster

  tagList(
    tags$embed(
      src = file.path("tmp_plots", file_to_show),
      type = "application/pdf",
      width = "100%",
      height = "800px"
    ), tags$div(style = "height: 30px;"))
})

  observeEvent(input$prev_plot, {
    if (plot_index() > 1) plot_index(plot_index() - 1)
  })

  observeEvent(input$next_plot, {
    if (plot_index() < length(plot_files())) plot_index(plot_index() + 1)
  })
}

shinyApp(ui, server)