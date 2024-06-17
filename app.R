if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()


library(VariantAnnotation)
library(biomaRt)
library(shiny)
library(dplyr)
library(httr)
library(jsonlite)
library(DT)

# Función para buscar información en ClinVar
fetch_clinvar_info <- function(variant_id) {
  esearch_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=", variant_id, "&retmode=json")
  esearch_response <- GET(esearch_url)
  
  if (status_code(esearch_response) == 200) {
    esearch_result <- content(esearch_response, "parsed", simplifyVector = TRUE)
    
    if (length(esearch_result$esearchresult$idlist) > 0) {
      clinvar_id <- esearch_result$esearchresult$idlist[[1]]
      esummary_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=", clinvar_id, "&retmode=json")
      esummary_response <- GET(esummary_url)
      
      if (status_code(esummary_response) == 200) {
        esummary_result <- content(esummary_response, "parsed", simplifyVector = TRUE)
        return(esummary_result$result[[as.character(clinvar_id)]])
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

ui <- fluidPage(
  titlePanel("Análisis de datos y genómica del cáncer"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Elige un archivo VCF", accept = c(".vcf", ".vcf.gz")),
      actionButton("loadData", "Cargar Datos")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Tabla de Variantes", DTOutput("variantsTable")),
        tabPanel("Información de la Variante", verbatimTextOutput("variantInfo")),
        tabPanel("Buscar en Google Scholar", uiOutput("externalLink")),
        tabPanel("Variantes Guardadas", DTOutput("savedVariantsTable"))
      )
    )
  )
)

server <- function(input, output, session) {
  variantes <- reactiveVal(NULL)
  saved_variants <- reactiveVal(data.frame())
  
  observeEvent(input$loadData, {
    req(input$file1)
    vcf_file <- input$file1$datapath
    genome <- "GRCh38"
    vcf <- readVcf(vcf_file, genome)
    
    variantes_df <- data.frame(
      chromosome = seqnames(rowRanges(vcf)),
      start = start(rowRanges(vcf)),
      end = end(rowRanges(vcf)),
      ref = ref(vcf),
      alt = sapply(alt(vcf), function(x) paste(x, collapse = "/")),
      id = rownames(info(vcf))
    )
    variantes(variantes_df)
  })
  
  output$variantsTable <- renderDT({
    req(variantes())
    datatable(variantes(), selection = "single", options = list(pageLength = 10))
  })
  
  observeEvent(input$variantsTable_rows_selected, {
    selected_row <- input$variantsTable_rows_selected
    if (length(selected_row) == 1) {
      variant_id <- variantes()[selected_row, "id"]
      info <- fetch_clinvar_info(variant_id)
      
      if (!is.null(info)) {
        info_text <- paste(
          "UID: ", info$uid, "\n",
          "Object Type: ", info$obj_type, "\n",
          "Accession: ", info$accession, "\n",
          "Accession Version: ", info$accession_version, "\n",
          "Title: ", info$title, "\n",
          "Variation Name: ", info$variation_name, "\n",
          "cDNA Change: ", info$cdna_change, "\n",
          "Variant Type: ", info$variant_type, "\n",
          "Germline Description: ", info$germline_classification$description, "\n",
          "Clinical Impact Description: ", info$clinical_impact_classification$description, "\n",
          "Oncogenicity Description: ", info$oncogenicity_classification$description, "\n",
          "Gene Symbol: ", info$genes$symbol, "\n",
          "Molecular Consequence: ", paste(info$molecular_consequence_list, collapse = ", "), "\n",
          "Protein Change: ", info$protein_change, "\n",
          sep = ""
        )
        output$variantInfo <- renderText({ info_text })
        
        # Generar el enlace a la búsqueda en Google Scholar para la variante
        google_scholar_url <- paste0("https://scholar.google.com/scholar?q=", "genetic variant ", variant_id)
        output$externalLink <- renderUI({
          tags$a(href = google_scholar_url, target = "_blank", "Buscar en Google Scholar")
        })
        
        # Agregar la variante seleccionada a la lista de variantes guardadas
        new_saved_variant <- variantes()[selected_row, ]
        saved_variants(rbind(saved_variants(), new_saved_variant))
      } else {
        output$variantInfo <- renderText({ "No se encontró información para esta variante." })
        output$externalLink <- renderUI({ NULL })
      }
    }
  })
  
  output$savedVariantsTable <- renderDT({
    datatable(saved_variants(), options = list(pageLength = 10))
  })
}

shinyApp(ui = ui, server = server)

