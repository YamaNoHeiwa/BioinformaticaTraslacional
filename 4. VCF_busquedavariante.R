if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("VariantAnnotation")

library(VariantAnnotation)
library(biomaRt)
library(jsonlite)
library(dplyr)
library(httr)
library(knitr)

# Define el path al archivo VCF y el genoma de referencia, en mi caso human chr1
# setwd()
setwd("~/Documentos/R/TFM/Beta/data/VCF/ClinVAR")
vcf_file <- "clinvar_1000_variants.vcf"
genome <- "GRCh38"
vcf <- readVcf(vcf_file, genome)
head(vcf)

variantes <- data.frame(
  chromosome = seqnames(rowRanges(vcf)),
  start = start(rowRanges(vcf)),
  end = end(rowRanges(vcf)),
  ref = ref(vcf),
  alt = sapply(alt(vcf), function(x) paste(x, collapse = "/")),
  id = rownames(info(vcf))  # Utilizar rownames(info(vcf)) para obtener el ID de la variante)
)
head(variantes)
variantes
write.csv(variantes, "variants_info2.csv", row.names = FALSE)
variantes <- read.csv("variants_info2.csv", stringsAsFactors = FALSE)
variantes$id

# Función para buscar información en ClinVar
fetch_clinvar_info <- function(variant_id) {
  # Construir la URL para la solicitud esearch
  esearch_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=", variant_id, "&retmode=json")
  # Hacer la solicitud GET
  esearch_response <- GET(esearch_url)
  
  # Verificar el estado de la respuesta
  if (status_code(esearch_response) == 200) {
    esearch_result <- content(esearch_response, "parsed", simplifyVector = TRUE)
    
    if (length(esearch_result$esearchresult$idlist) > 0) {
      # Obtener el ID de ClinVar
      clinvar_id <- esearch_result$esearchresult$idlist[[1]]
      
      # Construir la URL para obtener el resumen de ClinVar
      esummary_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=", clinvar_id, "&retmode=json")
      
      # Hacer la solicitud GET para el resumen
      esummary_response <- GET(esummary_url)
      
      # Verificar el estado de la respuesta
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

# Buscamos en CLINVAR
#variant_id <- "2205837"
variantes$id
info <- fetch_clinvar_info(variant_id)
print(info)

# Creamos un dataframe con la información
info_df <- data.frame(
  uid = info$uid,
  obj_type = info$obj_type,
  accession = info$accession,
  accession_version = info$accession_version,
  title = info$title,
  measure_id = info$variation_set$measure_id,
  variation_name = info$variation_set$variation_name,
  cdna_change = info$variation_set$cdna_change,
  variant_type = info$variation_set$variant_type,
  canonical_spdi = info$variation_set$canonical_spdi,
  scv = info$supporting_submissions$scv,
  rcv = info$supporting_submissions$rcv,
  germline_description = info$germline_classification$description,
  germline_last_evaluated = info$germline_classification$last_evaluated,
  germline_review_status = info$germline_classification$review_status,
  germline_trait_name = info$germline_classification$trait_set$trait_name,
  clinical_impact_description = info$clinical_impact_classification$description,
  clinical_impact_last_evaluated = info$clinical_impact_classification$last_evaluated,
  clinical_impact_review_status = info$clinical_impact_classification$review_status,
  oncogenicity_description = info$oncogenicity_classification$description,
  oncogenicity_last_evaluated = info$oncogenicity_classification$last_evaluated,
  oncogenicity_review_status = info$oncogenicity_classification$review_status,
  record_status = info$record_status,
  gene_symbol = info$genes$symbol,
  gene_id = info$genes$geneid,
  gene_strand = info$genes$strand,
  gene_source = info$genes$source,
  molecular_consequence = info$molecular_consequence_list,
  protein_change = info$protein_change,
  stringsAsFactors = FALSE
)

info_df

# uid: identificador de ClinVar
# obj_type/variant_type: tipo de variante
# title: descripción de la variante
# germline_description: Clasificación germinal
# gene_id: id del gen de NCBI
# molecular_consequence: consecuencia molecular de la variante
# protein_change: cambio de aminoácido resultante



