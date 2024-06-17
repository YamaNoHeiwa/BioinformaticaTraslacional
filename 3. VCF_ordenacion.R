
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
# BiocManager::install("VariantAnnotation")

library(VariantAnnotation)
library(biomaRt)
library(jsonlite)
library(dplyr)
library(httr)

# Define el path al archivo VCF y el genoma de referencia, en mi caso human chr1
# setwd()
setwd("~/Documentos/R/TFM/Beta/data/VCF/ClinVAR")

# vcf_file <- "data/VCF/fullgenome.vcf"
# URL https://ftp.ensembl.org/pub/release-112/variation/gvf/homo_sapiens/
# vcf_file = "clinvar.vcf"
# vcf_file
# genome <- "GRCh38" # release 112
# 
# # Leer el archivo VCF
# vcf <- readVcf(vcf_file, genome)
# 
# head(vcf)
# Extraer información básica de las variantes
# start: posición inicial de la variante en el cromosoma especificado
# end: el final (con un solo nucleótido suelen ser iguales)
# ref: secuencia de referencia, los nuclótidos en el genoma de referencia
# alt: secuencias alternativas para la variante, secuencias que reemplazan
# la secuencia de referencia en las muestras. Múltiples alternativas se separan
# con un "/".

# NO VOLVER A CARGAR, ME COSTÓ 8 HORAS
# variants_info <- data.frame(
#  chromosome = seqnames(rowRanges(vcf)),
# start = start(rowRanges(vcf)),
#  end = end(rowRanges(vcf)),
#  ref = ref(vcf),
#  alt = sapply(alt(vcf), function(x) paste(x, collapse = "/"))
# )
# head(variants_info)
# length(variants_info$chromosome) # 2.728.246 muestras
# variants_info # 2.728.046 variantes, tardó 8 horas en cargar
# Guardar el data.frame en un archivo CSV
#write.csv(variants_info, "variants_info.csv", row.names = FALSE)
# Leer el archivo CSV
variantes <- read.csv("variants_info.csv", stringsAsFactors = FALSE)
length(variantes)
# reducimos el tamaño por practicidad
variantes = variantes[1:200,]
variantes

# con stringsAsFactors = TRUE se lee como caracteres
#variantes <- read.csv("variants_info.csv", stringsAsFactors = FALSE)

create_variant_string <- function(chromosome, start, ref, alt) {
  variant_string <- paste(chromosome, start, start, ref, alt, sep = " ")
  return(variant_string)
}

variant_strings <- apply(variantes, 1, function(variant) {
  create_variant_string(variant["chromosome"], variant["start"], variant["ref"], variant["alt"])
})

url <- "https://rest.ensembl.org/vep/human/region"

response <- POST(url,
                 body = toJSON(list(variants = variant_strings)),
                 encode = "json",
                 add_headers("Content-Type" = "application/json", Accept = "application/json"))

if (http_status(response)$category == "Success") {
  result <- content(response, "parsed")
  print(result)
} else {
  print(paste("Error:", http_status(response)$message))
}

# Función 
extract_info <- function(variant) {
  most_severe_consequence <- ifelse(!is.null(variant$most_severe_consequence), variant$most_severe_consequence, NA)
  
  if (!is.null(variant$transcript_consequences) && length(variant$transcript_consequences) > 0) {
    transcript_consequences <- variant$transcript_consequences[[1]]
    polyphen_score <- ifelse(!is.null(transcript_consequences$polyphen_score), transcript_consequences$polyphen_score, NA)
    sift_score <- ifelse(!is.null(transcript_consequences$sift_score), transcript_consequences$sift_score, NA)
    impact <- ifelse(!is.null(transcript_consequences$impact), transcript_consequences$impact, NA)
  } else {
    polyphen_score <- NA
    sift_score <- NA
    impact <- NA
  }
  
  frequencies <- if (!is.null(variant$colocated_variants) && length(variant$colocated_variants) > 0) {
    freq <- variant$colocated_variants[[1]]$frequencies
    gnomadg <- ifelse(!is.null(freq$C$gnomadg), freq$C$gnomadg, NA)
    gnomade <- ifelse(!is.null(freq$C$gnomade), freq$C$gnomade, NA)
  } else {
    gnomadg <- NA
    gnomade <- NA
  }
  
  data.frame(
    input = variant$input,
    most_severe_consequence = most_severe_consequence,
    polyphen_score = polyphen_score,
    sift_score = sift_score,
    impact = impact,
    gnomadg = gnomadg,
    gnomade = gnomade,
    stringsAsFactors = FALSE
  )
}

# Extrae y mete la información en un data frame
results_df <- do.call(rbind, lapply(result, extract_info))

head(results_df)

# chromosome, start, end,ref alt
# . most_severe_consequence:
# . missense_variant: un cambio en una base de ADN.
# . synonymous_variant: no cambia el amino ácido codificado.
# . nonsense_variant: codon prematuro de parada
# . frameshift_variant: in/del
# . intron_variant: ocurre dentro de un intron
# . polyphen_score: predice el impacto de una sustitución de un aminoácido en
# la estructura y función de una proteína humana.
# . valores posibles: de 0 a , 0: probablemente benigno, 1: probablemente dañino. 
# . sift_score: un indicador de la herramienta SIFT, que predice si una sustitución
# de un aminoácido afecta a la función de la proteina basándose en la homología
# de secuencia y las propiedades físicas del aminoácido. 
# X < 0.05 predice que es deletéreo(mortífero), X > 0.05 predice que son tolerados
# - gnomADG: frecuencia alélica de la variante en gnomAD, 0.0009595 indica que
# la variante está en aprox. 0.096% de la población.
# - gnomADE: frecuencia alélica de la variante en gnomAD, 0.00104 indica que
# la variante está en aprox. 0.0104% del exoma de la población.


results_df$most_severe_consequence

# Función para ordenar las variantes por criterios
sort_variants <- function(data, criteria = "pathogenicity") {
  
  # patogenicidad
  if (criteria == "pathogenicity") {
    sorted_data <- data %>%
      arrange(most_severe_consequence, desc(polyphen_score), desc(sift_score))
    
  } # impacto 
  else if (criteria == "impact") {
    sorted_data <- data %>%
      arrange(desc(impact), desc(polyphen_score), desc(sift_score))
    
  } # frecuencia poblacional 
  else if (criteria == "frequency") {
    sorted_data <- data %>%
      arrange(desc(gnomadg), desc(gnomade))
    
  } else {
    stop("Invalid criteria. Please choose 'pathogenicity', 'impact', or 'frequency'.")
  }
  
  return(sorted_data)
}

# Usos
sorted_by_pathogenicity <- sort_variants(results_df, "pathogenicity")
sorted_by_impact <- sort_variants(results_df, "impact")
sorted_by_frequency <- sort_variants(results_df, "frequency")

# Mostrar las variantes ordenadas
head(sorted_by_pathogenicity)
head(sorted_by_impact)
head(sorted_by_frequency)






