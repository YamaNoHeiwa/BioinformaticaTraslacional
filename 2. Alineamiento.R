# Script para cargar el archivo sam

# Instalar Rsamtools si aún no está instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsamtools")

library(Rsamtools)

# Cargar archivo y convertirlo a formato BAM
#asBam("alineacion_SRR494102.sam", destination="alineacion_SRR494102.bam")
bamFile <- "~/Documentos/R/TFM/Beta/alineacion_SRR494102.bam.bam"
indexBam(bamFile)
alineaciones <- scanBam(bamFile)
alineaciones

# Contar cuantas alineaciones hay en un intervalo específico
param <- ScanBamParam(which=GRanges("chr1", IRanges(1, 1e6)))
conteos <- countBam(bamFile, param=param)
conteos

# Llamado de variantes

library(VariantAnnotation)
library(GenomicRanges)

# Leemos el archivo vcf con el genoma de referencia hg38
vcf <- readVcf("~/Documentos/R/TFM/Beta/variantes.vcf", "hg38")

# Filtrar variantes basado en calidad
vcf.filtrado <- vcf[qual(vcf) > 20]

head(vcf.filtrado)
metadata(vcf.filtrado)

# Extraer información específica del archivo de variantes filtrado
variantes <- data.frame(
  CHROM = seqnames(rowRanges(vcf.filtrado)),
  POS = start(rowRanges(vcf.filtrado)),
  REF = ref(vcf.filtrado),
  ALT = alt(vcf.filtrado),
  QUAL = qual(vcf.filtrado)
)
head(variantes)


####################################

# Anotaciones
# Lectura de archivos de BBDD

library(VariantAnnotation)

#dbsnp_vcf <- readVcf("~/Documentos/R/TFM/Beta/data/VCF/dbSNP/00-All_papu.vcf", "hg38")
clinvar_vcf <- readVcf("~/Documentos/R/TFM/Beta/data/VCF/ClinVAR/clinvar.vcf", "hg38")

alt_vals <- mapply(function(alt_list) {
  paste(as.character(alt_list), collapse = ",")
}, alt(clinvar_vcf))

# El campo 'FILTER' puede ser complejo. Extraemos y procesamos.
filter_vals <- mapply(function(filter_list) {
  if (is.na(filter_list)) {
    NA_character_
  } else {
    paste(as.character(filter_list), collapse = ";")
  }
}, filter(clinvar_vcf), SIMPLIFY = TRUE)

# Creamos el data frame para el objeto GRanges
clinvar_df <- data.frame(
  seqnames = seqnames(clinvar_vcf),
  start = start(clinvar_vcf),
  end = end(clinvar_vcf),
  strand = "*",
  ref = ref(clinvar_vcf),
  alt = alt_vals,
  filter = filter_vals,
  stringsAsFactors = FALSE
)

# Paso 3: Convertir el data frame en un objeto GRanges
clinvar_gr <- makeGRangesFromDataFrame(clinvar_df, keep.extra.columns = TRUE)


# Añadir los campos uno por uno, asegurándote de que cada uno

# Anotar con dbSNP
#snp_anno <- locateVariants(vcf.filtrado, dbsnp_vcf, AllVariants())

# Convertimos los archivos VCF a objetos GRanges
granges_gr <- GRanges(seqnames = seqnames(vcf.filtrado),
                       ranges = IRanges(start = start(vcf.filtrado), end = end(vcf.filtrado)),
                       strand = strand(vcf.filtrado),
                       REF = ref(vcf.filtrado),
                       ALT = alt(vcf.filtrado),
                       QUAL = qual(vcf.filtrado))

clinvar_gr <- GRanges(
  seqnames = seqnames(clinvar_vcf),
  ranges = ranges(clinvar_vcf),
  strand = rep("*", length(clinvar_vcf))
)

# Añadir las columnas de metadatos
mcols(clinvar_gr)$REF <- ref(clinvar_vcf)
mcols(clinvar_gr)$ALT <- alt(clinvar_vcf)
mcols(clinvar_gr)$QUAL <- as.numeric(as.character(qual(clinvar_vcf)))
mcols(clinvar_gr)$FILTER <- as.character(filter(clinvar_vcf))
mcols(clinvar_gr)$INFO <- info(clinvar_vcf)

# Ver el contenido del GRanges
head(granges_gr)
head(clinvar_gr)

overlaps <- findOverlaps(granges_vcf, clinvar_gr)



# Si tus datos VCF no tienen el prefijo 'chr', agrégalo
#seqlevels(granges_vcf) <- paste0("chr", seqlevels(granges_vcf))
granges_gr
clinvar_gr
# O si necesitas quitar el prefijo 'chr'
seqlevels(granges_vcf) <- gsub("chr", "", seqlevels(granges_vcf))

# Haz lo mismo para el objeto clinvar_gr para que coincidan
seqlevels(clinvar_gr) <- paste0("chr", seqlevels(clinvar_gr))
# O para quitar el prefijo 'chr'
seqlevels(clinvar_gr) <- gsub("chr", "", seqlevels(clinvar_gr))

overlaps <- findOverlaps(granges_vcf, clinvar_gr)

# Luego puedes extraer los datos de ClinVar basados en los solapamientos
clinvar_hits <- clinvar_gr[queryHits(overlaps)]

# Finalmente, combina esta información con tus datos de variantes
annotated_variants <- cbind(as.data.frame(granges_vcf[subjectHits(overlaps),]),
                            as.data.frame(clinvar_hits))

annotated_variants

# Anotar con ClinVar
clinvar_anno <- locateVariants(vcf.filtrado, clinvar_vcf, AllVariants())


