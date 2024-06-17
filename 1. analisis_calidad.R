
comprobarIntegridad = function(archivo){
  archivo_fastq <- archivo
  numero_de_lineas <- 1000
  lineas <- readLines(archivo_fastq, n = numero_de_lineas)
  
  # Verifica que el archivo tenga un número de líneas múltiplo de 4
  if (length(lineas) %% 4 == 0) {
    mensaje <- "El archivo tiene un formato preliminar correcto."
  } else {
    mensaje <- "El archivo NO tiene un formato preliminar correcto. Verifica el archivo."
  }
  
  # Verifica que las secuencias y las calidades tengan la misma longitud
  for (i in seq(1, length(lineas), by = 4)) {
    if (nchar(lineas[i+1]) != nchar(lineas[i+3])) {
      mensaje <- paste(mensaje, "\nDiscrepancia encontrada en las longitudes de secuencia y calidad en la línea", i+1)
      break
    }
  }
  
  print(mensaje)
}

# Control de calidad del ADN con Rqc
#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("Rqc")

library(Rqc)
# recortar líneas del FASTQ
# head -n 10000 SRR7108454.fastq > SRR7108454_10000.fastq

# Análisis del ADN explota sino tiene suficiente RAM la máquina que lo ejecuta
# por eso, he recortado el fichero fastq

setwd("~/Documentos/R/TFM/Beta")
archivo = "SRR494093.fastq"
archivo2 = "SRR7108454.fastq"
# Comprobamos el formato
comprobarIntegridad(archivo)

rqc(pattern = archivo)




