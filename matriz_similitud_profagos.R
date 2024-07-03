# Cargar paquetes necesarios
library(ggplot2)
library(readr)
library(pheatmap)

# Importar la matriz de ANI desde el fichero CSV 
matriz <- read.csv("matriz_similitud.csv", row.names = 1)
order_list <- readLines("orden.txt")

# Ordenar la matriz y convertirla en numérica
matriz_orden <- matriz[order_list, order_list]
matriz_numerica <- apply(matriz_orden, 2, function(x) gsub(",", ".", x))
matriz_numerica <- apply(matriz_numerica, 2, as.numeric)
matriz_numerica[is.na(matriz_numerica)] <- 0

# Crear un mapa de calor 
custom_colors <- colorRampPalette(colors = c("white", "lightblue1", "navy"))(100000)
figura <- pheatmap(matriz_numerica, 
         cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = FALSE, 
         fontsize = 6,
         color = custom_colors,
         angle_col = 90,
         labels_row = order_list)

# Guardar la figura en un fichero svg
ggsave("matriz_similitud.svg", plot = figura, width = 20, height = 20, bg = "transparent")



