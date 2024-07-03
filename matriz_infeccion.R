# Cargar paquetes necesarios
library(ggplot2)

# Cargar el fichero con los datos
matriz <- read.csv("matriz_infeccion.csv", header = TRUE, row.names = 1)

# Cambiar datos nulos por "0"
matriz[is.na(matriz)] <- 0
matriz <- as.matrix(matriz)

# Convertir la matriz en un dataframe largo
df <- as.data.frame(as.table(matriz))

# Renombrar las columnas
colnames(df) <- c("Cepa", "Profago", "Valor")

# Crear el heatmap con ggplot2
figura <- ggplot(df, aes(x = Profago, y = Cepa, fill = factor(Valor))) +
  theme_light() +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("white", "#857b7b"), name = "Presencia") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5), 
        axis.text.y = element_text(size=8)) +
  labs(x = "Profago", y = "Cepa") +
  guides(fill = FALSE)

ggsave("matriz_infeccion.svg", plot = figura, width = 10, height = 6, bg = "transparent")
