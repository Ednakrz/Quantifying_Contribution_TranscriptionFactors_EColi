'''
NAME
        Frecuencias_TF_porGO.R
VERSION
        1.0.0
AUTHOR
       Edna Karen Rivera Zagal

DESCRIPTION
    Este programa cuenta la frecuencia de cada factor de transcripción dominante en cada uno de renglones (GOs),
    de la matriz de contribución. Posteriormente se presentan los resultados en un histograma de la
    frecuencia de un factor de transcripción dominante en todos los GOs. Así mismo se crea un boxplot que
    muestra los outliers de los Factores de Transcripción dominantes en todos los GOs.

Input format
    file 1(contribution_matrix.txt)
    columns: TFs IDs (bnumber)
    indexes : GOs


Output format
    file 1 (Frecuencia_TFs_AllGos.tsv)
    (1) Dom_TF
    (2) Frecuencia

    Histograma (Cantidad de Gos, Tf dominante)
    Boxplot (Frecuecia, Gos)


USAGE
    R Frecuencias_TF_porGO.R


SEE ALSO
    none
GitHub link
'''

# Importar las librerias
library(ggplot2)
library(dplyr)
library(ggrepel)


# *********************Frecuencia TFs dominantentes ***************************#

# Procesar la columna Dom_TF para separar los TFs y contar su frecuencia
frecuencia_tfs <- Tfs_dominantes %>%
  # Separar los TFs en la columna Dom_TF
  tidyr::separate_rows(Dom_TF, sep = ", ") %>%
  # Contar la frecuencia de cada TF
  count(Dom_TF, name = "Frecuencia")

write.table(frecuencia_tfs, file = "Frecuencia_TFs_AllGos.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)



# **************************Histograma*****************************************#

# Crear el gráfico con etiquetas en las barras
ggplot(frecuencia_tfs, aes(x = Frecuencia, y = Dom_TF)) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Gráfico de barras
  geom_text(aes(label = Dom_TF), vjust = -0.5, color = "black", size = 4) +  # Etiquetas en las barras
  labs(
    title = "Frecuencia de TFs dominantes por GO",
    x = "Cantidad de GOs",
    y = "TF dominante"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar etiquetas del eje X


# Ordenar las barras de menor a mayor frecuencia
frecuencia_tfs <- frecuencia_tfs %>%
  arrange(Frecuencia)

# Crear el gráfico con etiquetas en las barras y ordenado
ggplot(frecuencia_tfs, aes(x = reorder(Dom_TF, Frecuencia), y = Frecuencia)) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Gráfico de barras
  geom_text(aes(label = Dom_TF), vjust = -0.5, color = "black", size = 4) +  # Etiquetas en las barras
  labs(
    title = "Frecuencia de TFs dominantes por GO",
    x = "TF dominante",
    y = "Frecuencia"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar etiquetas del eje X



# *********** Identificar outliers usando boxplot.stats ***********************#

outliers <- boxplot.stats(frecuencia_tfs$Frecuencia)$out

# Filtrar los TFs que son outliers
tfs_outliers <- frecuencia_tfs %>%
  filter(Frecuencia %in% outliers) %>%
  pull(Dom_TF)

# Crear el boxplot horizontal con etiquetas de outliers usando ggrepel
ggplot(frecuencia_tfs, aes(x = "", y = Frecuencia)) +
  geom_boxplot(fill = "steelblue", color = "black", outlier.color = "red", outlier.shape = 16) +  # Boxplot
  geom_text_repel(
    data = frecuencia_tfs %>% filter(Frecuencia %in% outliers),  # Filtrar outliers
    aes(label = Dom_TF),  # Mostrar el nombre del TF
    color = "red",  # Color de la etiqueta
    size = 4,  # Tamaño de la etiqueta
    direction = "x",  # Dirección de las etiquetas (horizontal)
    nudge_y = 0.5  # Ajustar la posición vertical de las etiquetas
  ) +
  coord_flip() +  # Cambiar la orientación a horizontal
  labs(
    title = "Distribución de frecuencias de TFs dominantes en GOs",
    y = "Frecuencia",
    x = "Gos"
  ) +
  theme_minimal()
