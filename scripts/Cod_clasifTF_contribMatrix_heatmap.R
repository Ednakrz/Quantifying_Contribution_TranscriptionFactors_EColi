'''
NAME
        ClasifTF_contribMatrix_heatmap.R
VERSION
        1.0.0
AUTHOR
       Edna Karen Rivera Zagal

DESCRIPTION
    Código estructurado para el análisis de factores de transcripción de una matriz
    de contribución. El script está estructurado en dos partes. La primera está centrada
    en la clasificación de los Factores de transcripción que se ecuentran en la
    matriz de contribución. La clasificación está dada por:

      0. TFs de valor de cero
      1. TFs mayores que cero pero distintos a TFs dominantes u outliers
      2. TFs outliers
      3. TFs dominantes
      4. TFs dominantes y outliers

    En la segunda parte, una vez que se obtiene la clasificación de los factores de
    transcripción, se obtiene un heatmap de la clasificación y el dendograma de dicha
    clasificación.

Input format
    file 1(contribution_matrix.txt)
    columns: TFs IDs (bnumber)
    indexes : GOs


Output format
    file 1 (df_clasificacionTFs.tsv)
    columns: TFs IDs (bnumber)
    indexes : GOs
    con la clasificación del 0 al 4

    file 2 (conteos_clasifTFs.tsv)
    (1) TFsTodocero
    (2) TFsConOutliers2o4
    (3) TFsConDominantes3o4

    Heatmap (Clasificación TFs, Gos)


USAGE
    R ClasifTF_contribMatrix_heatmap.R


SEE ALSO
    none
GitHub link
'''


# Cargar matriz de contribución

#matriz_contribución <- read.table()
#Cuatro_contribution_matrix_Spearman_Precise <- read.delim("C:/Users/equipo1/Desktop/Python/Genomica_Computacional/Genomica_Computacional_R/Cuatro_contribution_matrix_Spearman_Precise.txt", row.names=1)


### Cargar librerias
library(ComplexHeatmap)
library(circlize)
library(dplyr)

####################### Obtener los valores outliers ##########################################

# Contar las celdas distintas de cero en cada fila
conteo_distintos_cero <- rowSums(Cuatro_contribution_matrix_Spearman_Precise != 0)

# Ordenar las filas según el conteo de elementos distintos de cero
filas_ordenadas <- order(conteo_distintos_cero, decreasing = TRUE)

# Subconjunto del dataframe ordenado
df_ordenado <- Cuatro_contribution_matrix_Spearman_Precise[filas_ordenadas, ]

# Crear nuevo dataframe para guardar la información obtenida
tfOUTLIERS <- data.frame(GO = character(),
                         Cantidad_TFs_outliers = numeric(),
                         Nombre_TFs_outliers = character(),
                         Valor_outlier = character(),
                         Valor_outlier_ordenado = character(),
                         stringsAsFactors = FALSE)



for (i in (1:1152)){
  GO <- df_ordenado[i,]
  #
  idGO <- row.names(GO)
  # Eliminar columnas que tienen un valor de cero
  G0wo0 <- GO[, colSums(GO == 0) == 0]
  tGO <- t(G0wo0)
  vGO <- as.vector(tGO)

  # Calcular Q1 y Q3
  Q1 <- quantile(vGO, 0.25)
  Q3 <- quantile(vGO, 0.75)

  # Calcular IQR
  IQR <- Q3 - Q1

  # Calcular los límites
  limite_inferior <- Q1 - 1.5 * IQR
  limite_superior <- Q3 + 1.5 * IQR

  # Identificar los outliers
  outliers <- vGO[vGO < limite_inferior | vGO > limite_superior]
  outliers_ord <- sort(vGO[vGO < limite_inferior | vGO > limite_superior], decreasing = TRUE)

  # Convertir de nuevo la matriz transpuesta en un dataframe
  dfGO <-as.data.frame(tGO)

  if (length(outliers) != 0 ) {
    nombres_filas <- rownames(dfGO)[which(dfGO[,1] %in% outliers)]
    nombres_comas <- paste(nombres_filas, collapse = ",")
    val_outlier <- paste(outliers, collapse = ",")
    val_outlier_ord <- paste(outliers_ord, collapse = ",")
    # Llenar el dataframe con datos

    tfOUTLIERS <- rbind(tfOUTLIERS, data.frame(GO = idGO , Cantidad_TFs_outliers
                                               = length(outliers), Nombre_TFs_outliers = nombres_comas, Valor_outlier
                                               = val_outlier, Valor_outlier_ordenado = val_outlier_ord ))


  }

}

################### Clasificar los TFs dominantes ##############################

Cua_con_mat <- Cuatro_contribution_matrix_Spearman_Precise

# Eliminar la columna de b0000
Cua_con_mat <- Cua_con_mat[,-233]

# Eliminar filas donde todos los valores son cero
Cua_con_mat <- Cua_con_mat[rowSums(Cua_con_mat != 0) > 0, ]

# Encontrar los TFs con el valor máximo por fila
max_tfs <- apply(Cua_con_mat, 1, function(row) {
  max_value <- max(row)  # Encontrar el valor máximo en la fila
  tfs_with_max <- names(row[row == max_value])  # Obtener los nombres de los TFs con el valor máximo
  paste(tfs_with_max, collapse = ", ")  # Unir los nombres en una cadena separada por comas
})

# Crear un nuevo dataframe con los resultados
Tfs_dominantes <- data.frame(Max_TF = max_tfs, row.names = row.names(Cua_con_mat))


################# Crear el dataframe de clasificación. Gos x Tfs ###############

df_heatmap <- data.frame(matrix(ncol = ncol(Cuatro_contribution_matrix_Spearman_Precise), nrow = 0))
colnames(df_heatmap) <- colnames(Cuatro_contribution_matrix_Spearman_Precise)


################# Df de los TFs que sean cero  #################################

# Quitar el b0000 porque sino no existe un Go cuyos TFs sean todos cero
df_woGO_B0 <- subset(Cuatro_contribution_matrix_Spearman_Precise, select = -c(b0000))


for (i in (1:length(Cuatro_contribution_matrix_Spearman_Precise[,1]))){
  GO <- rownames(Cuatro_contribution_matrix_Spearman_Precise)[i]
  Gos_ceros <- df_woGO_B0[i,]
  todos_cero <- all(Gos_ceros == 0)
  if (todos_cero) {
    fila_original <- Cuatro_contribution_matrix_Spearman_Precise[i,]
    df_heatmap <- rbind(df_heatmap,fila_original)
  }
}

#colnames(df_heatmap) <- colnames(Cuatro_contribution_matrix_Spearman_Precise)

## Marcar si los valores del TF b0000 son dominantes, o dominantes y outliers


for (i in (1:nrow(df_heatmap))){
  GO <- df_heatmap[i,]
  # Calcular el primer y tercer cuartil (Q1 y Q3)
  Q1 <- quantile(GO, 0.25, na.rm = TRUE)
  Q3 <- quantile(GO, 0.75, na.rm = TRUE)

  # Calcular el rango intercuartílico (IQR)
  IQR <- Q3 - Q1

  # Definir los límites para identificar outliers
  limite_inferior <- Q1 - 1.5 * IQR
  limite_superior <- Q3 + 1.5 * IQR

  # Filtrar los TFs que son outliers
  outliers <- (GO < limite_inferior | GO > limite_superior)
  if(outliers[1,233]){
    df_heatmap[i,233] <- 4
  }else {
    if (df_heatmap[i,233] != 0){
      df_heatmap[i,233] <- 3 # En caso de que no tenga valor el TF B0000 que se quede como cero
    }

  }
}


##### Agregar los TFs dominantes con el valor de 3 como clasificatorio #########

## Esto se puede quitar porque ya cambie por un rownames
#Tfs_dominantes <- cbind(Gos = rownames(Tfs_dominantes), Tfs_dominantes)

# Normalizar nombres de columnas en df_heatmap
colnames(df_heatmap) <- trimws(tolower(colnames(df_heatmap)))


for (i in (1:nrow(Cuatro_contribution_matrix_Spearman_Precise))){
  GO <- rownames(Cuatro_contribution_matrix_Spearman_Precise)[i]
  for (y in (1:length(Tfs_dominantes[,1]))){
    GO_dom <- rownames(Tfs_dominantes)[y]
    if (GO == GO_dom){
      # Crear una nueva fila con el nombre en la primera columna y NA en las demás
      nueva_fila_GO <- c(rep(NA, ncol(df_heatmap)))
      # Agregar la nueva fila al dataframe
      df_heatmap <- rbind(df_heatmap, nueva_fila_GO)
      # Obtener el índice de la última fila agregada
      last_row_index <- nrow(df_heatmap)
      # agregar el nombre del GO al nombre del renglón
      rownames(df_heatmap)[last_row_index] <- GO_dom
      TFsdom <- Tfs_dominantes[y,1]
      cadena_dividida <- strsplit(TFsdom, split = ",")
      vector_tfsdom <- unlist(cadena_dividida)
      # Normalizar valores en vector_tfsdom
      vector_tfsdom <- trimws(tolower(vector_tfsdom))
      for (j in (1:length(vector_tfsdom))){
        df_heatmap[[vector_tfsdom[j]]][last_row_index] <- 3
      }
    }
  }
}


nombre_delGO <- "GO:0006355"
print(Cuatro_contribution_matrix_Spearman_Precise [nombre_delGO,])
print(df_heatmap[241,])


######### Agregar los TFs que son outliers con el identificador de 2 ###########

Gos_outliers <- tfOUTLIERS[,1]


for (i in (1:nrow(Cuatro_contribution_matrix_Spearman_Precise))){
  #GO <- df_heatmap[i,1] #i
  GO <- rownames(df_heatmap)[i]
  #index_GO <- which(df_heatmap$gos == GO)
  #index_GO <- which(rownames(df_heatmap)[i] == GO)
  for (y in (1:length(Gos_outliers))){
    GO_out <- Gos_outliers[y]
    if (GO == GO_out){
      Gos_uniq_outliers <- tfOUTLIERS[(tfOUTLIERS[,1] == GO),]
      TFsOutliers <- Gos_uniq_outliers[,3]
      cadena_dividida <- strsplit(TFsOutliers, split = ",")
      vector_ids <- unlist(cadena_dividida)
      vector_ids <- trimws(tolower(vector_ids))
      for (j in (1:length(vector_ids))){
        tf_out <- vector_ids[j]
        for (k in (1 : length(df_heatmap))){
          valor <- df_heatmap[i,k]
          if (isTRUE(df_heatmap[[tf_out]][i] == 3)){ ## agregué este para seleccionar los dominantes y outliers
            df_heatmap[[tf_out]][i] <- 4
          }
          if (is.na(df_heatmap[[tf_out]][i])){
            df_heatmap[[tf_out]][i] <- 2
          }
        }

      }
    }
  }
}


####### Agregar los valores de cero a las casillas que lo deban tener ##########

for (i in (1: nrow(Cuatro_contribution_matrix_Spearman_Precise))){
  GO <- Cuatro_contribution_matrix_Spearman_Precise[i,]
  GOid <- row.names(GO)
  for(j in(1:length(df_heatmap))){
    valor <- GO[1,j]
    if (valor == 0){
      df_heatmap[GOid,j] <- 0
    }
  }
}

###################### Convertir el resto de NAs a 1 ###########################

df_heatmap[is.na(df_heatmap)] <- 1


################## Guardar df_heatmap como archivo de texto ####################

# Obtener el nombre del dataframe principal como texto
nombre_df_principal <- deparse(substitute(Cuatro_contribution_matrix_Spearman_Precise))

# Crear nombre del archivo dinámico
nombre_archivo <- paste0("df_clasificacionTFs_", nombre_df_principal, ".txt")

# Guardar el dataframe
write.table(df_heatmap, file = nombre_archivo, sep = "\t", row.names = TRUE)


# Convertir todos los valores a numéricos (los no numéricos se harán NA)
valores <- as.numeric(unlist(df_heatmap))
conteos <- tabulate(valores + 1, nbins = 5)# +1 porque tabulate empieza en 1
names(conteos) <- 0:4
conteos

##################### Cambiar el bnumber del TF por su nombre ##################


# Convertir el nombre de las columnas del dataframe df_heatmap del identificador
# bnumb a su nombre original

TF_nomb <- RegulonDB_Network_bnumbersIDs[,1:2]

# Supongamos que tienes los siguientes dataframes:
# df: Tu dataframe principal con columnas identificadas por bnumb.
# df_info: Dataframe con dos columnas: bnumb y nombre.

# Paso 1: Eliminar duplicados en df_info
TF_nomb <- unique(TF_nomb)

# Dejar solo la primer opción, y eliminar los demás

# Paso 1: Quedarse solo con el primer nombre para cada bnumb
TF_nomb_unicos <- TF_nomb %>%
  group_by(X1.regulatorId) %>%          # Agrupar por bnumb
  slice(1) %>%                 # Tomar solo la primera fila de cada grupo
  ungroup()                    # Desagrupar

# Normalizar los formatos

# Paso 2: Crear el mapeo
mapeo <- setNames(TF_nomb_unicos$X2.regulatorName, TF_nomb_unicos$X1.regulatorId)

# Paso 3: Renombrar columnas, conservando las que no están en el mapeo
nombres_originales <- colnames(df_heatmap)
nombres_nuevos <- ifelse(
  nombres_originales %in% names(mapeo),  # Si el bnumb está en el mapeo
  mapeo[nombres_originales],             # Usar el nombre correspondiente
  nombres_originales                     # Si no, conservar el original (ej. "b0000")
)

colnames(df_heatmap) <- nombres_nuevos



############# Visualizar Heatmap del df clasificatorio #########################


max_df_heatmap <- as.matrix(df_heatmap)
# Cambiar el modo de la matriz a numeric
mode(max_df_heatmap) <- "numeric"

#write.table(max_df_heatmap, file = "matriz_clasificacionTFs.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


# Calcular distancias y clustering para columnas (eje X)
dist_columnas <- dist(t(max_df_heatmap))  # Transponer para columnas
hclust_columnas <- hclust(dist_columnas, method = "ward.D")

# Opcional: Ver el dendrograma
plot(hclust_columnas, main = "Dendrograma (ward.D)")

# Definir colores para cada número
colores <- c("0" = "gray", "1" = "green", "2" = "blue", "3" = "yellow", "4" = "red")

nom_GO <- c("0: TFs cero", "1: 0 < TF < Dominante/outlier","2: TF outlier", "3: TF dominante", "4: TF dominante y outlier")

# Crear la leyenda
leyenda <- Legend(at = 0:4,
                  labels = nom_GO,
                  title = "Clasificación TFs",
                  legend_gp = gpar(fill = colores))

hm <- Heatmap(
  max_df_heatmap,
  name = "Clasificación TFs",
  col = colores,
  show_heatmap_legend = FALSE,
  cluster_columns = hclust_columnas,  # Usar el clustering ward.D en columnas
  cluster_rows = TRUE,               # Clustering por defecto en filas (opcional)
  row_names_side = "right",
  column_names_side = "top",
  # Ajustar el tamaño de los nombres de las columnas (eje X)
  column_names_gp = gpar(fontsize = 8),  # Tamaño de fuente más pequeño (ej: 8pt)
  clustering_method_rows = "complete",  # Método para filas (opcional)
  clustering_distance_columns = "euclidean"  # Distancia usada en columnas
)

# Dibujar el heatmap con la leyenda
draw(hm, annotation_legend_list = leyenda)

help(heatmap)


########################## Cálculos ############################################

df_heatmap_WOB0 <- df_heatmap[,-233]
# Cálculo de TFs con solo ceros en todos los GOs
TFsTodocero <- sum(colSums(df_heatmap_WOB0 != 0) == 0)

# Cálculo de TFs con al menos un valor 2 o 4 (outliers)
TFsConOutliers2o4 <- sum(apply(df_heatmap_WOB0, 2, function(x) any(x %in% c(2, 4))))

# Cálculo de TFs con al menos un valor 3 o 4 (dominantes)
TFsConDominantes3o4 <- sum(apply(df_heatmap_WOB0, 2, function(x) any(x %in% c(3, 4))))

# Crear el dataframe resumen
resumen_TFs <- data.frame(
  TFsTodocero = TFsTodocero,
  TFsConOutliers2o4 = TFsConOutliers2o4,
  TFsConDominantes3o4 = TFsConDominantes3o4
)

# Mostrar el resultado
print(resumen_TFs)


TFs_todo_cero <- colnames(df_heatmap)[colSums(df_heatmap_WOB0 != 0) == 0]

# Mostrar los nombres
print(TFs_todo_cero)


write.table(resumen_TFs, file = "conteos_clasifTFs.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)







