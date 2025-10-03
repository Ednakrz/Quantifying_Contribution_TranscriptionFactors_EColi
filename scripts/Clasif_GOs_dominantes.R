'''
NAME
        Clasif_GOs_dominantes.R
VERSION
        1.0.0
AUTHOR
       Edna Karen Rivera Zagal

DESCRIPTION
    Código estructurado para Obtener la clasificación de los GOs en la matriz de
    contribución.

Input format
    file 1(contribution_matrix.txt)
    columns: TFs IDs (bnumber)
    indexes : GOs


Output format
    file 1 (df_Gos_dominantes_SpearmanPrecise_names.tsv)
    (1) TF
    (2) GOids_dominantes
    (3) Valor_dominante
    (4) Cantidad_GOids
    (5) Goids
    (6) GO_Levels

USAGE
    R Clasif_GOs_y_Pathways_dominantes.R


SEE ALSO
    none
GitHub link
'''

# Clasificación Gos dominantes por TF

# Crear un dataframe nuevo con los nombres por GO y por TF
# Matriz con los identificadores por TF y GO: Cuatro_contribution_matrix_Spearman_Precise.txt

# Librerías
library(dplyr)
library(tidyr)
library(stringr)


## ********* TF identificador a TF name ***************

# Cambiar el nombre de los Gos

# Extraemos los nombres de fila (los IDs)
ids <- rownames(Cuatro_contribution_matrix_Spearman_Precise)

# Creamos un vector de nombres GO con base en los IDs
nombres_GO <- ids_y_nombres_GO$V2[match(ids, ids_y_nombres_GO$V1)]

# Quitamos las filas donde no se encontró nombre (opcional)
# También podrías quedarte con todas y dejar NAs
filas_validas <- !is.na(nombres_GO)
Cuatro_contribution_matrix_Spearman_Precise_Goname <- Cuatro_contribution_matrix_Spearman_Precise[filas_validas, ]

# Asignamos los nombres GO como nuevos rownames
rownames(Cuatro_contribution_matrix_Spearman_Precise_Goname) <- nombres_GO[filas_validas]



# Código hasta abajo cod_clasif_TFdominante_GOdominante, untitle 1

# Suponiendo que TF_nomb_unicos tiene columnas: V1 = idTF, V2 = TFname
colnames(Cuatro_contribution_matrix_Spearman_Precise_Goname) <-
  TF_nomb_unicos$X2.regulatorName[match(colnames(Cuatro_contribution_matrix_Spearman_Precise_Goname), TF_nomb_unicos$X1.regulatorId)]

## ****** Encontrar los dominantes

# En este caso en particular cambiar la columna 233 a b000 porque no tenía nada
colnames(Cuatro_contribution_matrix_Spearman_Precise_Goname)[233] <- "b000"
colnames(Cuatro_contribution_matrix_Spearman_Precise_Goname) <- make.names(colnames(Cuatro_contribution_matrix_Spearman_Precise_Goname), unique = TRUE)


# Paso 1: Convertir rownames (GO IDs) a columna explícita
df_prep <- Cuatro_contribution_matrix_Spearman_Precise_Goname %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GO")  # Convierte row names a columna "GOid"

# Paso 2: Pivotear a formato largo (GOid, TF, valor)
df_long <- df_prep %>%
  pivot_longer(
    cols = -GO,
    names_to = "TF",
    values_to = "valor"
  )

# Paso 3: Encontrar GOid(s) dominante(s) por TF
df_resultado <- df_long %>%
  group_by(TF) %>%
  filter(valor == max(valor, na.rm = TRUE)) %>%
  summarise(
    GOids_dominantes = paste(GO, collapse = ", "),
    Valor_dominante = first(valor),
    Cantidad_GOids = n()  # Número de GOids dominantes (útil para identificar empates)
  ) %>%
  arrange(desc(Valor_dominante), TF)

write.table(df_resultado, file = "df_Gos_dominantes_SpearmanPrecise_names.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


## Cambiar el nombre de la fila de NA por B0000

df_Gos_dominantes_SpearmanPrecise_names_levels[169,1] <- "b000"
write.table(df_Gos_dominantes_SpearmanPrecise_names_levels, file = "df_Gos_dominantes_SpearmanPrecise_names_levels.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
