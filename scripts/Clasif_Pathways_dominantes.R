'''
NAME
        Clasif_Pathways_dominantes.R
VERSION
        1.0.0
AUTHOR
       Edna Karen Rivera Zagal

DESCRIPTION
    Código estructurado para Obtener la clasificación de los Pathways en la matriz de
    contribución.

Input format
    file 1(contribution_matrix.txt)
    columns: TFs IDs (bnumber)
    indexes : Pathways


Output format
    file 1 (df_Pathways_dominantes_SpearmanPrecise_names.tsv)
    (1) TF
    (2) Valor_dominante
    (3) Cantidad_Pathways
    (4) Pathways_dominantes
    (5) Pathway_id

USAGE
    R Clasif_Pathways_dominantes.R


SEE ALSO
    none
GitHub link
'''

# Clasificación Pathways dominantes por TF

# Librerías
library(dplyr)
library(tidyr)
library(stringr)

# Cambiar el nombre de los Pathways

# Extraemos los nombres de fila (los IDs)

contribution_matrix_Spearman_Precise_Pathway <- read.table(
  "contribution_matrix_Spearman_Precise_Pathway.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1
)

ids <- rownames(contribution_matrix_Spearman_Precise_Pathway)

# Creamos un vector de nombres GO con base en los IDs
nombres_Path <- genes_of_pathways_c$PathwayName[match(ids, genes_of_pathways_c$Pathway)]

# Quitamos las filas donde no se encontró nombre (opcional)
# También podrías quedarte con todas y dejar NAs
filas_validas <- !is.na(nombres_Path)
contribution_matrix_Spearman_Precise_Pathway_name <- contribution_matrix_Spearman_Precise_Pathway[filas_validas, ]

# Asignamos los nombres GO como nuevos rownames
rownames(contribution_matrix_Spearman_Precise_Pathway_name) <- nombres_Path[filas_validas]

## Cambiar al nombre del TF

# Suponiendo que TF_nomb_unicos tiene columnas: V1 = idTF, V2 = TFname
colnames(contribution_matrix_Spearman_Precise_Pathway_name) <-
  TF_nomb_unicos$X2.regulatorName[match(colnames(contribution_matrix_Spearman_Precise_Pathway_name), TF_nomb_unicos$X1.regulatorId)]

## ****** Encontrar los dominantes

# En este caso en particular cambiar la columna 233 a b000 porque no tenía nada
colnames(contribution_matrix_Spearman_Precise_Pathway_name)[233] <- "b000"
colnames(contribution_matrix_Spearman_Precise_Pathway_name) <- make.names(colnames(contribution_matrix_Spearman_Precise_Pathway_name), unique = TRUE)


# Paso 1: Convertir rownames (GO IDs) a columna explícita
df_prep <- contribution_matrix_Spearman_Precise_Pathway_name %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Pathway")  # Convierte row names a columna "GOid"

# Paso 2: Pivotear a formato largo (PathwayID, TF, valor)
df_long <- df_prep %>%
  pivot_longer(
    cols = -Pathway,
    names_to = "TF",
    values_to = "valor"
  )

# Paso 3: Encontrar GOid(s) dominante(s) por TF
df_resultado <- df_long %>%
  group_by(TF) %>%
  filter(valor == max(valor, na.rm = TRUE)) %>%
  summarise(
    Pathways_dominantes = paste(Pathway, collapse = "###"),  # separador seguro
    Valor_dominante = first(valor),
    Cantidad_Pathways = n()
  ) %>%
  arrange(desc(Valor_dominante), TF)

## Agregar la columna de los ids de los pathways de nuevo xd

df_resultado <- df_resultado %>%
  separate_rows(Pathways_dominantes, sep = "###") %>%
  left_join(
    genes_of_pathways_c %>% distinct(Pathway, PathwayName),
    by = c("Pathways_dominantes" = "PathwayName")
  ) %>%
  group_by(across(all_of(names(df_resultado)))) %>%
  summarise(
    Pathway_id = paste(unique(Pathway), collapse = ", "),
    .groups = "drop"
  )

## Colapsar todo por TF

df_resultado <- df_resultado %>%
  group_by(TF, Valor_dominante, Cantidad_Pathways) %>%
  summarise(
    Pathways_dominantes = paste(unique(Pathways_dominantes), collapse = ", "),
    Pathway_id = paste(unique(Pathway_id), collapse = ", "),
    .groups = "drop"
  )

write.table(df_resultado, file = "df_Gos_dominantes_SpearmanPrecise_Pathways_names.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)



