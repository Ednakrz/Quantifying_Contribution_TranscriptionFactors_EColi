'''
NAME
        Clasif_allTFsGOs.R
VERSION
        1.0.0
AUTHOR
       Edna Karen Rivera Zagal

DESCRIPTION
    Código estructurado para Obtener la clasificación de todos los Factores de
    Transcripción en la matriz de contribución.

Input format
    file 1(contribution_matrix.txt)
    columns: TFs IDs (bnumber)
    indexes : GOs


Output format
    file 1 (Clasif_TFxGOs_Goleves.tsv)
    (1) TF
    (2) GO
    (3) Valor
    (4) Goids
    (5) GO_Levels
    (6) Cantidad_GOids

USAGE
    R Clasif_allTFs.R


SEE ALSO
    none
GitHub link
'''

# Cargar las librerias para los analisis
library(dplyr)
library(tidyr)


# Obtencion de los datos
df_GO <- Cuatro_contribution_matrix_Spearman_Precise_Goname %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GO") %>%   # GO names como columna
  pivot_longer(
    cols = -GO,
    names_to = "TF",
    values_to = "valor"
  ) %>%
  filter(valor != 0) %>%  # nos quedamos solo con los GOs distintos de 0
  left_join(
    ids_y_nombres_GO %>% distinct(V1, V2),
    by = c("GO" = "V2")   # V2 es el nombre, match con GO
  ) %>%
  rename(GO_id = V1) %>%
  left_join(
    all_go_levels %>% distinct(GO.ID, GO.level),
    by = c("GO_id" = "GO.ID")
  ) %>%
  group_by(TF) %>%
  mutate(Cantidad_GOs = n_distinct(GO)) %>%
  ungroup()

df_GO <- df_GO %>%
  select(TF, GO, valor, GO_id, GO.level, Cantidad_GOs)

# Mantener los TFs contiguos

orden_TFs <- colnames(Cuatro_contribution_matrix_Spearman_Precise_Goname)
df_GO <- df_GO %>%
  mutate(TF = factor(TF, levels = orden_TFs)) %>%
  arrange(TF) %>%   # así todas las filas de un mismo TF quedan juntas
  select(TF, GO, valor, GO_id, GO.level, Cantidad_GOs)  # reordenar columnas

write.table(df_GO, file = "Clasif_TFxGOs_Goleves.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

