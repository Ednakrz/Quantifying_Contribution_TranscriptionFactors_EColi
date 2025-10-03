'''
NAME
        Clasif_allTFsPathways.R
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
    file 1 (Clasif_TFxPathways.tsv)
    (1) TF
    (2) Pathway
    (3) Valor
    (4) Pathway_id
    (5) Cantidad_Pathways

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
# Guardamos el orden original de los TFs
orden_TFs <- colnames(contribution_matrix_Spearman_Precise_Pathway_name)

df_pathways <- contribution_matrix_Spearman_Precise_Pathway_name %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Pathway") %>%   # Pathway names como columna
  pivot_longer(
    cols = -Pathway,
    names_to = "TF",
    values_to = "valor"
  ) %>%
  filter(valor != 0) %>%  # solo valores distintos de cero
  # renombramos la columna de ID antes de unir
  left_join(
    genes_of_pathways_c %>% distinct(Pathway, PathwayName) %>%
      rename(Pathway_id = Pathway),
    by = c("Pathway" = "PathwayName")
  ) %>%
  group_by(TF) %>%
  mutate(Cantidad_Pathways = n_distinct(Pathway)) %>%
  ungroup() %>%
  # asegurar orden de TFs según la matriz original
  mutate(TF = factor(TF, levels = colnames(contribution_matrix_Spearman_Precise_Pathway_name))) %>%
  arrange(TF) %>%
  # reordenar columnas
  select(TF, Pathway, valor, Pathway_id, Cantidad_Pathways)

write.table(df_pathways, file = "Clasif_TFxPathways.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

