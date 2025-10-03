#!/bin/bash

# Uso: bash Generator_analysis_ContriMatrix.sh GOs  o  bash Generator_analysis_ContriMatrix.sh Pathways
DATATYPE=$1

if [ -z "$DATATYPE" ]; then
    echo "Uso: bash $0 <GOs|Pathways>"
    exit 1
fi

mkdir -p intermediates results

echo ">>> Corriendo pipeline para $DATATYPE"

# Paso 2-4: R (diferentes scripts según el tipo de dato)
if [ "$DATATYPE" == "GOs" ]; then
    python scripts/4-contrib_matrix_generator.py -c data/ precise_Spearman.txt -r data/RegulonDB_Network_bnumbersIDs.txt -g data/filter_GOID_bnumb.txt
    R scripts/Frecuencias_TF_porGO
    R scripts/Clasif_GOs_dominantes
    R scripts/Cod_clasifTF_contribMatrix_heatmap
    R scripts/Clasif_allTFSGOs

elif [ "$DATATYPE" == "Pathways" ]; then
    python 4-contrib_matrix_generator_pathway.py -c PRECISE_SPEARMAN/ precise_Spearman.txt -r RegulonDB_Network_bnumbersIDs.txt -p PRECISE_SPEARMAN/genes_of_pathways.txt
    R scripts/Clasif_Pathways_dominantes
    R scripts/Cod_clasifTF_contribMatrix_heatmap
    R scripts/Clasif_allTFSPathways
else
    echo "Error: argumento debe ser GOs o Pathways"
    exit 1
fi

echo ">>> Pipeline completado para $DATATYPE ✅"
