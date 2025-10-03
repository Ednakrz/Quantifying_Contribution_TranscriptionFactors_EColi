'''
NAME
        4-contrib_matrix_generator.py
VERSION
        2.0.0
AUTHOR
        Victor Ulises Plascencia Perez / Edna Rivera 
DESCRIPTION
    Este programa genera una matriz de contribucion a partir de una matriz de coexpresion y una red de regulacion para los GOs involucrados en procesos biologicos.
    Genera un archivo con los genes que estan el archivo go, pero no en la matriz de coexpresion (lost_genes_coexpresion.txt)
    Genera uGenera un archivo con los genes que estan el archivo go y en la matriz de coexpresion, pero no en la red de regulacion (lost_genes_regulation.txt)
    
USAGE
    python 4-contrib_matrix_generator.py [-h] -c coexmat -r red -g gogene [--version]

ARGUMENTS
    -c -> ruta del archivo de coexpresion
    -r -> ruta con la red de regulacion
    -g -> archivo con los genes de cada GO

Input format 
    file 1 (coexpmat)
    columns: genes IDs (bnumber)
    indexes : genes IDs (bnumber)
    
    file 2 (red)
    (1) regulatorId. Transcription Factor (TF) identifier
    (2) regulatorName. Transcription Factor (TF) Name
    (3) regulatorGeneName. Gene(s) coding for the TF
    (4) regulatedId. Gene ID regulated by the TF (regulated Gene)
    (5) regulatedName. Gene regulated by the TF (regulated Gene)
    (6) function. Regulatory Function of the TF on the regulated Gene (+ activator, - repressor, +- dual, ? unknown)
    (7) confidenceLevel. RI confidence level based on its evidence (Values: Confirmed, Strong, Weak)
    
    file 2 (go gene)
    (1) GO
    (2) GeneName
    (3) GeneID
    
Output format
    file 1(contribution_matric.txt)
    columns: TFs IDs (bnumber)
    indexes : GOs
    
    file 2 (lost_genes_coexpresion.txt)
    (1) GO
    (2) LostGenes
    (3) NumGenes
    
    file 3 (lost_gene_regulation.txt)
    (1) GO
    (2) LostGenes
    (3) NumGenes

Execution example

python .\src\4-contrib_matrix_generator.py -c .\data\pearson_colombos.txt -r .\results\RegulonDB_Network_bnumbersIDs.txt -g .\results\GOID_bnumb_BP.txt

SEE ALSO
    none 
GitHub link     
'''

# Importar las librerias necesarias para los analisis
import pandas as pd
import numpy as np
import random as rd
import argparse  


# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(description="Genera un archivo con la red de regulacion de RegulonDB con bnumber de genes en lugar de ECKID")
parser.add_argument('-c','--coexpmat', metavar='coexmat', required= True, help='Ruta del archivo de coexpresion')
parser.add_argument('-r', '--red', metavar='red', required=True, help='Ruta del archivo con la red de regulacion con bnumbers')
parser.add_argument('-g', '--go', metavar='gogene', required=True, help='Ruta del archivo con los genes de cada GO')
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

# Leer los argumentos desde la terminal
args = parser.parse_args()
coexpmat = args.coexpmat  
red = args.red
go_gene = args.go

## Leer los archivos
# Gene Network
# Archivo de GOs_geneRegulates
GOs_file = pd.read_table(go_gene, sep='\t')  # para leer archivos de texto en un DataFrame
# Importar red de regulacion con bnumber
reg_network = pd.read_table(red, sep='\t')
# Importar matriz de coexpresion
coexp_matrix = pd.read_table(coexpmat, sep=',')


# Hacemos una lista con los GOs unicos y los TFs unicos
GOs = GOs_file.loc[:,'GO'].drop_duplicates().tolist()
regulators = reg_network.loc[:,'1)regulatorId'].drop_duplicates()


### Analisis mayores (core del programa)
# Dataframe para almacenar los datos
contrib_mat = pd.DataFrame(columns= regulators, index=GOs)
lost_genes_coex = pd.DataFrame(columns=['LostGenes', 'NumGenes'])
lost_genes_reg = pd.DataFrame(columns=['LostGenes', 'NumGenes'])


for go in GOs:
    hyp_regulator = pd.DataFrame() 
    go_section = GOs_file[GOs_file.loc[:,'GO'] == go] # Extraer la matriz de coexpresion de los genes del GO
    indexes_coexp_mat = coexp_matrix.index #Obtener los indices de la Matriz de Coexpresion
    go_coexp_mat = coexp_matrix.loc[indexes_coexp_mat.isin(go_section.loc[:,'GeneID']),indexes_coexp_mat.isin(go_section.loc[:,'GeneID'])].fillna(0)
    indexes_go_coexp_mat = go_coexp_mat.index #genes del go en la matriz de regulacion
    lost_genes_coexp_num = len(go_section) - len(go_coexp_mat.index) # Cantidad e identidad de genes del GO que no estan presentes en la coexpresion Matrix 
    if lost_genes_coexp_num > 0:
        lost_gene_go = go_section.loc[:,'GeneID'][~go_section.loc[:,'GeneID'].isin(go_coexp_mat.index)].tolist() #Identificar y Almacenar los Genes Perdidos el ~ es para selecciónar los que no cumplen con esa condición 
        lost_genes_coex.loc[go, 'LostGenes'] = lost_gene_go  #Renglon con el valor de go? aunque no lo hayan designado antes? o le asigna el nombre de go a la primer fila? y así sucesivamente ?
        lost_genes_coex.loc[go, 'NumGenes'] = lost_genes_coexp_num
        
    
    # Verificar si la matrix de coexp esta vacia
    
    if go_coexp_mat.empty:
        continue

    go_reg_network = reg_network[reg_network.loc[:,'4)regulatedId'].isin(go_coexp_mat.index)] 
    regulator_genes = go_reg_network.loc[:,'1)regulatorId']
    regulated_genes = go_reg_network.loc[:,'4)regulatedId']
   
   
    # Determinar la cantidad e identidad de genes de la matriz de coexpresion ausentes en la red
    lost_genes_reg_num =  len(go_coexp_mat) - len(regulated_genes.unique())
    
    if lost_genes_reg_num > 0: 
        lost_genes_reg_go = indexes_go_coexp_mat[~indexes_go_coexp_mat.isin(regulated_genes.unique())]
        hyp_regulator = pd.DataFrame({'4)regulatedId': lost_genes_reg_go})  
        hyp_regulator.loc[:, '1)regulatorId'] = 'b0000'  
        lost_genes_reg.loc[go, 'LostGenes'] = lost_genes_reg_go  # Guardar los genes que no esan en la red de regulacion
        lost_genes_reg.loc[go, 'NumGenes'] = lost_genes_reg_num
        
    
    # Generacion de matriz de regulacion
    TF_gene_DF = pd.concat([regulator_genes, regulated_genes], axis=1) #Une los dataframes por columnas 
   
   
    if lost_genes_reg_num > 0:
        hyp_row = hyp_regulator[['1)regulatorId','4)regulatedId']]
        TF_gene_DF = pd.concat([TF_gene_DF,hyp_row], ignore_index=True)
        
    hyp_regulator = 0 
    TF_gene_DF['relation'] = 1 # Crea una nueva columna llamada 'relation' en el DataFrame TF_gene_DF y le asigna el valor 1 para indicar la relación entre reguladores y genes regulados.
    reg_matrix = pd.pivot_table(TF_gene_DF, index='4)regulatedId',columns='1)regulatorId', values='relation',fill_value=0) #Una matriz que representa la relación entre reguladores y genes regulados, donde los valores son 1 si hay una relación y 0 si no la hay.
    
    # Pasar de una matrix triangular a una simetrica
    go_coexp_mat = go_coexp_mat.T + go_coexp_mat
    
    
    # Operaciones matemáticas
    tfs = reg_matrix.columns
    go_coexp_mat = go_coexp_mat.to_numpy()
    np.fill_diagonal(go_coexp_mat, 0)
    reg_matrix = reg_matrix.to_numpy()
    go_coexp_mat_float = go_coexp_mat.astype(float)
    go_int_mat = go_coexp_mat@reg_matrix  
    # Redondear
    vectoGo = np.sum(go_int_mat, axis=0)
    
    vectoGo = np.around(vectoGo, decimals=2) #Redondea los valores a 2 decimales
    go_tf_vector = pd.Series(vectoGo, index=tfs) 
    contrib_mat.loc[go,tfs] = go_tf_vector 
    vectoGo = 0 
    
# Guardar los archivos generados
contrib_mat.round(1).fillna(0).to_csv('results/Cuatro_contribution_matrix_Spearman_Precise.txt', sep='\t') 
lost_genes_coex.to_csv('results/lost_genes_coexpresion_Spearman_Precise.txt', sep='\t')
lost_genes_reg.to_csv('results/lost_genes_regulation_Spearman_Precise.txt', sep='\t')

