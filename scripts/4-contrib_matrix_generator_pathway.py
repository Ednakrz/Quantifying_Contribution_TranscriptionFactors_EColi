'''
NAME
        4-contrib_matrix_generator.py
VERSION
        1.0.0
AUTHOR
        Victor Ulises Plascencia Perez / Edna Rivera 
DESCRIPTION
    Este programa genera una matriz de contribucion a partir de una matriz de coexpresion y una red de regulacion para los Pathways involucrados en procesos biologicos.
    Genera un archivo con los genes que estan el archivo de Pathways, pero no en la matriz de coexpresion (lost_genes_coexpresion.txt)
    Genera uGenera un archivo con los genes que estan el archivo Pathways y en la matriz de coexpresion, pero no en la red de regulacion (lost_genes_regulation.txt)
    
USAGE
    python 4-contrib_matrix_generator.py [-h] -c coexmat -r red -p pathway [--version]

ARGUMENTS
    -c -> ruta del archivo de coexpresion
    -r -> ruta con la red de regulacion
    -p -> archivo con los genes de cada Pathway

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
    
    file 2 (Pathway)
    (1) Pathway
    (2) GeneName
    (3) GeneID
    (4) PathwayName
    
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

python .\src\4-contrib_matrix_generator.py -c .\data\pearson_colombos.txt -r .\results\RegulonDB_Network_bnumbersIDs.txt -f .\results\genes_of_patwhays_c.txt

SEE ALSO
    none 
GitHub link     
'''

# Importar las librerias necesarias para los analisi
import pandas as pd 
import numpy as np
import random as rd
import argparse   
 

# Definir argumentos opcionales y posicionales
parser = argparse.ArgumentParser(description="Genera un archivo con la red de regulacion de RegulonDB con bnumber de genes en lugar de ECKID")
parser.add_argument('-c','--coexpmat', metavar='coexmat', required= True, help='Ruta del archivo de coexpresion')
parser.add_argument('-r', '--red', metavar='red', required=True, help='Ruta del archivo con la red de regulacion con bnumbers')
parser.add_argument('-p', '--pat', metavar='patgene', required=True, help='Ruta del archivo con los genes de cada Pathway')
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

# Leer los argumentos desde la terminal
args = parser.parse_args()
coexpmat = args.coexpmat  
red = args.red
pat_gene = args.pat

## Leer los archivos
# Gene Network
Pats_file = pd.read_table(pat_gene, sep='\t')  
# Importar red de regulacion con bnumber
reg_network = pd.read_table(red, sep='\t')
# Importar matriz de coexpresion
coexp_matrix = pd.read_table(coexpmat, sep=',')


# Hacemos una lista con los Pathways unicos y los TFs unicos
Pats = Pats_file.loc[:,'Pathway'].drop_duplicates().tolist()
regulators = reg_network.loc[:,'1)regulatorId'].drop_duplicates()


### Analisis mayores (core del programa)
# Dataframe para almacenar los datos
contrib_mat = pd.DataFrame(columns= regulators, index=Pats)
lost_genes_coex = pd.DataFrame(columns=['LostGenes', 'NumGenes'])
lost_genes_reg = pd.DataFrame(columns=['LostGenes', 'NumGenes'])
for pat in Pats:
    hyp_regulator = pd.DataFrame() 
    # Extraer la matriz de coexpresion de los genes del GO
    pat_section = Pats_file[Pats_file.loc[:,'Pathway'] == pat] 
    
    indexes_coexp_mat = coexp_matrix.index 
    pat_coexp_mat = coexp_matrix.loc[indexes_coexp_mat.isin(pat_section.loc[:,'GeneID']),indexes_coexp_mat.isin(pat_section.loc[:,'GeneID'])].fillna(0)
    indexes_pat_coexp_mat = pat_coexp_mat.index 
    
    # Determinar la cantidad e identidad de genes del GO que no estan presentes en la matriz de coexpresion
    lost_genes_coexp_num = len(pat_section) - len(pat_coexp_mat.index)
    if lost_genes_coexp_num > 0:
        lost_gene_pat = pat_section.loc[:,'GeneID'][~pat_section.loc[:,'GeneID'].isin(pat_coexp_mat.index)].tolist()  
        lost_genes_coex.loc[pat, 'LostGenes'] = lost_gene_pat  
        lost_genes_coex.loc[pat, 'NumGenes'] = lost_genes_coexp_num
        
    
    # Verificar si la matrix de coexp esta vacia
    if pat_coexp_mat.empty:
        continue
    
    # Extraer los genes regulados y los reguladores
    pat_reg_network = reg_network[reg_network.loc[:,'4)regulatedId'].isin(pat_coexp_mat.index)] 
    regulator_genes = pat_reg_network.loc[:,'1)regulatorId']
    regulated_genes = pat_reg_network.loc[:,'4)regulatedId']


    # Determinar la cantidad e identidad de genes que de la matriz de coexpresion ausentes en la red
    lost_genes_reg_num =  len(pat_coexp_mat) - len(regulated_genes.unique())
    if lost_genes_reg_num > 0: 
        lost_genes_reg_pat = indexes_pat_coexp_mat[~indexes_pat_coexp_mat.isin(regulated_genes.unique())]
        hyp_regulator = pd.DataFrame({'4)regulatedId': lost_genes_reg_pat})  
        hyp_regulator.loc[:, '1)regulatorId'] = 'b0000' 
        # Guardar los genes que no esan en la red de regulacion
        lost_genes_reg.loc[pat, 'LostGenes'] = lost_genes_reg_pat
        lost_genes_reg.loc[pat, 'NumGenes'] = lost_genes_reg_num
    
    # Generacion de matriz de regulacion
    TF_gene_DF = pd.concat([regulator_genes, regulated_genes], axis=1) 
    
    if lost_genes_reg_num > 0:
        hyp_row = hyp_regulator[['1)regulatorId','4)regulatedId']]
        TF_gene_DF = pd.concat([TF_gene_DF,hyp_row], ignore_index=True)
        
    hyp_regulator = 0 
    TF_gene_DF['relation'] = 1 
    reg_matrix = pd.pivot_table(TF_gene_DF, index='4)regulatedId',columns='1)regulatorId', values='relation',fill_value=0) 
    

    # Operaciones matemáticas
    tfs = reg_matrix.columns
    pat_coexp_mat = pat_coexp_mat.to_numpy()
    np.fill_diagonal(pat_coexp_mat, 0)
    reg_matrix = reg_matrix.to_numpy()
    pat_coexp_mat_float = pat_coexp_mat.astype(float)
    pat_int_mat = pat_coexp_mat@reg_matrix 
    vectoPat = np.sum(pat_int_mat, axis=0)
    vectoPat = np.around(vectoPat, decimals=2) 
    pat_tf_vector = pd.Series(vectoPat, index=tfs) 
    contrib_mat.loc[pat,tfs] = pat_tf_vector 
    vectoPat = 0 
        
# Guardar los archivos generados
contrib_mat.round(1).fillna(0).to_csv('results/contribution_matrix_Spearman_Precise_Pathway.txt', sep='\t') 
lost_genes_coex.to_csv('results/lost_genes_coexpresion_Spearman_Precise_Pathway.txt', sep='\t')
lost_genes_reg.to_csv('results/lost_genes_regulation_Spearman_Precise_Pathway.txt', sep='\t')


