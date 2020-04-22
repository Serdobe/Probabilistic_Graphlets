# -*- coding: utf-8 -*-

"""
# ===========================================
#                             _           
#     /\                     | |          
#    /  \   _ __  _ __  _   _| |_   _ ___ 
#   / /\ \ | '_ \| '_ \| | | | | | | / __|
#  / ____ \| | | | | | | |_| | | |_| \__ \
# /_/    \_\_| |_|_| |_|\__,_|_|\__,_|___/
#                                         
# ===========================================
"""   


from goatools import obo_parser

import pandas as pd
import Bio.UniProt.GOA as GOA

def take_my_Y(synonym_list):
    for gene in synonym_list:
        if gene.startswith('Y'):
            return gene

        
GO_file = '/home/sergio/workspace_Eclipse/Lucky_GOGO/Results/Python_Srcipts/go-basic.obo'
sc_GAF_file = '/home/sergio/workspace_Eclipse/Lucky_GOGO/Results/Python_Srcipts/sgd.gaf'
go_dag = obo_parser.GODag(GO_file)


with open(sc_GAF_file, 'rt') as fp:
    sc_gaf = pd.DataFrame(annotation for annotation in GOA.gafiterator(fp))

    
sc_gaf = sc_gaf[sc_gaf['Evidence'].isin(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'])]
sc_gaf = sc_gaf[sc_gaf['Aspect']=='F']
sc_gaf['Yeast_ID'] = [take_my_Y(gene) for gene in sc_gaf['Synonym']]


sergio_df = sc_gaf[['Yeast_ID', 'GO_ID']].drop_duplicates()
new_sergio_df = pd.DataFrame(columns=['Yeast_ID','GO_ID'])

for _, row in sergio_df.iterrows():
    new_sergio_df = new_sergio_df.append(row)
    parents = pd.DataFrame([(row.Yeast_ID,parent_id) for parent_id in go_dag[row.GO_ID].get_all_parents()
                           if go_dag[parent_id].namespace=='molecular_function'])
    parents.columns = ['Yeast_ID','GO_ID']
    new_sergio_df = new_sergio_df.append(parents)

new_sergio_df['Level'] = [go_dag[term].level for term in new_sergio_df['GO_ID']]
new_sergio_df = new_sergio_df[new_sergio_df.Level != 0]

save_path_Table= "/home/sergio/workspace_Eclipse/Lucky_GOGO/Results/Python_Srcipts/Final_Results/Annotation Script/Molecular Function/Annotation_Molecular_Function.csv"
            
new_sergio_df.to_csv(save_path_Table, sep = "\t", header = True, index = False)
















