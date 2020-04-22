from goatools import obo_parser
import pandas as pd
import Bio.UniProt.GOA as GOA

def take_my_Y(synonym_list):
    for gene in synonym_list:
        if gene.startswith('Y'):
            return gene

        
GO_file = '/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Python_Srcipts/go-basic.obo'

# In this case we start from a txt file with all the links between GO and Genes:

BP_path = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_BP.csv"
CC_path = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_CC.csv"
MF_path = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_MF.csv"

BP_Data = pd.read_csv(BP_path, sep = "\t", header = 0)
CC_Data = pd.read_csv(CC_path, sep = "\t", header = 0)
MF_Data = pd.read_csv(MF_path, sep = "\t", header = 0)

go_dag = obo_parser.GODag(GO_file)

# Using the parser we can determine the levels of our GO terms:

BP_Data['Level'] = [go_dag[term].level for term in BP_Data['GO_ID']]
CC_Data['Level'] = [go_dag[term].level for term in CC_Data['GO_ID']]
MF_Data['Level'] = [go_dag[term].level for term in MF_Data['GO_ID']]

BP_Data = BP_Data[BP_Data.Level != 0]
CC_Data = CC_Data[CC_Data.Level != 0]
MF_Data = MF_Data[MF_Data.Level != 0]

# Save:

save_BP = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_BP_Levels.csv"
save_CC = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_CC_Levels.csv"
save_MF = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_MF_Levels.csv"
       
BP_Data.to_csv(save_BP, sep = "\t", header = True, index = False)
CC_Data.to_csv(save_CC, sep = "\t", header = True, index = False)
MF_Data.to_csv(save_MF, sep = "\t", header = True, index = False)


