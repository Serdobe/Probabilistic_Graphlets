## MAIN Function

# General Packages:

import pandas as pd
from multiprocessing import Pool

# Personal Packages

import multiprocessing
import Distance_Script
import Enrichment_Script_All
import Plots_script_Enrichment
import Script_Parallel_Nodes_Count

# Description:

'''
This script contains all the steps to develop the whole protocol to compare
binary against probabilistics network. To use it is mandatory to load each of the
libraries in python and some in C++.
'''

# Main:

Name_Network_Original = "Coex_Net_Original_Low_Rank_10%.txt"

file_Directory = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/"
save_directory = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/" 
Run_path = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Node_Parallel/Release/"


# 1st STEP (Creation and computation of the orbits from a file) #

# First Create the Networks and Store the Names:

Script_Parallel_Nodes_Count.Creat_Network(Name_Network_Original, save_directory)

# Load the names to use them in the enrichment:

Names_Genes = pd.read_csv(save_directory +"_" + "Gene_Names_" + Name_Network_Original,
                              sep = " ", header = None)  # Not sure if separation should be a comma (10/01/20))

# A) Then we compute in parallel the Orbits (For Probabilistic and Binary):

Name_Network_Prob = ['Work_Network' + "_" + "Prob_" + Name_Network_Original]
Name_Network_Bin = ['Work_Network' + "_" + "Bin_" + Name_Network_Original]

total_nodes = (len(open(save_directory + "_" + "Gene_Names_" + Name_Network_Original).readlines(  ))) - 1 #Because the first line

Script_Parallel_Nodes_Count.Compute_Orbits_Spliting_Nodes(Name_Network_Bin[0], total_nodes,
                                                          file_Directory, save_directory,
                                                          Run_path, threads = 7)

Script_Parallel_Nodes_Count.Compute_Orbits_Spliting_Nodes(Name_Network_Prob[0], total_nodes,
                                                          file_Directory, save_directory,
                                                          Run_path, threads = 7)

# 2nd STEP (Calculation of the distance):

# Prepare matrix to the calculations:

GDV_Prob = Distance_Script.Prepare_Matrix_GDV(Name_Network_Prob[0], save_directory)
GDV_Bin = Distance_Script.Prepare_Matrix_GDV(Name_Network_Bin[0], save_directory)

Distance_Script.Prepare_To_Distance(GDV_Prob, Name_Network_Prob[0], save_directory)
Distance_Script.Prepare_To_Distance(GDV_Bin, Name_Network_Bin[0], save_directory)

# Calculate Tijana distance (in this case):

Distance_Script.Run_Tijana_Distance(Name_Network_Prob[0], save_directory, Name_Network_Original)
Distance_Script.Run_Tijana_Distance(Name_Network_Bin[0], save_directory, Name_Network_Original)

# From here, the user will have enough information to perform an enrichment analyses and 
# all the rest of the analyses exposed in the paper.

# 3rd STEP (Enrichment):

# Annotation MF, BP and CC are, in this moment, prepared for Yeastract gene names

# Paths to the annotations:

path_BP = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_BP_Levels.csv" 
path_CC = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_CC_Levels.csv"
path_MF = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_MF_Levels.csv"


# We do it in parallel, to check the arguments go to the original function:

sub_command1 = (Name_Network_Prob[0], 
                save_directory, 
                path_BP, Name_Network_Original,"BP", 
                                                   5000000, 
                                                   1,
                                                   10, 
                                                   100, 
                                                   66,
                                                   "Bin_VS_Prob")

sub_command2 = (Name_Network_Prob[0], save_directory, path_CC, Name_Network_Original,"CC", 
                                                   5000000, 
                                                   1,
                                                   10, 
                                                   100, 
                                                   66,
                                                   "Bin_VS_Prob")

sub_command3 = (Name_Network_Prob[0], save_directory, path_MF, Name_Network_Original,"MF", 
                                                   5000000, 
                                                   1,
                                                   10, 
                                                   100, 
                                                   66,
                                                   "Bin_VS_Prob")

sub_command4 = (Name_Network_Bin[0], save_directory, path_BP, Name_Network_Original,"BP", 
                                                   5000000, 
                                                   1,
                                                   10, 
                                                   100, 
                                                   66,
                                                   "Bin_VS_Prob")

sub_command5 = (Name_Network_Bin[0], save_directory, path_CC, Name_Network_Original,"CC", 
                                                   5000000, 
                                                   1,
                                                   10, 
                                                   100, 
                                                   66,
                                                   "Bin_VS_Prob")

sub_command6 = (Name_Network_Bin[0], save_directory, path_MF, Name_Network_Original,"MF", 
                                                   5000000, 
                                                   1,
                                                   10, 
                                                   100, 
                                                   66,
                                                  "Bin_VS_Prob")
  

Process = [sub_command1,
                   sub_command2,
                   sub_command3,
                   sub_command4,
                   sub_command5,
                   sub_command6]
         
n_cores = multiprocessing.cpu_count()    
pool1 = Pool(processes  = n_cores)  
pool1.starmap(Enrichment_Script_All.Enrichment_Analyses_GO_terms, Process)
            
pool1.close()
pool1.join() 

