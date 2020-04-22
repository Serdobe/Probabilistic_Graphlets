
############################################
# Script to make the plots for the article #
############################################

# Description

'''
This script contains all the fucntions to generate the plots of the enrichment analyses
including Jaccard index and Go term mean level
'''

import random
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import comb
from collections import Counter
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from pyclustering.cluster.kmedoids import kmedoids
from statsmodels.stats.multitest import multipletests  
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

from ast import literal_eval 



# Plots:


def To_Jaccard_Index(Name_Prob, 
                   Name_Bin,
                   Repetitions, 
                   Rand_K_Compare,
                   save_directory,
                   enrichment):
    
    copy1 = pd.read_csv(save_directory + "_Rand_Information_"+ "_" + enrichment+ "_" + Name_Prob + ".txt",
                              sep = ",", header = None, skiprows=[0])
    
    copy2 = pd.read_csv(save_directory + "_Rand_Information_"+ "_" + enrichment+ "_" + Name_Bin + ".txt",
                              sep = ",", header = None, skiprows=[0])
    
    Result = pd.DataFrame(columns=["k", "Jaccard"])
    
    
    for k in range(1 , 100 , 5):
        
        selection_1 = copy1[copy1[0] == k]
        selection_2 = copy2[copy2[0] == k]
        
        for repetition in range(1, Repetitions + 1):
            
            repetition_pull_1 = selection_1.head(k)
            repetition_pull_2 = selection_2.head(k)
            
            repetition_pull_1.index = range(len(repetition_pull_1))
            repetition_pull_2.index = range(len(repetition_pull_2))
            
            # GO terms of the selection:
            
            # 1 Group:
            
            GO_Terms1 = pd.DataFrame(repetition_pull_1[2])
            GO_Terms1.columns = ["Term"]
            
            GO_Terms1['Term'] = GO_Terms1['Term'].apply(literal_eval)
            
            GO_terms_iter1 = []
    
            for index, rows in GO_Terms1.iterrows():
                
                GO_terms_iter1.extend([rows.Term])
                        
            GO_terms_iter_No_Blank1 = [x for x in GO_terms_iter1 if x != []]
            GO_terms_flat1 = [item for sublist in GO_terms_iter_No_Blank1 for item in sublist]
    
            # 2 Group:
            
            GO_Terms2 = pd.DataFrame(repetition_pull_2[2])
            GO_Terms2.columns = ["Term"]
            GO_Terms2['Term'] = GO_Terms2['Term'].apply(literal_eval)
            GO_terms_iter2 = []
    
            for index, rows in GO_Terms2.iterrows():
                
                GO_terms_iter2.extend([rows.Term])
                        
            GO_terms_iter_No_Blank2 = [x for x in GO_terms_iter2 if len(x) > 0]
            GO_terms_flat2 = [item for sublist in GO_terms_iter_No_Blank2 for item in sublist] 
            
            # Unique:
            
            GO_Array_1 = set(GO_terms_flat1)
            GO_Array_2 = set(GO_terms_flat2)
            
            # Save Result:
            
            # Control the division by 0:
            
            if(len(GO_Array_1) == 0 & len(GO_Array_2) == 0):
            
                Result_itera = pd.DataFrame({"k": [k],
                                             "Jaccard" : [None]})
                
                Result = Result.append(Result_itera)
                
            else:
                
                Result_itera = pd.DataFrame({"k": [k],
                                             "Jaccard" : [len(GO_Array_1.intersection(GO_Array_2)) / len(GO_Array_1.union(GO_Array_2))]})
                
                Result = Result.append(Result_itera)

            # Remove 
            
            selection_1 = selection_1.iloc[k:]
            selection_2 = selection_2.iloc[k:]
            
            # If not double plot discomment the following line:
            
    #Print_Jaccard(Result, Name_Bin, enrichment, save_directory)
    
    # For doble plot:
    
    return(Result)
    
        
    

    

def Analyze_GO_Terms(Terms_directory, 
                     save_directory, 
                     Name_Prob,
                     Name_Bin,
                     enrichment,
                     Repetitions = 10):
    
    
    Terms_GO = pd.read_csv(Terms_directory, sep = "\t")
    
    copy1 = pd.read_csv(save_directory + "_Rand_Information_"+ "_" + enrichment+ "_" + Name_Prob + ".txt",
                              sep = ",", header = None, skiprows=[0])
    
    copy2 = pd.read_csv(save_directory + "_Rand_Information_"+ "_" + enrichment+ "_" + Name_Bin + ".txt",
                              sep = ",", header = None, skiprows=[0])
    
    Result = pd.DataFrame(columns=["k", "Mean", "Network"])

    for k in range(1 , 100 , 5):
        
        selection_1 = copy1[copy1[0] == k]
        selection_2 = copy2[copy2[0] == k]
        
        for repetition in range(1, Repetitions + 1):
            
            repetition_pull_1 = selection_1.head(k)
            repetition_pull_2 = selection_2.head(k)
            
            repetition_pull_1.index = range(len(repetition_pull_1))
            repetition_pull_2.index = range(len(repetition_pull_2))
            
            # GO terms of the selection:
            
            # 1 Group:
            
            GO_Terms1 = pd.DataFrame(repetition_pull_1[2])
            GO_Terms1.columns = ["Term"]
            
            GO_Terms1['Term'] = GO_Terms1['Term'].apply(literal_eval)
            
            GO_terms_iter1 = []
    
            for index, rows in GO_Terms1.iterrows():
                
                GO_terms_iter1.extend([rows.Term])
                        
            GO_terms_iter_No_Blank1 = [x for x in GO_terms_iter1 if x != []]
            GO_terms_flat1 = [item for sublist in GO_terms_iter_No_Blank1 for item in sublist]
    
            # 2 Group:
            
            GO_Terms2 = pd.DataFrame(repetition_pull_2[2])
            GO_Terms2.columns = ["Term"]
            GO_Terms2['Term'] = GO_Terms2['Term'].apply(literal_eval)
            GO_terms_iter2 = []
    
            for index, rows in GO_Terms2.iterrows():
                
                GO_terms_iter2.extend([rows.Term])
                        
            GO_terms_iter_No_Blank2 = [x for x in GO_terms_iter2 if len(x) > 0]
            GO_terms_flat2 = [item for sublist in GO_terms_iter_No_Blank2 for item in sublist] 
            
            # Unique:
            
            GO_Array_1 = set(GO_terms_flat1)
            GO_Array_2 = set(GO_terms_flat2)
            
            # Intersection of terms:
            
            GO_Intersection = GO_Array_1.intersection(GO_Array_2)
            
            GO_Array_1 = pd.DataFrame(GO_Array_1)
            GO_Array_2 = pd.DataFrame(GO_Array_2)
            
            # Only in probabilities:
            
            GO_Array_1_only = GO_Array_1[~GO_Array_1.isin(GO_Intersection)]
            GO_Array_1_only = GO_Array_1_only.dropna()
            
            # Only in binary:
            
            GO_Array_2_only = GO_Array_2[~GO_Array_2.isin(GO_Intersection)]
            GO_Array_2_only = GO_Array_2_only.dropna()
            
            # Subset the terms:
            
            GO_Array_Intersection =  Terms_GO[Terms_GO.GO_ID.isin(GO_Intersection)].reset_index(drop=True)
            GO_Array_Intersection = GO_Array_Intersection.drop_duplicates(subset='GO_ID', keep="last")
            
            if len(GO_Array_1_only) == 0:
                
                GO_Array_1 = pd.DataFrame({"Level" : [0]})
                
            elif len(GO_Array_1_only) != 0:
            
                GO_Array_1 = Terms_GO[Terms_GO["GO_ID"].isin(GO_Array_1_only[0])].reset_index(drop=True)
                GO_Array_1 = GO_Array_1.drop_duplicates(subset='GO_ID', keep="last")
                
            if len(GO_Array_2_only) == 0:
                
                GO_Array_2 = pd.DataFrame({"Level" : [0]})
                
            elif len(GO_Array_2_only) != 0:
            
                GO_Array_2 = Terms_GO[Terms_GO["GO_ID"].isin(GO_Array_2_only[0])].reset_index(drop=True)
                GO_Array_2 = GO_Array_2.drop_duplicates(subset='GO_ID', keep="last")
            
            # Get the means:
            
            Intersection_mean = GO_Array_Intersection["Level"].mean()
            GO_Array_1_Value = GO_Array_1["Level"].mean()
            GO_Array_2_Value = GO_Array_2["Level"].mean() 
            
            # Control to avoid bugs:
            
            if GO_Array_1_Value == 0:
                
                GO_Array_1_Value = 0
                
            if GO_Array_2_Value == 0:
                
                GO_Array_2_Value = 0
                
            if np.isnan(Intersection_mean) == True:
                
                Intersection_mean = 0
            
            # Save the info in the result data frame:
            
            Result_itera = pd.DataFrame({"k": [k], 
                                        "Mean": [Intersection_mean] ,
                                        "Network" : ["Inter"]})
            
            Result = Result.append(Result_itera, ignore_index = True) 
            
            Result_itera = pd.DataFrame({"k": [k], 
                                        "Mean": [GO_Array_1_Value] ,
                                        "Network" : ["Prob"]})
    
            Result = Result.append(Result_itera,  ignore_index = True) 
            
            Result_itera = pd.DataFrame({"k": [k], 
                                        "Mean": [GO_Array_2_Value] ,
                                        "Network" : ["Bin"]})
    
            Result = Result.append(Result_itera,  ignore_index = True) 
            
            # Remove 
            
            selection_1 = selection_1.iloc[k:]
            selection_2 = selection_2.iloc[k:]  
    
    Print_GO_Terms_information(Result, Name_Bin, enrichment, save_directory )

def Print_GO_Terms_information(Result, Name_Bin, enrichment, save_directory):
       
    # Plot arguments (Explained in the CoEx module):
    
    Result.Mean = Result.Mean.astype(float)
    Result.k = Result.k.astype(float)
    
    Result = Result[Result.Mean != 0]
    fig, ax = plt.subplots(figsize=(6,3))
    ax.xaxis.grid(True, which='minor')

    ax.tick_params("both", labelsize= 15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)
    
    ax.xaxis.set_major_locator(MultipleLocator(10))
    
    #ax.axis([0, 100, 3.5, 5])
    
    ax.xaxis.set_tick_params(width=10)
    ax.yaxis.set_tick_params(width=10)
    
    ax.set_xlabel('# k',fontsize = 15) 
    ax.set_ylabel('# Level', fontsize = 15)
            
    ax = sns.lineplot(x="k", y="Mean", hue="Network", legend = False, data= Result, ax=ax, hue_order=['Prob', 'Bin', 'Inter'])
        
    ax.set_xlabel('# k',fontsize = 15) 
    ax.set_ylabel('# Levels', fontsize = 15)                                           
    
    fig.savefig(save_directory + "_GO_Levels_Natasa_PAPER" + Name_Bin + enrichment + ".png")
    

def Print_Results_Comparison(Network_Name_bin,
                             Network_Name_prob,
                             save_directory,option, comparion,
                             Enrichment = "BP"):
    
    # Files #:
    
    if Enrichment == "BP":
        
        x1 = 40
        y1 = 100
        x2 = 2
        y2 = 26
        x3 = 0
        y3 = 25
        
    if Enrichment == "CC":
        
        x1 = 50
        y1 = 100
        x2 = 10
        y2 = 38
        x3 = 10
        y3 = 40
    
    if Enrichment == "MF":
        
        x1 = 40
        y1 = 100
        x2 = 5
        y2 = 30
        x3 = 5
        y3 = 18
    

    # GO Enrichment:
    
    path_GO_prob = [save_directory + "_Enrichment_GO_" + str(Enrichment) + "_" + str(Network_Name_prob) + ".txt"]
    path_GO_bin = [save_directory + "_Enrichment_GO_" + str(Enrichment) + "_" + str(Network_Name_bin) + ".txt"]
    

    Prob_file_GO = pd.read_csv(path_GO_prob[0], sep = ",", header = 0)
    Other_file_GO = pd.read_csv(path_GO_bin[0], sep = ",", header = 0)
    
    # Clusters Enrichment:
       
    path_Clust_prob = [save_directory + "_Enrichment_Cluster_" + str(Enrichment) + "_" + str(Network_Name_prob) + ".txt"]
    path_Clust_bin = [save_directory + "_Enrichment_Cluster_" + str(Enrichment) + "_" + str(Network_Name_bin) + ".txt"]
    
    Prob_file_Cluster = pd.read_csv(path_Clust_prob[0], sep = ",", header = 0)
    Other_file_Cluster = pd.read_csv(path_Clust_bin[0], sep = ",", header = 0)
    
    # Genes Enrichment:
    
    path_gene_prob = [save_directory + "_Enrichment_Genes_" + str(Enrichment) + "_" + str(Network_Name_prob) + ".txt"]
    path_gene_bin = [save_directory + "_Enrichment_Genes_" + str(Enrichment) + "_" + str(Network_Name_bin) + ".txt"]
    
    Prob_file_Genes = pd.read_csv(path_gene_prob[0], sep = ",", header = 0)
    Other_file_Genes = pd.read_csv(path_gene_bin[0], sep = ",", header = 0)
    
    #  Printing #
    
    # Clusters Plots:
    
    # Data frame preparation for plot:
    
    Prob_Category = ["Prob"] * Prob_file_GO.size
    Other_Category = ["Bin"] * Other_file_GO.size
    
    Prob_file_Cluster['Network'] = pd.Series(Prob_Category)
    Other_file_Cluster['Network'] = pd.Series(Other_Category)

    frames = [Prob_file_Cluster, Other_file_Cluster]
    
    Data_Base_Cluster = pd.concat(frames)
    Data_Base_Cluster = Data_Base_Cluster[Data_Base_Cluster.Enrichment_Proportions != 0]
    Data_Base_Cluster.Enrichment_Proportions = Data_Base_Cluster.Enrichment_Proportions.astype(float)
    
    # Plot arguments (Explained in the CoEx module):
    
    fig, ax = plt.subplots(figsize=(6,3))
    ax.xaxis.grid(True, which='minor')
 
    ax.tick_params("both", labelsize= 15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)
    
    ax.xaxis.set_tick_params(width=10)
    ax.yaxis.set_tick_params(width=10)
    
        
    ax.axis([0, 100, x1, y1]) # Cluster Y xs
    
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.set_xlabel('# k',fontsize = 10) 
    ax.set_ylabel('# %', fontsize = 10)
    
    ax = sns.lineplot(x="K_Number", y="Enrichment_Proportions", hue="Network" , legend = False, data=Data_Base_Cluster, ax=ax)
    
    ax.set_xlabel('# k',fontsize = 15) 
    ax.set_ylabel('# %', fontsize = 15)
    
    fig.savefig(save_directory + "_Line_Clusters_plot_PAPER" + str(option) + Network_Name_bin + comparison + ".png")

    # GO Plots:

    # Data frame preparation for plot:
    
    Prob_Category_GO = ["Prob"] * Prob_file_GO.size
    Other_Category_GO = ["Bin"] * Other_file_GO.size
    
    Prob_file_GO['Network'] = pd.Series(Prob_Category_GO)
    Other_file_GO['Network'] = pd.Series(Other_Category_GO)
    
    frames_GO = [Prob_file_GO, Other_file_GO]
    
    Data_Base_GO = pd.concat(frames_GO)
    Data_Base_GO = Data_Base_GO[Data_Base_GO.GO_Proportions != 0]
    
    Data_Base_GO.GO_Proportions = Data_Base_GO.GO_Proportions.astype(float)
    
    # Plot arguments (Explained in the CoEx module):    
    
    fig, ax = plt.subplots(figsize=(6,3))
    ax.xaxis.grid(True, which='minor')
 
    ax.tick_params("both", labelsize= 15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)
    
    ax.xaxis.set_tick_params(width=10)
    ax.yaxis.set_tick_params(width=10)
    
    ax.axis([0, 100, x2, y2]) # GO axes
    
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.set_xlabel('# k',fontsize = 10) 
    ax.set_ylabel('# %', fontsize = 10)
    
    ax = sns.lineplot(x="K_Number", y="GO_Proportions", hue="Network" , legend = False, data=Data_Base_GO, ax=ax)
    
    ax.set_xlabel('# k',fontsize = 15) 
    ax.set_ylabel('# %', fontsize = 15)
    
    fig.savefig(save_directory + "_GO_Clusters_plot_PAPER" + str(option) + Network_Name_bin + comparison + ".png")

    # Gene Plots:    

    # Data frame preparation for plot:
    
    Prob_Category = ["Prob"] * Prob_file_Genes.size
    Other_Category = ["Bin"] * Other_file_Genes.size
    
    Prob_file_Genes['Network'] = pd.Series(Prob_Category)
    Other_file_Genes['Network'] = pd.Series(Other_Category)
    
    frames = [Prob_file_Genes, Other_file_Genes]
    
    Data_Base_Genes = pd.concat(frames)
    Data_Base_Genes= Data_Base_Genes[Data_Base_Genes.Total_enriched != 0]
    Data_Base_Genes.Total_enriched = Data_Base_Genes.Total_enriched.astype(float)
    
    # Plot arguments (Explained in the CoEx module): 

    fig, ax = plt.subplots(figsize=(6,3))
    ax.xaxis.grid(True, which='minor')
 
    ax.tick_params("both", labelsize= 15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)
    
    ax.xaxis.set_tick_params(width=10)
    ax.yaxis.set_tick_params(width=10)
    
    ax.axis([0, 100, x3, y3]) # GO axes
    
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.set_xlabel('# k',fontsize = 10) 
    ax.set_ylabel('# %', fontsize = 10)
    
    ax = sns.lineplot(x="K_Option", y="Total_enriched", hue="Network" ,legend = False, data=Data_Base_Genes, ax=ax)
    
    ax.set_xlabel('# k',fontsize = 15) 
    ax.set_ylabel('# %', fontsize = 15)
    
    fig.savefig(save_directory + "_Genes_plot_PAPER" + str(option) + Network_Name_bin + comparison + ".png")
    
def Print_Jaccard(Matrix, Name_Bin, enrichment, save_directory):
    
    Jaccard_Data = Matrix.dropna()

    # Plot arguments (Explained in the CoEx module):
    
    fig, ax = plt.subplots(figsize=(6,3))
    ax.xaxis.grid(True, which='minor')
 
    ax.tick_params("both", labelsize= 15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)
    
    ax.xaxis.set_tick_params(width=10)
    ax.yaxis.set_tick_params(width=10)
    
    ax.axis([0, 100, 0.4, 0.73])
    
    ax.xaxis.set_major_locator(MultipleLocator(10))
    

    ax.set_xlabel('# k',fontsize = 15) 
    ax.set_ylabel('# Index', fontsize = 15)
    
    ax = sns.lineplot(x="k", y="Jaccard", data= Jaccard_Data, ax=ax)
        
    ax.set_xlabel('# k',fontsize = 15) 
    ax.set_ylabel('# Index', fontsize = 15)
    
    fig.savefig(save_directory + "_Jaccard_PAPER" + Name_Bin + enrichment + ".png")


def Print_Jaccard_Two(Matrix_1, Matrix_2, Matrix_3, Name_Bin,  enrichment, save_directory,  option = "2"):
    
    Jaccard_Data_1 = Matrix_1.dropna()
    Jaccard_Data_2 = Matrix_2.dropna()
    Jaccard_Data_3 = Matrix_3.dropna()
    
    # Labels:
    
    if option == "2":
        
        Category_1 = ["Network_1"] * Jaccard_Data_1.size
        Category_2 = ["Network_2"] * Jaccard_Data_2.size
        
        Jaccard_Data_1['Group'] = pd.Series(Category_1)
        Jaccard_Data_2['Group'] = pd.Series(Category_2)
        
    
        frames = Jaccard_Data_1.append(Jaccard_Data_2)
        
        fig, ax = plt.subplots(figsize=(6,3))
        ax.xaxis.grid(True, which='minor')
     
        ax.tick_params("both", labelsize= 15)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(2.5)
        ax.spines['bottom'].set_linewidth(2.5)
        
        ax.xaxis.set_tick_params(width=10)
        ax.yaxis.set_tick_params(width=10)
        
        ax.axis([0, 100, 0.5, 0.8])
        
        ax.xaxis.set_major_locator(MultipleLocator(10))
        
    
        ax.set_xlabel('# k',fontsize = 15) 
        ax.set_ylabel('# Index', fontsize = 15)
                      
        
        ax = sns.lineplot(x="k", y="Jaccard", hue = "Group",legend = False, 
                          data= frames, ax=ax, hue_order = ["Network_1", "Network_2"], palette = ["b", "r"])
            
        ax.set_xlabel('# k',fontsize = 15) 
        ax.set_ylabel('# Index', fontsize = 15)
        
        fig.savefig(save_directory + "_Jaccard_PAPER_Both" + Name_Bin + enrichment + ".png")
        
    if option == "3":
        
        Category_1 = ["Network_1"] * Jaccard_Data_1.size
        Category_2 = ["Network_2"] * Jaccard_Data_2.size
        Category_3 = ["Network_3"] * Jaccard_Data_3.size
        
        Jaccard_Data_1['Group'] = pd.Series(Category_1)
        Jaccard_Data_2['Group'] = pd.Series(Category_2)
        Jaccard_Data_3['Group'] = pd.Series(Category_3)
        
        
        frames = Jaccard_Data_1.append(Jaccard_Data_2)
        frames = frames.append(Jaccard_Data_3)
        
        fig, ax = plt.subplots(figsize=(6,3))
        ax.xaxis.grid(True, which='minor')
     
        ax.tick_params("both", labelsize= 15)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(2.5)
        ax.spines['bottom'].set_linewidth(2.5)
        
        ax.xaxis.set_tick_params(width=10)
        ax.yaxis.set_tick_params(width=10)
        
        ax.axis([0, 100, 0.4, 0.8])
        
        ax.xaxis.set_major_locator(MultipleLocator(10))
        
    
        ax.set_xlabel('# k',fontsize = 15) 
        ax.set_ylabel('# Index', fontsize = 15)
                        
        ax = sns.lineplot(x="k", y="Jaccard", hue = "Group",legend = False, 
                          data= frames, ax=ax, hue_order = ["Network_1", "Network_2", "Network_3"], palette = ["b", "r", "g"])
            
        ax.set_xlabel('# k',fontsize = 15) 
        ax.set_ylabel('# Index', fontsize = 15)
        
        fig.savefig(save_directory + "_Jaccard_PAPER_Both" + Name_Bin + enrichment + ".png")
        

    # Plot arguments (Explained in the CoEx module):
    


# MAIN:
    
save_directory = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/GI_Results/low/"

Name_Network_Original = "Network_Low_Confidence_Named.txt"
Name_Network_Prob = ['Work_Network' + "_" + "Prob_" + Name_Network_Original]
Name_Network_Bin = ['Work_Network' + "_" + "Bin_" + Name_Network_Original]

path_BP = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Python_Srcipts/Final_Results/Annotation Script/Biological Process/Table_Yeast2Goat_BF.csv" 
path_CC = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_BP_Levels.csv"
path_MF = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Enrich_Annot/Entrez_Annotation_BP_Levels.csv"


Print_Results_Comparison(Name_Network_Bin[0],
                             Name_Network_Prob[0],
                             save_directory,
                             Enrichment = "BP",
                             option = "Low_Confidence_BP",
                             comparison = "Bin_VS_Prob")

Print_Results_Comparison(Name_Network_Bin[0],
                             Name_Network_Prob[0],
                             save_directory,
                             Enrichment = "CC",
                             option = "Low_Confidence_CC",
                             comparison = "Bin_VS_Prob")

Print_Results_Comparison(Name_Network_Bin[0],
                             Name_Network_Prob[0],
                             save_directory,
                             Enrichment = "MF",
                             option = "Low_Confidence_MF",
                             comparison = "Bin_VS_Prob")


To_Jaccard_Index(Name_Network_Prob[0], 
                                       Name_Network_Bin[0],
                                       10, 
                                       66,
                                       save_directory,
                                       "BP")
To_Jaccard_Index(Name_Network_Prob[0], 
                                       Name_Network_Bin[0],
                                       10, 
                                       66,
                                       save_directory,
                                       "CC")
To_Jaccard_Index(Name_Network_Prob[0], 
                                       Name_Network_Bin[0],
                                       10, 
                                       66,
                                       save_directory,
                                       "MF")


# BP

Terms_directory = path_BP

Analyze_GO_Terms(Terms_directory,
                                       save_directory,
                                       Name_Network_Prob[0],
                                       Name_Network_Bin[0],
                                       "BP",
                                       Repetitions = 10)

# CC

Terms_directory =  path_CC


Analyze_GO_Terms(Terms_directory,
                                       save_directory,
                                       Name_Network_Prob[0],
                                       Name_Network_Bin[0],
                                       "CC",
                                       Repetitions = 10)


# MF

Terms_directory = path_MF

Analyze_GO_Terms(Terms_directory,
                                       save_directory,
                                       Name_Network_Prob[0],
                                       Name_Network_Bin[0],
                                       "MF",
                                       Repetitions = 10)




# Jaccard index for various plots (if we run this code make sure the comment of the plots is done in the function)


path_BP = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Human_PPI/Annotation_Human_Ensembl/Table_Enrich_Ensembl_BP.csv" 
path_CC = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Human_PPI/Annotation_Human_Ensembl/Table_Enrich_Ensembl_CC.csv"
path_MF = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Human_PPI/Annotation_Human_Ensembl/Table_Enrich_Ensembl_MF.csv"

# Network 1:

save_directory_1 = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Data_Rank_Method/Threshold_Correct/High/"
Name_Network_Original_1 = "Coex_Net_Original_High_Rank_1%.txt"
Name_Network_Prob_1 = ['Work_Network' + "_" + "Prob_" + Name_Network_Original_1]
Name_Network_Bin_1 = ['Work_Network' + "_" + "Bin_" + Name_Network_Original_1]

Matrix_1_BP = To_Jaccard_Index(Name_Network_Prob_1[0], Name_Network_Bin_1[0], 10, 66, save_directory_1, 
                            "BP")
Matrix_1_CC = To_Jaccard_Index(Name_Network_Prob_1[0], Name_Network_Bin_1[0], 10, 66, save_directory_1, 
                            "CC")
Matrix_1_MF = To_Jaccard_Index(Name_Network_Prob_1[0], Name_Network_Bin_1[0], 10, 66, save_directory_1, 
                            "MF")

# Network 2:

save_directory_2 = "//home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Data_Rank_Method/Threshold_Correct/Medium/"
Name_Network_Original_2 = "Coex_Net_Original_Medium_Rank_5%.txt"
Name_Network_Prob_2 = ['Work_Network' + "_" + "Prob_" + Name_Network_Original_2]
Name_Network_Bin_2 = ['Work_Network' + "_" + "Bin_" + Name_Network_Original_2]

Matrix_2_BP = To_Jaccard_Index(Name_Network_Prob_2[0], Name_Network_Bin_2[0], 10, 66, save_directory_2, 
                            "BP")

Matrix_2_CC = To_Jaccard_Index(Name_Network_Prob_2[0], Name_Network_Bin_2[0], 10, 66, save_directory_2, 
                            "CC")

Matrix_2_MF = To_Jaccard_Index(Name_Network_Prob_2[0], Name_Network_Bin_2[0], 10, 66, save_directory_2, 
                            "MF")

# Network 3:

save_directory_3 = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/CoEx_Results/Data_Rank_Method/Threshold_Correct/Low/"
Name_Network_Original_3 = "Coex_Net_Original_Low_Rank_10%.txt"
Name_Network_Prob_3 = ['Work_Network' + "_" + "Prob_" + Name_Network_Original_3]
Name_Network_Bin_3 = ['Work_Network' + "_" + "Bin_" + Name_Network_Original_3]

Matrix_3_BP = To_Jaccard_Index(Name_Network_Prob_3[0], Name_Network_Bin_3[0], 10, 66, save_directory_3, 
                            "BP")
Matrix_3_CC = To_Jaccard_Index(Name_Network_Prob_3[0], Name_Network_Bin_3[0], 10, 66, save_directory_3, 
                            "CC")
Matrix_3_MF = To_Jaccard_Index(Name_Network_Prob_3[0], Name_Network_Bin_3[0], 10, 66, save_directory_3, 
                            "MF")

# Plot:

Print_Jaccard_Two(Matrix_1_BP, Matrix_2_BP, Matrix_3_BP, Name_Network_Original_1, "BP", save_directory_1, "3")
Print_Jaccard_Two(Matrix_1_CC, Matrix_2_CC, Matrix_3_CC,  Name_Network_Original_1, "CC", save_directory_1, "3")
Print_Jaccard_Two(Matrix_1_MF, Matrix_2_MF, Matrix_3_MF,  Name_Network_Original_1, "MF", save_directory_1, "3")
