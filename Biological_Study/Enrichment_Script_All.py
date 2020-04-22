# ENRICHMENT SCRIPT

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

from ast import literal_eval 


# DESCRIPTION.

'''
This script contains all the functions which are necessary to perfrom an enrichment analyses
'''

# Function to calculate the percentage of Clusters with at least one GO term enriched:

def Function_Calculate_Clust_Perc(Data_frame, total_K):
    
    Result =  pd.DataFrame(columns=["K_Number", "Enrichment_Proportions"])
    
    for k in range(1,total_K,5):
        
        container = Data_frame[Data_frame.Option == (k)]
        result_it = container[container.Enriched_GO != 0]
        result_it = len(result_it.index) * 100 / (k)
        
        if result_it != 0:
            
            Iterator_DB = pd.DataFrame({'K_Number' : [k], 
                                        'Enrichment_Proportions' : [result_it]})
            Result = Result.append(Iterator_DB)
    
        else:
            
            Iterator_DB = pd.DataFrame({'K_Number' : [k], 
                                        'Enrichment_Proportions' : [0]})

            Result = Result.append(Iterator_DB)
        
    return(Result)
    
# Function to Calculate Percentage GO terms enriched in each k option:

def Function_GO_enriched(Data_frame, total_GO_Terms, total_K):
    
    Result = pd.DataFrame(columns=["K_Number", "GO_Proportions"])
    
    for k in range(1,total_K,5):
        
        container = Data_frame[Data_frame.Option == (k)]
        
        if container.size != 0:
            
            GO_terms = pd.DataFrame(container["Term"])
            GO_terms_iter = []
            for index, rows in GO_terms.iterrows(): 
                GO_terms_iter.extend([rows.Term])
                
            GO_terms_iter_No_Blank = [x for x in GO_terms_iter if x != []]
            GO_terms_flat = [item for sublist in GO_terms_iter_No_Blank for item in sublist]
            
            Total_GO = len(Counter(GO_terms_flat).keys())*100/total_GO_Terms

            Iterator = pd.DataFrame({'K_Number' : [k], 
                                       'GO_Proportions' : [Total_GO]})   
            Result = Result.append(Iterator)
            
        elif container.size == 0:
            Iterator = pd.DataFrame({'K_Number' : [k], 
                                        'GO_Proportions' : [0]})   
            Result = Result.append(Iterator)
            
    return(Result)
    
# Function to calculate the genes that are enriched per each cluster
    
def Function_Gene_enriched(Data_frame_Genes, Data_frame_Results, total_K, Gene_Names):
    
    total_genes = Gene_Names.size
    genes_enriched = pd.DataFrame(columns=["K_Option", "Total_enriched"])
    
    for k in range(1, total_K, 5):
        
        container_genes = Data_frame_Genes[Data_frame_Genes.Option == (k)]
        container_terms_enriched = Data_frame_Results[Data_frame_Results.Option == (k)]
        
        Count = 0
        
        for cluster in range(1, k):
            
            genes_cluster = container_genes[container_genes.Cluster == cluster]
            GO_cluster = container_terms_enriched[container_terms_enriched.Cluster == cluster]
            
            if GO_cluster["Term"].size == 0:
                
                Count = Count + 0
                
            elif GO_cluster["Term"].size != 0:
        
            # Take the GO terms (to flat variable):
            
                GO_cluster =  pd.DataFrame(GO_cluster.Term)
                GO_terms_iter = []
                
                for index, rows in GO_cluster.iterrows(): 
                    GO_terms_iter.extend([rows.Term])
                GO_terms_iter_No_Blank = [x for x in GO_terms_iter if x != []]
                GO_terms_flat = [item for sublist in GO_terms_iter_No_Blank for item in sublist]
            
            # Transform the data into a data frame:
            
                GO_erniched_cluster_k = pd.DataFrame({"GO_Term" : GO_terms_flat})
            
            # Merge to count the number:
            
                genes = genes_cluster[genes_cluster.GO_Term.isin(GO_erniched_cluster_k["GO_Term"])]
                
            # Count
            
            Count = Count + genes["Gene"].nunique()
            
        Total = ((Count * 100)/total_genes)
            
        Result = pd.DataFrame({"K_Option": k, "Total_enriched": [Total]})
        genes_enriched = genes_enriched.append(Result)
    
    return(genes_enriched)
            
# Function to choose the cluster method:
    
def Cluster_Option(option, genes, Distance,
                   method_clust = "Kmedoids"):
    
    if method_clust == "Kmedoids":
        
            medoids_random= random.choices(list(range(len(genes))), k= option)
            medoids_result = kmedoids(np.array(Distance), medoids_random , data_type = "distance_matrix")
            medoids_result = medoids_result.process()
            clusters_i = medoids_result.get_clusters()
        
            return clusters_i
            
    else:
            
        print("No more avaliable distances yet")
        
# Function to perform the enrichment:
      
def Enrichment_Analyses_GO_terms(Name_Network,  
                                 save_directory,
                                 Annotation_Directory,
                                 Original_Network_Name,
                                 enrichment = "NO", 
                                 MaxSize = 500, 
                                 MinSize = 5,
                                  Repetitions = 10, 
                                  total_K = 100, 
                                  Rand_K_Compare = 66,
                                  comparison = "Bin_VS_Prob"
                                ):
    
    ## Preparation of the Network to Enrich ##
    
    # All the information of the data base:
    
    GO_Complete_Experimental = pd.read_csv(Annotation_Directory, sep = "\t",
                              header = 0)
    
    GO_Complete_Experimental.drop("Level", inplace = True, axis = 1) 
    GO_Complete_Experimental.columns = ["Gene", "GO_Term"] 
    
    # The genes of our Network:
    
    genes_Network = pd.read_csv(save_directory + "_Gene_Names_" + str(Original_Network_Name) ,
                              sep = ",", header = None, skiprows=[0])
    genes_Network = pd.DataFrame({"Gene" : genes_Network[1]})
  
    # Annotation of the genes of our Network:
    
    GO_Join_Network = GO_Complete_Experimental[GO_Complete_Experimental.Gene.isin(genes_Network.Gene)]
    
    # Filters:
    
    Cut_By = pd.DataFrame(GO_Join_Network.groupby('GO_Term')['Gene'].nunique(dropna=True))
    Cut_By = Cut_By[Cut_By.Gene > MinSize]
    Cut_By = Cut_By[Cut_By.Gene < MaxSize] 
    
    GO_Join_Network = GO_Join_Network[GO_Join_Network.GO_Term.isin(Cut_By.index)]
    
    ## Distances Preparation for clustering ##
    
    # Loading:
    
    Distance = pd.read_csv(save_directory + "_Result_Tijana_Final_" + Name_Network,
                                   sep=" ", header = 0, )
    Distance.set_index('1', inplace=True)

    ## Variables to store the results ##

    # Result DataFrames:
    
    Results_GO_Enrichment_Final = pd.DataFrame(columns=["K_Option", "Terms_Enriched"])
    Results_Cluster_Enrichment_Final = pd.DataFrame(columns=["K_Option", "Cluster_Enriched"])
    Results_Genes_Percent_Final = pd.DataFrame(columns=["K_Option", "Total_enriched"])
    
    # Results for Rand_Index:
    
    Results_Rand_Index = pd.DataFrame(columns=["Option","Cluster", "Terms"])
    
    # Results for Gene:
    
    Results_Gene_GO_final = pd.DataFrame(columns=["Gene", "GO_Term", "Cluster", "Option", "Repetition"]) 
    
    ## Enrichment Analyses ##:
           
    for Statistics in range(Repetitions):
        
        print("Repetition number",Statistics)
        
        # Per each repetition we should reload the variables:
        
        Results_Enrichment = pd.DataFrame(columns=["Option", "Enriched_GO",
                                                      "Cluster", "Num_Genes_Annotated", "Term"]) 
    
        Results_Gene_GO = pd.DataFrame(columns=["Gene", "GO_Term", "Cluster", "Option","Repetition"]) 

     
        for option in range(1,total_K,5):
        
            # First we start with the computing of the clusters
            
            medoids_ini = option
            
            print("K number", option)
            
            clusters = Cluster_Option(option, genes_Network, Distance,
                                                        method_clust = "Kmedoids")
            
            clusters_Data_frame = pd.DataFrame(columns=["Gene", "Cluster", "Option"])

            for i in range(len(clusters)):
                
                gene_selection = genes_Network.iloc[clusters[i]]
                cluster_Repetition = np.repeat(i + 1, len(gene_selection))
                Option_Repetition = np.repeat(option, len(gene_selection))
                Iterator_DB = pd.DataFrame({'Gene' : gene_selection["Gene"], 'Cluster' : cluster_Repetition, 
                                                "Option" : Option_Repetition })
                clusters_Data_frame = clusters_Data_frame.append(Iterator_DB )
                
            
            # Number of genes annotated and how many of them are in each category:
            
            Total_Annotated_Genes = GO_Join_Network["Gene"].nunique()
            Number_Genes_per_GO = GO_Join_Network.groupby('GO_Term')['Gene'].nunique(dropna=True)
            
            # With this information we go to the cluster:
            
            # For each cluster:
            
            for cluster in range(option):
                
                # Cluster Selection:
                
                selection_cluster = clusters_Data_frame[clusters_Data_frame.Cluster == (cluster + 1)]

                # Annotation of genes in the cluster with GO:
                
                GO_Selected = GO_Join_Network[GO_Join_Network.Gene.isin(selection_cluster.Gene)]

                # Put inside of the external variable to keep info:
                
                Genes_itera = pd.DataFrame({"Gene": GO_Selected.Gene,
                                                         "GO_Term":GO_Selected.GO_Term,
                                                         'Cluster':cluster + 1,
                                                         "Option":option,
                                                         "Repetition":Statistics})

                
                Results_Gene_GO = Results_Gene_GO.append(Genes_itera)

                # Number Genes with a concrete GO term in the cluster:
                
                k_selection = GO_Selected.groupby('GO_Term')['Gene'].nunique()
                
                K_and_k_data_frame = pd.merge(Number_Genes_per_GO, k_selection , 
                                    on='GO_Term', how='right')

                # Total genes with annotation in the cluster:

                Total_Annotated_Cluster = GO_Selected["Gene"].nunique()

                # For the results of the enrichment in the cluster:
                
                Results_Enrichment_Cluster = pd.DataFrame(columns=["GO_Term", "p_value", 'Cluster', 'Option']) 

                # For each GO term in the cluster:
                
                # Probabilistic:
                
                for GO_term in range(len(K_and_k_data_frame)):
                    
                    Enrich_GO = K_and_k_data_frame.iloc[GO_term]
                    M =  Total_Annotated_Genes 
                    k = Enrich_GO[0] 
                    N = Total_Annotated_Cluster
                    X = Enrich_GO[1]
        
                    p_value = hypergeom.sf(X - 1, M, k, N) 
                    
                    results_db_iter = pd.DataFrame({'GO_Term' : K_and_k_data_frame.index[GO_term], 
                                                        'p_value':[p_value], 'Cluster' : [cluster + 1], 
                                                        'Option' : [medoids_ini]})
    
                    Results_Enrichment_Cluster = Results_Enrichment_Cluster.append(results_db_iter)
                    
                
                if len(K_and_k_data_frame) !=0: # To avoid errors with empty clusters
        
                    
                    p_value_Correction = multipletests(Results_Enrichment_Cluster["p_value"] ,
                                                   alpha=0.01, method='fdr_bh',  is_sorted=False, 
                                                   returnsorted=False)
                    
                    Results_Enrichment_Cluster['p_value_Correction'] = p_value_Correction[1]
                   
                    count_GO_Enriched = sum(Results_Enrichment_Cluster['p_value_Correction'] < 0.05)

                    names_GO_Enriched = Results_Enrichment_Cluster[Results_Enrichment_Cluster.p_value_Correction < 0.05]
                    names_GO_Enriched = list(names_GO_Enriched["GO_Term"])
                    
                    Results_itera_Enrich = pd.DataFrame({'Option':[medoids_ini], "Enriched_GO":[count_GO_Enriched] , 
                                                         'Cluster':[cluster + 1], "Num_Genes_Annotated": [Total_Annotated_Cluster],
                                                         "Term" : [names_GO_Enriched]}) 
    
                    Results_Enrichment = Results_Enrichment.append(Results_itera_Enrich)
                    
                elif len(K_and_k_data_frame) ==0:
                    
                   Results_itera_Enrich = pd.DataFrame({'Option':[medoids_ini], "Enriched_GO":[0] , 
                                                         'Cluster':[cluster + 1], "Num_Genes_Annotated": [Total_Annotated_Cluster],
                                                         "Term": [[]]})
    
                   Results_Enrichment = Results_Enrichment.append(Results_itera_Enrich)
                    
             
                # Calculations per each iteration:
                
        GO_percent = Function_GO_enriched(Results_Enrichment, GO_Join_Network["GO_Term"].nunique(), total_K)
        Cluster_percent = Function_Calculate_Clust_Perc(Results_Enrichment, total_K)
        Genes_percent = Function_Gene_enriched(Results_Gene_GO,Results_Enrichment, total_K, genes_Network)
        
        Results_GO_Enrichment_Final = Results_GO_Enrichment_Final.append(GO_percent)
        Results_Cluster_Enrichment_Final = Results_Cluster_Enrichment_Final.append(Cluster_percent)
        Results_Genes_Percent_Final = Results_Genes_Percent_Final.append(Genes_percent)
        
        # For Rand:
        
        Rand = pd.DataFrame({"Option" : Results_Enrichment.Option,
                             "Cluster" : Results_Enrichment.Cluster,
                             "Terms" : Results_Enrichment.Term})
    
        Results_Rand_Index = Results_Rand_Index.append(Rand)
        
        # For Cluster.
        
        Results_Gene_GO_final = Results_Gene_GO_final.append(Results_Gene_GO)
        
    # Only a concrete k value to compare with The Rand index:
    
    #Results_Rand_Index = Results_Rand_Index[Results_Rand_Index.Option == Rand_K_Compare] 

    Results_GO_Enrichment_Final.drop("K_Option", inplace = True, axis = 1) 
    Results_GO_Enrichment_Final.drop("Terms_Enriched", inplace = True, axis = 1)              
    Results_Cluster_Enrichment_Final.drop("Cluster_Enriched", inplace = True, axis = 1)
    Results_Cluster_Enrichment_Final.drop("K_Option", inplace = True, axis = 1)

    
    # Save files just in case:
    
    Results_GO_Enrichment_Final.to_csv(save_directory + "_Enrichment_GO"+ "_" + enrichment+ "_" +Name_Network + ".txt", header = True, index = False)
    Results_Cluster_Enrichment_Final.to_csv(save_directory + "_Enrichment_Cluster"  + "_"+ enrichment+ "_" + Name_Network+ ".txt", header = True, index = False)
    Results_Genes_Percent_Final.to_csv(save_directory + "_Enrichment_Genes" + "_" +enrichment+ "_"  + Name_Network+ ".txt", header = True, index = False)
    
    # Save info about Genes and Clusters:
    
    Results_Gene_GO_final.to_csv(save_directory + "_Cluster_Information_"+ "_" + enrichment+ "_" +Name_Network + ".txt", header = True, index = False)
    
    # Save Rand:
    
    Results_Rand_Index.to_csv(save_directory + "_Rand_Information_"+ "_" + enrichment+ "_" +Name_Network + ".txt", header = True, index = False)
         
# Descriptives #
    
'''
Here we introduce two different statistics which can help us to know more about our clustering.
The Rand index compares the clusters (how similar they are, in other words, if the methods are
dividing the data in the same way). On the other hand the Jaccard index allows as to know if the
information which is captured is the same or not (more about biological information).
'''    

# Function to compute the Rand index: #
    
# First creation of contigency:
    
def To_Matrix_Rand(Name_Prob, 
                   Name_Bin,
                   Repetitions, 
                   Rand_K_Compare,
                   save_directory,
                   enrichment):
    
    Result = pd.DataFrame(np.zeros((Rand_K_Compare + 1, Rand_K_Compare + 1)))
    
    copy1 = pd.read_csv(save_directory + "_Rand_Information_"+ "_" + enrichment+ "_" + Name_Prob[0] + ".txt",
                              sep = ",", header = None, skiprows=[0])
    
    copy2 = pd.read_csv(save_directory + "_Rand_Information_"+ "_" + enrichment+ "_" + Name_Bin[0] + ".txt",
                              sep = ",", header = None, skiprows=[0])
    
    for repetition in range(1, Repetitions + 1):
        
        repetition_pull_1 = copy1.head(Rand_K_Compare)
        repetition_pull_2 = copy2.head(Rand_K_Compare)
        
        for cluster in range(1, len(repetition_pull_1) + 1):
            
            GO_terms_1 = repetition_pull_1.Terms[repetition_pull_1.Cluster == cluster]
            
            for cluster2 in range(1, len(repetition_pull_2) + 1):
                
                GO_terms_2 = repetition_pull_2.Terms[repetition_pull_2.Cluster == cluster2]
                
                if repetition == 1:

                    Result.loc[cluster, cluster2] = (Result.loc[cluster, cluster2] + 
                              len((set(GO_terms_1[0]) & set(GO_terms_2[0]))))
                
                elif repetition > 1:
               
                    Result.loc[cluster, cluster2] = (Result.loc[cluster, cluster2] + 
                              len((set(GO_terms_1[0]) & set(GO_terms_2[0])))/2)
                          
        copy1 = copy1.iloc[Rand_K_Compare:]
        copy2 = copy2.iloc[Rand_K_Compare:]
        
    Result = Result.drop(0)
    Result = Result.drop(0, axis=1)
    
    # Here we calculate the true positive and all the things necessary for the Rand index
    
    vComb = np.vectorize(myComb)
    
    tp,fp,tn,fn =  get_tp_fp_tn_fn(np.array(Result), vComb)
    
    # Give the results
    
    print("Rand index: %f" % (float(tp + tn) / (tp + fp + fn + tn)))
            
# Function to calculate the true positives and negatives for Rand Index:
    
def myComb(a,b):
    
  return comb(a,b,exact=True)
    
def get_tp_fp_tn_fn(cooccurrence_matrix, vComb):
    
  tp_plus_fp = vComb(cooccurrence_matrix.sum(0, dtype=int),2).sum()
  tp_plus_fn = vComb(cooccurrence_matrix.sum(1, dtype=int),2).sum()
  tp = vComb(cooccurrence_matrix.astype(int), 2).sum()
  fp = tp_plus_fp - tp
  fn = tp_plus_fn - tp
  tn = comb(cooccurrence_matrix.sum(), 2) - tp - fp - fn

  return [tp, fp, tn, fn]
        
# Function to compute the Jaccard Index: #    

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
            
    Print_Jaccard(Result, Name_Bin, enrichment, save_directory)
    
    
# Print Jaccard Index:
            
def Print_Jaccard(Matrix, Name_Bin, enrichment, save_directory):
    
    Jaccard_Data = Matrix.dropna()

    # Plot arguments (Explained in the CoEx module):
    
    fig, ax = plt.subplots(figsize=(12,8))
    plt.grid(0)

    ax.tick_params("both", labelsize= 24)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)
    ax.xaxis.set_tick_params(width=2.5)
    ax.yaxis.set_tick_params(width=2.5)
    
    plt.figtext(.5,.9,'Jaccard Index ' + enrichment, fontsize=30, ha='center')
    ax.set_xlabel('# k',fontsize = 24) 
    ax.set_ylabel('# Index', fontsize = 24)
    
    ax = sns.lineplot(x="k", y="Jaccard", data= Jaccard_Data, ax=ax)
        
    ax.set_xlabel('# k',fontsize = 24) 
    ax.set_ylabel('# Index', fontsize = 24)
    
    fig.savefig(save_directory + "_Jaccard_" + Name_Bin + enrichment + ".png")
        
        
def To_Matrix_Rand_2(Name_Prob, 
                   Name_Bin,
                   Repetitions, 
                   Rand_K_Compare,
                   save_directory,
                   enrichment):
       
    copy1 = pd.read_csv(save_directory + "_Rand_Information_"+ "_" + enrichment+ "_" + Name_Prob + ".txt",
                              sep = ",", header = None, skiprows=[0])
    
    copy2 = pd.read_csv(save_directory + "_Rand_Information_"+ "_" + enrichment+ "_" + Name_Bin + ".txt",
                              sep = ",", header = None, skiprows=[0])
    
    copy1[2] = copy1[2].apply(literal_eval)
    copy2[2] = copy2[2].apply(literal_eval)
    
    Result_Final = pd.DataFrame(columns=["k", "Rand"])
    
    for k in range(1 , 10 , 5):
        
        selection_1 = copy1[copy1[0] == k]
        selection_2 = copy2[copy2[0] == k]
        
        for repetition in range(1, Repetitions + 1):
            
            Result = pd.DataFrame(np.zeros((k + 1, k + 1)))
        
            repetition_pull_1 = selection_1.head(k)
            repetition_pull_2 = selection_2.head(k)
            
            repetition_pull_1.columns = ["1", "Cluster","Terms"]
            repetition_pull_2.columns = ["1", "Cluster","Terms"]
            
            for cluster in range(1, len(repetition_pull_1) + 1):
                
                GO_terms_1, = repetition_pull_1.Terms[repetition_pull_1.Cluster == cluster]
                
                for cluster2 in range(1, len(repetition_pull_2) + 1):
                    
                    GO_terms_2, = repetition_pull_2.Terms[repetition_pull_2.Cluster == cluster2]
                    
                    Result.loc[cluster, cluster2] = len((set(GO_terms_1) & 
                              set(GO_terms_2)))
                
            selection_1 = selection_1.iloc[k:]
            selection_2 = selection_2.iloc[k:]
            
            Result = Result.drop(0, axis=1)
            Result = Result.drop(0, axis=0)
            
            vComb = np.vectorize(myComb)
            tp,fp,tn,fn =  get_tp_fp_tn_fn(np.array(Result), vComb)
        
            Rand = (float(tp + tn) / (tp + fp + fn + tn))
        
            Result_itera = pd.DataFrame({"k": [k],
                                             "Rand" : [Rand]})
        
            Result_Final = Result_Final.append(Result_itera)
 
           
Terms_directory = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Human_PPI/Annotation_Human_Ensembl/Table_Enrich_Ensembl_BP.csv"
save_directory = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Human_PPI/High/"
Name_Prob = "Work_Network_Prob_Human_High_PPI.txt"
Name_Bin = "Work_Network_Bin_Human_High_PPI.txt"
enrichment = "BP"
Repetitions = 10

        
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
    
    #Result = Result[Result.Network != "Inter"] If we dont want to include the interjection
    
    fig, ax = plt.subplots(figsize=(12,8))
    plt.grid(0)

    ax.tick_params("both", labelsize= 24)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)
    ax.xaxis.set_tick_params(width=2.5)
    ax.yaxis.set_tick_params(width=2.5)
    
    plt.figtext(.5,.9,'Level GO ' + enrichment, fontsize=30, ha='center')
    ax.set_xlabel('# k',fontsize = 24) 
    ax.set_ylabel('# Level', fontsize = 24)
            
    ax = sns.lineplot(x="k", y="Mean", hue="Network", data= Result, ax=ax, hue_order=['Prob', 'Bin', 'Inter'])
        
    ax.set_xlabel('# k',fontsize = 24) 
    ax.set_ylabel('# Levels', fontsize = 24)
                  
    #tes_Prob = Result[Result.Network == "Prob"]
    #tes_Bin = Result[Result.Network == "Bin"]
    #tes_inter = Result[Result.Network == "Inter"]
    
    #mean_Prob = ["Prob " + str(round(np.mean(tes_Prob["Mean"]) , 3))]
    #mean_Bin =  ["Bin " + str(round(np.mean(tes_Bin["Mean"]),3))]
    #mean_inter =  ["Inter " + str(round(np.mean(tes_inter["Mean"]),3))]
                  
    #new_labels = [mean_Prob[0], mean_Bin[0], mean_inter[0]]
    #[new_labels[2], new_labels[1], new_labels[0]]
                                 
    ax.legend(loc='upper left', shadow = True, fontsize = 15, markerfirst=False)
    
    fig.savefig(save_directory + "_GO_Levels_Natasa_2" + Name_Bin + enrichment + ".png")
    
    

    

