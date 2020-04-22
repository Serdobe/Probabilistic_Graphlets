'''
This Script Contains all the functions which are necessary to apply the correlation Distance
from a GDV matrix.
'''
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


from scipy.spatial.distance import squareform

# Function to Prepare The Matrix GDV:

def Prepare_Matrix_GDV_Corr(Network, save_directory):
    
    # Read GDV:
    
    prob_file = pd.read_csv(save_directory + "GDV_Work_Network_Prob_" + Network)
    bin_file  = pd.read_csv(save_directory + "GDV_Work_Network_Bin_" + Network)
    
    prob_file = prob_file.drop("NODE", axis = 1)
    bin_file  = bin_file.drop("NODE" , axis = 1)
    
    # Noise vector is necessary to Avoid Problems with the Spearman Correlation
    
    Noise_Vector = pd.DataFrame({"0" : [1],
                    "1" : [1], 
                    "2" : [1] ,
                    "3" : [1] ,
                    "4" : [1] ,
                    "5" : [1] ,
                    "6" : [1] ,
                    "7" : [1] ,
                    "8" : [1] ,
                    "9" : [1] ,
                    "10": [1] ,
                    "11": [1] ,
                    "12": [1] ,
                    "13": [1], 
                    "14": [1]})
    
    prob_file = prob_file.append(Noise_Vector)
    bin_file = bin_file.append(Noise_Vector)
    
    # Calculate the Correlation and Save the Matrix:
    
    GDV_corr_Prob = Spearman_Correlation(prob_file) # Spearman Correlation Probabilistic
    GDV_corr_Bin = Spearman_Correlation(bin_file) # Spearman Correlation Binary
    
    GDV_corr_Prob.to_csv(save_directory + "GDV_Correlation_Prob_" + str(Network), header = True, index = False)
    GDV_corr_Bin.to_csv(save_directory + "GDV_Correlation_Bin_" + str(Network), header = True, index = False)
    
    print("Done")
       
# Function to calculate Spearman Correlation from GDV:
    
def Spearman_Correlation(GDV):

    GDV_corr = GDV.rank().corr(method="pearson") # Spearman (because is the rank)
    
    return(GDV_corr)
    
# Function to get the upper triangle (is not ordered check if this is necessary or not):
    
def upper_tri_masking(A):
    m = A.shape[0]
    r = np.arange(m)
    mask = r[:,None] < r
    return A[mask]
    
# Function to create the final Matrix:
    
def Final_Matrix_Corr(Network, save_directory):
    
    file_prob = pd.read_csv(save_directory + "GDV_Correlation_Prob_" + str(Network), header = 0)
    file_bin = pd.read_csv(save_directory  + "GDV_Correlation_Bin_"  + str(Network), header = 0)
    
    diagonal_prob = np.fill_diagonal(file_prob.values,0)
    diagonal_bin = np.fill_diagonal(file_bin.values,0)
    
    diagonal_prob = squareform(file_prob)
    diagonal_bin =  squareform(file_bin)
    
    return(pd.DataFrame(diagonal_prob).transpose() ,pd.DataFrame(diagonal_bin).transpose())
    
    
# Function to Plot the Correlatio_Matrix.
    
Network = "EmpiricalCoExEta_0_Mu_0_5000_Repetition_14_0.003_RG_.txt" # bin RG
Network = "EmpiricalGIEta_0_Mu_0_5000_Repetition_20_0.003_RG_.txt"   # GI  RG
Network = "EmpiricalCoExEta_0_Mu_0_5000_Repetition_20_0.003_RG_.txt" # CoEx RG
Network = "EmpiricalPPIEta_0_Mu_0_5000_Repetition_20_0.003_ER_.txt" # ER Bin


save_directory = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Real_Distributions/Considering_Random/ER_1/"

def Correlation_Matrix_Print(Network, save_directory):
    
    file = pd.read_csv(save_directory + "GDV_Correlation_Prob_" + str(Network), header = 0)
    file = pd.read_csv(save_directory + "GDV_Correlation_Bin_" + str(Network), header = 0)
    
    # Probabilistics
    
    fig, ax = plt.subplots(figsize=(12,8))
    plt.grid(0)
    
    sns.set(font_scale=1.8)
    ax = sns.clustermap(file, 
            xticklabels=file.columns,
            yticklabels=file.columns,
            row_cluster = True, col_cluster = True)
    
    fig.savefig(save_directory + "_Correlation_Plot_PAPER_ER_Bin" + Network + ".png")
    
    
    

    ax.tick_params("both", labelsize= 20)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)
    ax.xaxis.set_tick_params(width=2.5)
    ax.yaxis.set_tick_params(width=2.5)
    
    plt.figtext(.5,.9,'Heat Map GCM Prob', fontsize=30, ha='center')
    ax.set_xlabel('# Orbits',fontsize = 20) 
    ax.set_ylabel('# Orbits', fontsize = 20)
                  
    fig.savefig(save_directory + "_Correlation_Plot_Prob_" + Network + ".png")
                  
    # Binary
    
    fig2, ax2 = plt.subplots(figsize=(12,8))
    plt.grid(0)
    
    ax2 = sns.clustermap(file_bin, 
            xticklabels=file_bin.columns,
            yticklabels=file_bin.columns)
    
    fig2.savefig(save_directory + "_Correlation_Plot_Bin_" + Network + ".png")

    ax2.tick_params("both", labelsize= 20)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_linewidth(2.5)
    ax2.spines['bottom'].set_linewidth(2.5)
    ax2.xaxis.set_tick_params(width=2.5)
    ax2.yaxis.set_tick_params(width=2.5)
    
    plt.figtext(.5,.9,'Heat Map GCM Bin', fontsize=30, ha='center')
    ax2.set_xlabel('# Orbits',fontsize = 20) 
    ax2.set_ylabel('# Orbits', fontsize = 20)
                  
    fig2.savefig(save_directory + "_Correlation_Plot_Bin_" + Network + ".png")
         
    
    









