'''
This Script is for the last part of the protocol. Here we get the final correlation Matrices and calculate the distances.
With this distances we can perform a lot of different analyses
'''

from itertools import product

from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import precision_recall_curve, auc, roc_curve, roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os

from scipy.stats import mannwhitneyu
from statsmodels.distributions.empirical_distribution import ECDF

# Function To compare Distributions (Uniform Vs PPI):
'''
Finish this comparisons with all the distributions
'''
def Compare_Distributions(save_path):
    
    PPI = pd.read_csv("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Real_Distributions/Empirical_PPI/Work_Network_Prob_Low_PPI.txt",
                                sep = " ", header = None)
    CoEx = pd.read_csv("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Real_Distributions/Empirical_CoEx/Work_Network_Prob_Coex_Net_Original_Medium_Rank_5%.txt",
                                sep = " ", header = None)
    GI = pd.read_csv("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Real_Distributions/Empirical_GI/Work_Network_Prob_Network_Medium_Confidence_Named.txt",
                                sep = " ", header = None)
    Uniform = pd.read_csv("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Work_Network_Prob_UniformEta_0.0350076103500761_Mu_0.27_5000_Repetition_0_0.003_ER_.txt",
                                sep = " ", header = None)
    Beta_1 = pd.read_csv("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Work_Network_Prob_BetaEta_0.001_Mu_0.9_5000_Repetition_0_0.003_ER_.txt",
                                sep = " ", header = None)
    Beta_2 = pd.read_csv("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Real_Distributions/Beta_Mean_0.27_Var_0.0069/Work_Network_Prob_BetaEta_0.0350076103500761_Mu_0.27_5000_Repetition_0_0.003_ER_.txt",
                                sep = " ", header = None)
    Beta_3 = pd.read_csv("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Real_Distributions/Beta_Mean_0.78_Var_0.0028/Work_Network_Prob_BetaEta_0.01631701631701632_Mu_0.78_5000_Repetition_0_0.003_ER_.txt",
                                sep = " ", header = None)
    
    ### Table to get Means and Variances: ###
    
    # Probabilities Distributions:
    
    PPI        =  PPI.iloc[:, 2]
    CoEx       =  CoEx.iloc[:, 2]
    GI         =  GI.iloc[:, 2]
    Uniform    =  Uniform.iloc[:, 2]
    Beta_1     =  Beta_1.iloc[:, 2]
    Beta_2     =  Beta_2.iloc[:, 2]
    Beta_3     =  Beta_3.iloc[:, 2]
    
    # Means:
    
    means_final = []
    
    means_final +=  [np.mean(PPI)]
    means_final +=  [np.mean(CoEx)]
    means_final += [np.mean(GI)]
    means_final += [np.mean(Uniform)]
    means_final += [np.mean(Beta_1)]
    means_final += [np.mean(Beta_2)]
    means_final += [np.mean(Beta_3)]
    
    # Variance:
    
    variance_final = []
    
    variance_final += [np.var(PPI)]
    variance_final += [np.var(CoEx)]
    variance_final += [np.var(GI)]
    variance_final += [np.var(Uniform)]
    variance_final += [np.var(Beta_1)]
    variance_final += [np.var(Beta_2)]
    variance_final += [np.var(Beta_3)]
    
    # Networks:
    
    Distributions = ["PPI", "CoEx", "GI", "Uniform", "Beta_1", "Beta_2", "Beta_3"]
    
    # Final dataframe and save it:
    
    Final_Mean_Variance_db = pd.DataFrame({"Dist"      : Distributions,
                                            "Mean"     : means_final  ,
                                            "Variance" : variance_final})
    
    Final_Mean_Variance_db.to_csv(save_path + "Table_Mean_Variance.txt", index = None, sep = " ")
    

    ### Table to compare Distributions with Mann Whitney test: ####
    
    Options1 = [0,1,2,3,4,5,6]
    Options2 = [0,1,2,3,4,5,6]
    
    Result_P_values = pd.DataFrame(index=np.arange(7), columns=np.arange(7))
    
    for option1, option2 in product(Options1, Options2):
        
        # For First Comparison:
        
        if option1 == 0:
            test1 = PPI
            name = "PPI"
        if option1 == 1:
            test1 = CoEx
            name = "CoEx"
        if option1 == 2:
            test1 = GI
            name = "GI"
        if option1 == 3:
            test1 = Uniform
            name = "Uniform"
        if option1 == 4:
            test1 = Beta_1
            name = "Beta_1"
        if option1 == 5:
            test1 = Beta_2
            name = "Beta_2"
        if option1 == 6:
            test1 = Beta_3
            name = "Beta_3"
        
        # For second Comparison:

        if option2 == 0:
            test2 = PPI
            name2 = "PPI"
        if option2 == 1:
            test2 = CoEx
            name2 = "CoEx"
        if option2 == 2:
            test2 = GI
            name2 = "GI"
        if option2 == 3:
            test2 = Uniform
            name2 = "Uniform"
        if option2 == 4:
            test2 = Beta_1
            name2 = "Beta_1"
        if option2 == 5:
            test2 = Beta_2
            name2 = "Beta_2"
        if option2 == 6:
            test2 = Beta_3
            name2 = "Beta_2"
        
        _, p_value= mannwhitneyu(test1, test2, use_continuity=True, alternative='two-sided')
        
        Result_P_values.iloc[option1,option2] = p_value
        
        # Plot comparisons:
        
        cdf1 = ECDF(test1)
        cdf2 = ECDF(test2)

        x = np.linspace(0,1,2**10)

        y1 = cdf1(x)
        y2 = cdf2(x)

        plt.plot(x,y1)
        plt.plot(x,y2)
        
        plt.savefig(save_path + "Distribution_Cumulative_Comparisons_" + str(name) + "_" + str(name2) + ".png")
        
        
    
    Result_P_values.columns = ["PPI", "CoEx", "GI", "Uniform", "Beta_1", "Beta_2", "Beta_3"]
    Result_P_values.index = ["PPI", "CoEx", "GI", "Uniform", "Beta_1", "Beta_2", "Beta_3"]
        
    Result_P_values.to_csv(save_path + "Table_ManWhitsen.txt", sep = " ")
        

# Function to prepare the Matrices:

def Prepare_The_Matrices(Option1, Option2, Option3, empirical = False):
    
    # Read File:
    
    file_Bin1 = pd.read_csv(Option1 + "Final_Correlation_Bin.txt")
    file_Prob1 = pd.read_csv(Option1 + "Final_Correlation_Prob.txt")
    
    file_Bin2 = pd.read_csv(Option2 + "Final_Correlation_Bin.txt")
    file_Prob2 = pd.read_csv(Option2 + "Final_Correlation_Prob.txt")
        
    file_Bin3 = pd.read_csv(Option3 + "Final_Correlation_Bin.txt")
    file_Prob3 = pd.read_csv(Option3 + "Final_Correlation_Prob.txt")
    
    # Categories Distribution:
    
    if empirical == False:
    
        Group_Bin = file_Bin["Variance"].map(str) + file_Bin["Mean"].map(str)
        Group_Prob = file_Prob["Variance"].map(str) + file_Prob["Mean"].map(str)
        Final_Cat_Bin_Dist = Group_Bin.reset_index()
        Final_Cat_Prob_Dist = Group_Prob.reset_index()
        
    else:
        
        # Model 1:
        
        file_Bin1   = file_Bin1.drop_duplicates()
        file_Prob1  = file_Prob1.drop_duplicates()
        Group_Bin1  = file_Bin1.Real.fillna("Uniform")
        Group_Prob1 = file_Prob1.Real.fillna("Uniform")
        Group_Bin1  = Group_Bin1.reset_index(drop = True)
        Group_Prob1 = Group_Prob1.reset_index(drop = True)
        file_Bin1   = file_Bin1.reset_index(drop = True)
        file_Prob1  = file_Prob1.reset_index(drop = True)
        
        # Model 2:
        
        file_Bin2   = file_Bin2.drop_duplicates()
        file_Prob2  = file_Prob2.drop_duplicates()
        Group_Bin2  = file_Bin2.Real.fillna("Uniform")
        Group_Prob2 = file_Prob2.Real.fillna("Uniform")
        Group_Bin2  = Group_Bin2.reset_index(drop = True)
        Group_Prob2 = Group_Prob2.reset_index(drop = True)
        file_Bin2   = file_Bin2.reset_index(drop = True)
        file_Prob2  = file_Prob2.reset_index(drop = True)
        
        # Model 3:
        
        file_Bin3   = file_Bin3.drop_duplicates()
        file_Prob3  = file_Prob3.drop_duplicates()
        Group_Bin3  = file_Bin3.Real.fillna("Uniform")
        Group_Prob3 = file_Prob3.Real.fillna("Uniform")
        Group_Bin3  = Group_Bin3.reset_index(drop = True)
        Group_Prob3 = Group_Prob3.reset_index(drop = True)
        file_Bin3   = file_Bin3.reset_index(drop = True)
        file_Prob3  = file_Prob3.reset_index(drop = True)
        
        # Put together:
        
        file_Bin  = file_Bin1.append(file_Bin2)
        file_Bin  = file_Bin.append(file_Bin3)
        file_Prob = file_Prob1.append(file_Prob2)
        file_Prob = file_Prob.append(file_Prob3)
        
        Group_Bin = Group_Bin1.append(Group_Bin2)
        Group_Bin = Group_Bin.append(Group_Bin3)
        Group_Prob = Group_Prob1.append(Group_Prob2)
        Group_Prob = Group_Prob.append(Group_Prob3)
        

    # Add New Category Distribution:
    
    if empirical == False:
    
        Final_Cat_Bin_Dist.columns = ["T", "G"]
        Final_Cat_Prob_Dist.columns = ["T", "G"]
        file_Bin["Density"] = Final_Cat_Bin_Dist["G"]
        file_Prob["Density"] = Final_Cat_Prob_Dist["G"]
    
        # Delete non used Variables:
    
        file_Bin = file_Bin.drop(['Nodes', 'Densities', 'Distribution','Repetitions', 'Variance', 'Mean'], axis=1)
        file_Prob = file_Prob.drop(['Nodes', 'Densities', 'Distribution','Repetitions', 'Variance', 'Mean'], axis=1)
    
        # Relabel Distributions:

        Number2Labels = {
        
                "0.00.0"  : "Unif",
                "0.99900000000000010.5": "Var : 0.999 Mean : 0.5",
                "0.99900000000000010.1": "Var : 0.999 Mean : 0.1",
                "0.99900000000000010.9": "Var : 0.999 Mean : 0.9",
                "0.50.5"  : "Var : 0.5 Mean : 0.5",
                "0.50.1"  : "Var : 0.5 Mean : 0.1",
                "0.50.9"  : "Var : 0.5 Mean : 0.9",
                "0.0010.5"  : "Var : 0.001 Mean : 0.5",
                "0.0010.1"  : "Var : 0.001 Mean : 0.1",
                "0.0010.9"  : "Var : 0.001 Mean : 0.9",
        
            }
     
        file_Bin.Density = file_Bin.Density.apply(Number2Labels.get)
        file_Prob.Density = file_Prob.Density.apply(Number2Labels.get)
        
    else:
        
        file_Bin = file_Bin.drop(['Nodes', 'Densities', 'Distribution','Repetitions', 'Variance', 'Mean', 'Real'], axis=1)
        file_Prob = file_Prob.drop(['Nodes', 'Densities', 'Distribution','Repetitions', 'Variance', 'Mean', 'Real'], axis=1)
        file_Bin["Density"] = Group_Bin
        file_Prob["Density"] = Group_Prob
        
    # Return Prepared Information:
    
    return(file_Bin, file_Prob)
    
# Funtion To Calculate the Distances:
       
def Similarity_Infor_for_One_Model(file, Option, Network, empirical = False, group = "ALL"):
    
    # Save the information in a different Variable Before Distance:
    
    Info_Categories_Dens= file["Network"] + file["Density"]
    
    # Drop The Variables:
    
    file_Distance = file.drop(["Network", "Density"], axis=1)
    
    # Calculate the distance:
    
    distances = pdist(file_Distance.values, metric='euclidean')
    dist_matrix = squareform(distances)
    Data_Frame_Distance = pd.DataFrame(dist_matrix)
    
    # Go to similarity Matrix (Normalize):
    
    Data_Frame_Similarity = 1 - Data_Frame_Distance
    Data_Frame_Similarity = Data_Frame_Similarity / Data_Frame_Similarity.max()[0]
    Data_Frame_Similarity_np = np.array(Data_Frame_Similarity)
    
    # Save Similarity Matrix:
    
    Data_Frame_Similarity.to_csv( Option + "Correlation_Distance"+ Network + ".txt", header = True, index = False)
    
    # Y_Score (without diagonal):
    
    y_score = Data_Frame_Similarity_np[np.triu_indices_from(Data_Frame_Similarity_np, k = 1)]
    
    # Group Labels (Y_True):
    
    if empirical == False:
    
        Number2Numbers= {
        
                "Unif": 1,
                "Var : 0.999 Mean : 0.5" : 2,
                "Var : 0.999 Mean : 0.1" : 3,
                "Var : 0.999 Mean : 0.9" : 4,
                "Var : 0.5 Mean : 0.5" : 5,
                "Var : 0.5 Mean : 0.1" : 6,
                "Var : 0.5 Mean : 0.9" : 7,
                "Var : 0.001 Mean : 0.5" : 8,
                "Var : 0.001 Mean : 0.1" : 9,
                "Var : 0.001 Mean : 0.9" : 10,
                }
    
        y_true = Info_Categories_Dens
        y_true = Info_Categories_Dens.apply(Number2Numbers.get)
        my_y_true = ~pdist(np.array(y_true).reshape(-1,1), 'cityblock').astype(bool)
        
    else:
        
        if group == "ALL":
            
            Number2Numbers= {
                    "_ER_Uniform" : 2,
                    "_ER_PPI"     : 2,
                    "_ER_GI"      : 3,
                    "_ER_CoEx"    : 4,
                    "_BA_Uniform" : 5,
                    "_BA_PPI"     : 5,
                    "_BA_GI"      : 6,
                    "_BA_CoEx"    : 7,
                    "_RG_Uniform" : 8,
                    "_RG_PPI"     : 8,
                    "_RG_GI"      : 10,
                    "_RG_CoEx"    : 11
                    }
        
            y_true = Info_Categories_Dens
            y_true = Info_Categories_Dens.apply(Number2Numbers.get)
            my_y_true = ~pdist(np.array(y_true).reshape(-1,1), 'cityblock').astype(bool)
            
        if group == "Density":
            
            Number2Numbers= {
                    "Uniform" : 1,
                    "GI"      : 2,
                    "CoEx"    : 3,
                    "PPI"     : 1
                    }
            
            y_true = file["Density"]
            y_true = file["Density"].apply(Number2Numbers.get)
            my_y_true = ~pdist(np.array(y_true).reshape(-1,1), 'cityblock').astype(bool)
            
        if group == "Topology":
            
            Number2Numbers= {
                    "_ER_"   : 1,
                    "_BA_"   : 2,
                    "_RG_"   : 3,
                    }
            
            y_true = file["Network"]
            y_true = file["Network"].apply(Number2Numbers.get)
            my_y_true = ~pdist(np.array(y_true).reshape(-1,1), 'cityblock').astype(bool)
            
    # Return the information:
    
    return(y_score, my_y_true)
    
# Function to Make the Plots:
    
def Precision_Recall_plot(x_bin,
                          x_prob,
                          y_bin,
                          y_prob,
                          Option,
                          save_directory):
    
    # Groups:
    
    Bin_group = pd.Series(['Bin']).repeat(len(x_bin))
    Prob_group = pd.Series(['Prob']).repeat(len(x_prob))
    
    Final_Groups = Bin_group.append(Prob_group)
    
    # Data frame:
    
    precision_Bin = pd.Series(y_bin)
    precision_Prob = pd.Series(y_prob)
    precision_Bin = precision_Bin.append(precision_Prob)
    
    recall_Bin = pd.Series(x_bin)
    recall_Prob = pd.Series(x_prob)
    recall_Bin = recall_Bin.append(recall_Prob)
    
    Final_Groups = Final_Groups.reset_index(drop = True)
    recall_Bin = recall_Bin.reset_index(drop = True)
    precision_Bin = precision_Bin.reset_index(drop = True)
    
    Plot_df = pd.DataFrame()
    Plot_df["Recall"] = recall_Bin
    Plot_df["Precision"] = precision_Bin
    Plot_df["Group"] = Final_Groups
    
    Plot_df = Plot_df.reset_index(drop=True)
    
    # Plot:
     
    fig, ax = plt.subplots(figsize=(12,8))
    plt.grid(0)
    paper_rc = {'lines.linewidth': 1, 'lines.markersize': 10}                  
    sns.set_context("paper", rc = paper_rc) 
    
    sns.lineplot(x = "Recall", y = "Precision", hue = "Group", data = Plot_df )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)
    ax.xaxis.set_tick_params(width=2.5)
    ax.yaxis.set_tick_params(width=2.5)
    
    plt.figtext(.5,.9,'Precision-Recall', fontsize=30, ha='center')
    ax.set_xlabel('# Recall',fontsize = 20) 
    ax.set_ylabel('# Precision', fontsize = 20)
                            
    fig.savefig(save_directory + "_Correlation_Plot_Prob_" + Option + ".png")
    
# ROC Curve Plot:
    
def Plot_ROC(fpr_bin, tpr_bin,
             fpr_prob, tpr_prob,
             save_path, Model):
    
    # Groups:

    Bin_group = pd.Series(['Bin']).repeat(len(fpr_bin))
    Prob_group = pd.Series(['Prob']).repeat(len(fpr_prob))
    
    Final_Groups = Bin_group.append(Prob_group) 
    Final_Groups = Final_Groups.reset_index(drop = True)
    
    # DataFrame:
    
    fpr_prob = pd.Series(fpr_prob) 
    fpr_bin = pd.Series(fpr_bin)
    tpr_prob = pd.Series(tpr_prob) 
    tpr_bin = pd.Series(tpr_bin)
    
    
    Final_Data = pd.DataFrame()
    Final_Data["fpr"] = fpr_bin.append(fpr_prob)
    Final_Data["tpr"] = tpr_bin.append(tpr_prob)
    Final_Data = Final_Data.reset_index(drop = True)
    Final_Data["Network"] = Final_Groups
    
    Final_Data_Bin = Final_Data[Final_Data.Network == "Bin"]
    Final_Data_Prob = Final_Data[Final_Data.Network == "Prob"]
    
    fig, ax = plt.subplots(figsize=(12,8))
    lw = 2
    plt.plot(Final_Data_Bin["fpr"], Final_Data_Bin["tpr"], c = "orange",
         lw=lw)
    plt.plot(Final_Data_Prob["fpr"], Final_Data_Prob["tpr"], c = "green",
         lw=lw)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([-0.05, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=20)
    plt.ylabel('True Positive Rate', fontsize=20)
    plt.title('ROC', fontsize=20)
    plt.legend(["Bin", "Prob", "Random"], loc="lower right", prop={'size': 15})
    plt.savefig(save_path + "ROC_Curve_From_" + Model + ".png")
    plt.show()
   
# Main:   
    
# Work Directory:

os.chdir("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Real_Distributions/Considering_Random/")
save_directory = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Real_Distributions/Considering_Random/"

# Preparing the Data: 

Bin_Corr, Prob_Corr = Prepare_The_Matrices(Option1 = "ER_1/",
                                           Option2 = "BA_1/",
                                           Option3 = "RG_1/",
                                           empirical = True)
    
# Getting the rest of information:

## GROUP = ALL ##

# Getting the variables for the ROC Curve:

y_score_Bin, y_true_Bin   = Similarity_Infor_for_One_Model(Bin_Corr,  Option = "ALL", Network = "Bin", empirical = True)   
y_score_Prob, y_true_Prob = Similarity_Infor_for_One_Model(Prob_Corr, Option = "ALL", Network = "Prob", empirical = True)  

# Precision-Recall:

precision_Bin, recall_Bin, thresholds_Bin = precision_recall_curve(y_true_Bin, y_score_Bin)
precision_Prob, recall_Prob, thresholds_Prob = precision_recall_curve(y_true_Prob, y_score_Prob)

Precision_Recall_plot( recall_Bin,    recall_Prob, 
                       precision_Bin, precision_Prob, 
                       Option = "ALL",    save_directory = save_directory)

auc(recall_Bin, precision_Bin)   # 0.3673  Precision-Recall
auc(recall_Prob, precision_Prob) # 0.8711  Precision-Recall

# ROC Curve:

fpr_Bin, tpr_ER_Bin, thresholds_ER_Bin = roc_curve(y_true_Bin, y_score_Bin)
fpr_Prob, tpr_ER_Prob, thresholds_ER_Prob = roc_curve(y_true_Prob, y_score_Prob)

roc_auc_score(y_true_Bin, y_score_Bin)   # 0.8807  Area under ROC
roc_auc_score(y_true_Prob, y_score_Prob) # 0.9803  Area Under ROC

Plot_ROC(fpr_Bin, tpr_ER_Bin, fpr_Prob, tpr_ER_Prob, save_directory, Model= "ALL")

## GROUP = Density ##

# Getting the variables for the ROC Curve:

y_score_Bin, y_true_Bin   = Similarity_Infor_for_One_Model(Bin_Corr,  Option = "ALL", Network = "Bin", empirical = True, group = "Density")   
y_score_Prob, y_true_Prob = Similarity_Infor_for_One_Model(Prob_Corr, Option = "ALL", Network = "Prob", empirical = True, group = "Density")  

# Precision-Recall:

precision_Bin, recall_Bin, thresholds_Bin = precision_recall_curve(y_true_Bin, y_score_Bin)
precision_Prob, recall_Prob, thresholds_Prob = precision_recall_curve(y_true_Prob, y_score_Prob)

Precision_Recall_plot( recall_Bin,    recall_Prob, 
                       precision_Bin, precision_Prob, 
                       Option = "ALL_Group_Density",    save_directory = save_directory)

auc(recall_Bin, precision_Bin)   # 0.3711  Precision-Recall
auc(recall_Prob, precision_Prob) # 0.5363  Precision-Recall

# ROC Curve:

fpr_Bin, tpr_ER_Bin, thresholds_ER_Bin = roc_curve(y_true_Bin, y_score_Bin)
fpr_Prob, tpr_ER_Prob, thresholds_ER_Prob = roc_curve(y_true_Prob, y_score_Prob)

roc_auc_score(y_true_Bin, y_score_Bin)   # 0.4981 Area under ROC
roc_auc_score(y_true_Prob, y_score_Prob) # 0.5424 Area Under ROC

Plot_ROC(fpr_Bin, tpr_ER_Bin, fpr_Prob, tpr_ER_Prob, save_directory, Model= "ALL_Group_Density")

## GROUP = Topology ##

# Getting the variables for the ROC Curve:

y_score_Bin, y_true_Bin   = Similarity_Infor_for_One_Model(Bin_Corr,  Option = "ALL", Network = "Bin", empirical = True, group = "Topology")   
y_score_Prob, y_true_Prob = Similarity_Infor_for_One_Model(Prob_Corr, Option = "ALL", Network = "Prob", empirical = True, group = "Topology")  

# Precision-Recall:

precision_Bin, recall_Bin, thresholds_Bin = precision_recall_curve(y_true_Bin, y_score_Bin)
precision_Prob, recall_Prob, thresholds_Prob = precision_recall_curve(y_true_Prob, y_score_Prob)

Precision_Recall_plot( recall_Bin,    recall_Prob, 
                       precision_Bin, precision_Prob, 
                       Option = "ALL_Group_Topology",    save_directory = save_directory)

auc(recall_Bin, precision_Bin)   # 1  Precision-Recall
auc(recall_Prob, precision_Prob) # 0.999  Precision-Recall

# ROC Curve:

fpr_Bin, tpr_ER_Bin, thresholds_ER_Bin = roc_curve(y_true_Bin, y_score_Bin)
fpr_Prob, tpr_ER_Prob, thresholds_ER_Prob = roc_curve(y_true_Prob, y_score_Prob)

roc_auc_score(y_true_Bin, y_score_Bin)   # 1 Area under ROC
roc_auc_score(y_true_Prob, y_score_Prob) # 1 Area Under ROC

Plot_ROC(fpr_Bin, tpr_ER_Bin, fpr_Prob, tpr_ER_Prob, save_directory, Model= "ALL_Group_Topology")












