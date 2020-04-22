from itertools import product

import pandas as pd
import Creation_Network_Script
import Correlation_Distance_Functions
import Script_Parallel_Nodes_Count_Dist

'''
At this point , this code is prepared to compare two different things (which wast its main purpose)
the first is empirical vs varios beta. In this case eliminate the emp_option from the loops is necessary.
If you want to compare real distributions, then you have to add it. The distributions chosen as representative
for CoEx, GI, and PPI are described in te Creation_Network_Script.py
'''

# For Empirical Vs Various Beta Comparison

Distribution = ["Unif", "Beta"]
Network = ["_RG_"]
Nodes = [5000]
Densities = [0.091] # Take care of RG we have to change it
Repetitions = range(0, 50) 

ETA_MU = [[0.78, 0.01631701631701632], [0.27, 0.0350076103500761], [0.513, 0.3202164663312399]] # For the last one we only care about the shape
emp_option = ["PPI", "GI", "CoEx"]

# For Empirical Vs Uniform Comparison:


Network = ["_ER_"]
Distribution = ["Unif", "Empirical"]
Nodes = [5000]
Densities = [0.003]
Repetitions = range(0, 4) # 50 of each distribution at the end (each time create a Network)
ETA = [0]
MU = [0]
emp_option = ["PPI", "GI", "CoEx"]

# Save Directory:

save_directory = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/"

# Generation of the synthetics networks 

for repet, nodes, density, (mu,eta),  Emp_option in product(Repetitions, Nodes, Densities, ETA_MU, emp_option):
    
    
    Network_ER = Creation_Network_Script.get_Hyper_Geometric(nodes, density)
   
    Network_Name = ["Eta_" + str(eta) + "_Mu_" + str(mu)+ "_" +str(nodes) + "_Repetition_" + str(repet) + "_" + str(density)]
    
    Network_Name1 = [Network_Name[0] + "_RG_.txt"]

    Creation_Network_Script.Apply_Prob(Network_ER, save_directory, Network_Name1[0], eta, mu, empirical = False, option =  Emp_option)

  
# DO All the counts: 
# Here we use the same pypline as exposed in the biological analyses section.

file_Directory = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/"
save_directory = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/"
Run_path = "/home/sergio/workspace_Eclipse/Lucky_GOGO_Node_Parallel/Release/"
    
for repet, nodes, density, distribution, network, (mu,eta), Emp_option in product(Repetitions, Nodes, Densities, Distribution, Network, ETA_MU,emp_option):
    
    if distribution == "Unif":
        option = "Uniform" 
    if distribution == "Beta":
        option = "Beta"
    if distribution == "Empirical":
        option = "Empirical" + Emp_option
 
    Network_Name = [str(option)+"Eta_" + str(eta) + "_Mu_" + str(mu)+  "_"  + str(nodes) + "_Repetition_" + str(repet) + "_" + 
                    str(density) + str(network) + ".txt"]
    
    
    Name_Network_Prob = ["Work_Network_Prob_" + Network_Name[0]]
    Name_Network_Bin = ["Work_Network_Bin_" + Network_Name[0]]
    
    Network_Name2 = ["Eta_" + str(eta) + "_Mu_" + str(mu) + "_"  +str(nodes) + "_Repetition_" + str(repet) + "_" + 
                str(density) + str(network) + ".txt"]
    
    total_nodes = (len(open(save_directory + "_" + "Gene_Names" + Network_Name2[0]).readlines(  ))) - 1 
    
    Script_Parallel_Nodes_Count_Dist.Compute_Orbits_Spliting_Nodes(Name_Network_Bin[0], total_nodes,
                                                          file_Directory, save_directory,
                                                          Run_path, threads = 4)

    Script_Parallel_Nodes_Count_Dist.Compute_Orbits_Spliting_Nodes(Name_Network_Prob[0], total_nodes,
                                                          file_Directory, save_directory,
                                                          Run_path, threads = 4)
    
    GDV_Prob = Script_Parallel_Nodes_Count_Dist.Prepare_Matrix_GDV(Name_Network_Prob[0], save_directory)
    GDV_Bin = Script_Parallel_Nodes_Count_Dist.Prepare_Matrix_GDV(Name_Network_Bin[0], save_directory)

# Correlations:

Final_Correlation_Matrix_Bin = pd.DataFrame()
Final_Correlation_Matrix_Prob = pd.DataFrame()

Names = []
    
for repet, nodes, density, distribution, network , (mu,eta), Emp_option in product(Repetitions, Nodes, Densities, Distribution, Network, ETA_MU, emp_option):
    
    if distribution == "Unif":
        option = "Uniform" 
    if distribution == "Beta":
        option = "Beta"
    if distribution == "Empirical":
        option = "Empirical" + Emp_option
    
    Network_Name = [str(option)+"Eta_" + str(eta) + "_Mu_" + str(mu)+  "_"  + str(nodes) + "_Repetition_" + str(repet) + "_" + 
                    str(density) + str(network) + ".txt"]
    
    Names += Network_Name
    
    # Prepare the correlation Matrix:
    
    Network_Name   = ["UniformBA"] 
    
    Correlation_Distance_Functions.Prepare_Matrix_GDV_Corr(Network_Name[0],save_directory)
    Corr_Prob, Corr_Bin = Correlation_Distance_Functions.Final_Matrix_Corr(Network_Name[0],save_directory)
    
    # Add the name of the Network:
    
    if distribution == "Unif":
        i_eta = 0
        i_mu = 0
    if distribution == "Beta":
        i_eta = eta
        i_mu = mu
    if distribution == "Empirical":
        i_eta = 1
        i_mu = 1
        
    Corr_Prob["Network"] = str(network)
    Corr_Prob["Nodes"] = str(nodes)
    Corr_Prob["Densities"] = str(density)
    Corr_Prob["Distribution"] = str(distribution)
    Corr_Prob["Repetitions"] = str(repet)
    Corr_Prob["Variance"] = str(i_eta)
    Corr_Prob["Mean"] = str(i_mu)
    
    Corr_Bin["Network"] = str(network)
    Corr_Bin["Nodes"] = str(nodes)
    Corr_Bin["Densities"] = str(density)
    Corr_Bin["Distribution"] = str(distribution)
    Corr_Bin["Repetitions"] = str(repet)
    Corr_Bin["Variance"] = str(i_eta)
    Corr_Bin["Mean"] = str(i_mu)
    
    if distribution == "Empirical":
        
        Corr_Prob["Real"] = str(Emp_option)
        Corr_Bin["Real"] = str(Emp_option)
    
    # Add line to the final Data Frame of Results:
    
    Final_Correlation_Matrix_Prob = Final_Correlation_Matrix_Prob.append(Corr_Prob, sort = False)
    Final_Correlation_Matrix_Bin =  Final_Correlation_Matrix_Bin.append(Corr_Bin, sort = False)

# Save the results:
    
Final_Correlation_Matrix_Prob.to_csv(save_directory + "Final_Correlation_Prob.txt", header = True, index = False)
Final_Correlation_Matrix_Bin.to_csv(save_directory + "Final_Correlation_Bin.txt", header = True, index = False)    
    
# Once this is finished, the GCD and pGCD is obtained, the next step is to calculate the distances.