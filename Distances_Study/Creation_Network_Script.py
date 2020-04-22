'''

This Script contains all the stepts to create a network with different probabilities distribution

'''

# Packages:

import numpy as np
import scipy.stats
import pandas as pd
import seaborn as sns
import networkx as nx
from scipy.stats import uniform
import matplotlib.pyplot as plt


# Network Creation:

# Functions:

# Function to choose the Radius for Rangom Geometric:
# Basically is just checking with which raidus we get what we want.

def Radius_Random_Geometric(Nodes, tries = 100):
    
    for r in np.linspace(0, 1,tries):
        G = nx.random_geometric_graph(Nodes, radius=r, dim=3)
        rho = nx.density(G)
        print(f"radius  = {r  :.3f}")
        print(f"density = {rho:.3f}")
        print()
        
# Function to Create Radom Geometric:

def get_Hyper_Geometric(Node, Radius):
    
    G = nx.random_geometric_graph(Node, Radius, dim=3)
    return(G)

# Function to Create a Network with BA model controling the density:

def get_m_for_barabasi_albert(Nodes,density):
    
    Edges = int(round(Nodes/2 - np.sqrt(Nodes**2/4 - Nodes*(Nodes-1)*density/2)))
    G = nx.barabasi_albert_graph(Nodes,Edges)
    
    return(G)
    
# Function to Create Erdos-Renyi Graph:
    
def get_Erdos_Renyi(Nodes, density):
    
    G = nx.erdos_renyi_graph(Nodes, density)
    return(G)
       
# Function to Add the Probabilities Based on Distributions:
    
    
def Apply_Prob(Network, save_directory, Network_Name, ETA, MU, empirical = False, option = "PPI"):
    
    n = Network.number_of_edges()
    start = 0
    width = 1
    
    # Uniform:
    
    data_Uniform = uniform.rvs(size=n, loc = start, scale=width)
    
    # Beta:
    
    if empirical == False:
    
        mu  = MU # Mean
        eta = ETA # Normalized Variance (0 to 1)  0.2
        
        var = eta*mu*(1-mu)
        nu = mu*(1-mu)/var-1
        alpha = nu*  mu
        beta  = nu*(1-mu)
        
        data_Beta = scipy.stats.beta.rvs(alpha, beta, size= n )
    
    # Rename the Network:
    
    Network = nx.to_pandas_edgelist(Network)
    x1 = 1
    
    Network.columns = ["node1", "node2"]
    
    # Uniform:
    
    Network["weight"] = data_Uniform
    Network_Uniform_Prob = Network.copy()
    Network_Uniform_Bin =  Network_Uniform_Prob.copy()
    Network_Uniform_Bin.weight =  np.repeat(x1,len(Network))
    
    Network_Uniform_Prob = nx.from_pandas_edgelist(Network_Uniform_Prob,'node1','node2','weight')
    Network_Uniform_Bin  = nx.from_pandas_edgelist(Network_Uniform_Bin,'node1','node2','weight')
    
    # Empirical:
    
    if empirical == True:
        
        # Selection of the distribution and getting the probabilities from each of them:
        
        if option == "PPI":
            file = pd.read_csv("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Topology_Probabilities/PPI/Work_Network_Prob_Low_PPI.txt",
                               sep = " ", header = None)
            probs = np.random.choice(file[2] , size= n , replace=True, p=None)
        if option == "GI":
            file = pd.read_csv("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Real_Distributions/Empirical_GI/Work_Network_Prob_Network_Medium_Confidence_Named.txt",
                               sep = " ", header = None)
            probs = np.random.choice(file[2] , size= n , replace=True, p=None)
        if option == "CoEx":
            file = pd.read_csv("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances/Tests/Real_Distributions/Empirical_CoEx/Work_Network_Prob_Coex_Net_Original_Medium_Rank_5%.txt",
                               sep = " ", header = None)
            probs = np.random.choice(file[2] , size= n , replace=True, p=None)
            
        Network["weight"] = probs
        Network_Real_Prob = Network.copy()
        Network_Real_Bin =  Network_Real_Prob.copy()
        Network_Real_Bin.weight =  np.repeat(x1,len(Network))
    
        Network_Real_Prob = nx.from_pandas_edgelist(Network_Real_Prob,'node1','node2','weight')
        Network_Real_Bin  = nx.from_pandas_edgelist(Network_Real_Bin,'node1','node2','weight')
        
    # Beta:
    
    if empirical == False:
    
        Network["weight"] = data_Beta
        Network_Beta_Prob = Network.copy()
        Network_Beta_Bin =  Network_Beta_Prob.copy()
        Network_Beta_Bin.weight =  np.repeat(x1,len(Network))
    
        Network_Beta_Prob = nx.from_pandas_edgelist(Network_Beta_Prob,'node1','node2','weight')
        Network_Beta_Bin  = nx.from_pandas_edgelist(Network_Beta_Bin,'node1','node2','weight')
    
    # Rename:
    
    Genes = Network_Uniform_Prob.nodes()
    gene2int = {gene:int for int,gene in enumerate(np.random.permutation(Genes),1)}
    gene2int_df = pd.DataFrame({'gene':list(gene2int.keys()), 'int':list(gene2int.values())})
    
    gene2int_df.to_csv(save_directory +"_" + "Gene_Names" + Network_Name)
    
    # Relabel:
    
    Network_Uniform_Prob = nx.relabel_nodes(Network_Uniform_Prob, gene2int)
    Network_Uniform_Bin  = nx.relabel_nodes(Network_Uniform_Bin, gene2int)
    
    # Save Beta:
    
    if empirical == False:
        
        Network_Beta_Prob    = nx.relabel_nodes(Network_Beta_Prob, gene2int)
        Network_Beta_Bin     = nx.relabel_nodes(Network_Beta_Bin, gene2int)
        
        with open(save_directory + 'Work_Network' + "_" + "Prob_Beta" + Network_Name, 'w') as f:
            for n1,n2 in sorted(map(sorted,Network_Beta_Prob.edges())):
                f.write(f"{n1} {n2} {Network_Beta_Prob[n1][n2]['weight']}\n")
         
        with open(save_directory + 'Work_Network' + "_" + "Bin_Beta" + Network_Name, 'w') as fl:
            for n1,n2 in sorted(map(sorted,Network_Beta_Bin.edges())):
                fl.write(f"{n1} {n2} {Network_Beta_Bin[n1][n2]['weight']}\n")
    
    # Save Empirical:
                    
    if empirical == True:
        
        Network_Real_Prob    = nx.relabel_nodes(Network_Real_Prob, gene2int)
        Network_Real_Bin     = nx.relabel_nodes(Network_Real_Bin, gene2int)
        
        with open(save_directory + 'Work_Network' + "_" + "Prob_Empirical"+ option + Network_Name, 'w') as f:
            for n1,n2 in sorted(map(sorted,Network_Real_Prob.edges())):
                f.write(f"{n1} {n2} {Network_Real_Prob[n1][n2]['weight']}\n")
         
        with open(save_directory + 'Work_Network' + "_" + "Bin_Empirical"+ option  + Network_Name, 'w') as fl:
            for n1,n2 in sorted(map(sorted,Network_Real_Bin.edges())):
                fl.write(f"{n1} {n2} {Network_Real_Bin[n1][n2]['weight']}\n")
        
    # Save Uniform:
    
    with open(save_directory + 'Work_Network' + "_" + "Prob_Uniform"  + Network_Name, 'w') as f:
        for n1,n2 in sorted(map(sorted,Network_Uniform_Prob.edges())):
            f.write(f"{n1} {n2} {Network_Uniform_Prob[n1][n2]['weight']}\n")
         
    with open(save_directory + 'Work_Network' + "_" + "Bin_Uniform"  +Network_Name, 'w') as fl:
        for n1,n2 in sorted(map(sorted,Network_Uniform_Bin.edges())):
            fl.write(f"{n1} {n2} {Network_Uniform_Bin[n1][n2]['weight']}\n")
         

# Function to PLot the Distributions:
    
def Plot_Distribution(probabilities, save_directory, Network_Name, option):
    
        probabilities = np.array(probabilities).astype(np.float)

        plt.figure(figsize=(9,5))
        plt.rcParams["axes.labelsize"] = 15

        a = sns.distplot(probabilities, kde=False)
        
        a.axes.set_title("Distribution", fontsize=18)
        a.set(xlabel='# Prob')
        a.tick_params(labelsize=12)
        
        fig = a.get_figure()
        fig.savefig(save_directory + "Distribution_" + str(Network_Name) + "_" + str(option) + ".png")
        
        plt.show()    
        


        

























        










