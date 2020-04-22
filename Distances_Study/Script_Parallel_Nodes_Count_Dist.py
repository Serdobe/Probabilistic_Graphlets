import os
import subprocess 
import numpy as np
import pandas as pd
import networkx as nx  
import multiprocessing
from multiprocessing import Pool

# DESCRIPTION:

'''

This library includes all the function which are necessary to compute the orbits
and the distance from an edge list with gene names. All is prepared to be paralelized
per nodes. 


'''

# This file is for the orbit calculation but parallelizing the nodes:

def Creat_Network(Network_Name, save_directory):
    
    original_pd = pd.read_csv(Network_Name, sep = " ")
    original_pd.columns = ["node1", "node2", "weight"]
    x1 = 1
    
    original_pd_prob = original_pd
    original_pd_bin = original_pd.copy()
    original_pd_bin.weight =  np.repeat(x1,len(original_pd))
    
    original_pd_prob = nx.from_pandas_edgelist(original_pd_prob,'node1','node2','weight')
    original_pd_bin =  nx.from_pandas_edgelist(original_pd_bin,'node1','node2','weight')
    
    # with this we first create the names for both:
    
    Genes = original_pd_prob.nodes()

    gene2int = {gene:int for int,gene in enumerate(np.random.permutation(Genes),1)}
    gene2int_df = pd.DataFrame({'gene':list(gene2int.keys()), 'int':list(gene2int.values())})

    gene2int_df.to_csv(save_directory +"_" + "Gene_Names_" + Network_Name)
    
    # Now we create the two files of the Networks:
    
    G_prob = nx.relabel_nodes(original_pd_prob, gene2int)
    G_bin = nx.relabel_nodes(original_pd_bin, gene2int)
    
    with open(save_directory + 'Work_Network' + "_" + "Prob_" + Network_Name, 'w') as f:
     for n1,n2 in sorted(map(sorted,G_prob.edges())):
         f.write(f"{n1} {n2} {G_prob[n1][n2]['weight']}\n")
         
    with open(save_directory + 'Work_Network' + "_" + "Bin_" + Network_Name, 'w') as fl:
     for n1,n2 in sorted(map(sorted,G_bin.edges())):
         fl.write(f"{n1} {n2} {G_bin[n1][n2]['weight']}\n")
         
# Function to call different times the Terminal:

def run_process(process):                                                             
    os.system(format(process)) 
    
# Function to divide the nodes in different chunks:
    
def divide_Nodes(Total_Number_Nodes, threads):
    
    n_cores = multiprocessing.cpu_count() * threads # Depending of the threads we split more or less
    block_size = int(Total_Number_Nodes/n_cores)
    
    data_splits = pd.DataFrame(columns=["From", "To"])
    
    for block in range(n_cores):
        
    
        if block == n_cores:
             From = block * block_size
             To   = Total_Number_Nodes
             
        else:
             From = block*block_size
             To   = (block+1)*block_size
             
        junk_vector =  pd.DataFrame({"From": [From],
                                     "To": [To]})
    
        data_splits = data_splits.append(junk_vector)
    data_splits.index = pd.Series(range(n_cores))
        
    return(data_splits)

# Function to run different times G-tries:    
     
def Compute_Orbits_Spliting_Nodes(Network_Name, Total_Number_Nodes, file_Directory, save_directory, Run_path, threads = 1):
    
    os.chdir("/home/sergio/workspace_Eclipse/Lucky_GOGO_Extra/Results/Testeos_Varios/Network_Distances") # Cambiado para las distancias
    
    # The groups are made taking into account CPUs and threads.
    
    groups = divide_Nodes(Total_Number_Nodes, threads)
    
    commands = pd.DataFrame(columns=["Command_G0, Command_G3", "Command_G4"])

    for sentence in range(len(groups)):
        
        sub_command_G0 = [str(Run_path) + "Lucky_GOGO_Node_Parallel -s 3" + " "  +
                             "-m gtrie undir3.gt " + " " +
                             " -g" + " " +str(file_Directory)  + str(Network_Name) +
                             " -p " +  str(save_directory) +"Total_G0_"+ "Split_" 
                             + str(sentence) + "_"  + str(Network_Name) + 
                             " -fp "+  str(save_directory) +
                             "Orbits_G0_" + "Split_" + 
                             str(sentence) + "_" + str(Network_Name)+ " -g0 "]
        
        sub_command_G3 = [str(Run_path) + "Lucky_GOGO_Node_Parallel -s 3" + " "  +
                             "-m gtrie undir3.gt " + " " +
                             " -g" + " " +str(file_Directory)  + str(Network_Name) +
                             " -p " +  str(save_directory) +"Total_G3_" + "Split_" 
                             + str(sentence) + "_"   + str(Network_Name)+ 
                             " -fp "+  str(save_directory) +
                             "Orbits_G3_"+ "Split_" 
                             + str(sentence) + "_"  + str(Network_Name)+ " -paral " +
                             " -From " + str(groups.From[sentence]) + " -To " + str(groups.To[sentence])]
        
        sub_command_G4 = [str(Run_path) + "Lucky_GOGO_Node_Parallel -s 4" + " "  +
                             "-m gtrie undir4.gt " + " " +
                             " -g" + " " +str(file_Directory)  + str(Network_Name) +
                             " -p " +  str(save_directory) +"Total_G4_"+ "Split_" 
                             + str(sentence) + "_"  + str(Network_Name) + 
                             " -fp "+  str(save_directory) +
                             "Orbits_G4_"+ "Split_" 
                             + str(sentence) + "_"  + str(Network_Name)+ " -paral" + 
                             " -From " + str(groups.From[sentence]) + " -To " + str(groups.To[sentence])]
        
        result_ite = pd.DataFrame({"Command_G0": sub_command_G0[0], 
                                   "Command_G3": sub_command_G3[0], 
                                   "Command_G4": sub_command_G4[0]}, index=[sentence])
    
        commands = commands.append(result_ite)
    
    # Multiprocess: Is necessary to decide the order of the computing to not saturate the cores:
    
    # First the G0 and G3 (should both finish quicK):
    
    processes1  = []
    
    for split in range(len(commands)):
        
        processes1.append(commands.Command_G3[split])
    
    # Now we split the process by CPUs:
      
    n_cores = multiprocessing.cpu_count() 
    
    pool1 = Pool(processes  = n_cores)
    pool1.map(run_process, processes1)
    pool1.close()
    pool1.join() 
    
    processes2  = []
    
    print("Finish with G3")
    
    # For G4:

    for split in range(len(commands)):
        
          processes2.append(commands.Command_G4[split])
          
    pool2 = Pool(processes  = n_cores)
    pool2.map(run_process, processes2)
    pool2.close()
    pool2.join() 
    
    print("Finish with G4")

    # Finally only one for the G0:
    
    processes3 = []
    processes3.append(commands.Command_G0[0])
    
    pool3 = Pool(processes  = 1) 
    pool3.map(run_process, processes3) 
    pool3.close()
    pool3.join() 
    
    print("End of G0")
        
    # Now we have to join the results of each chunk:
    
    Sum_Orbits_3_Nodes(Network_Name, save_directory, len(commands))
    Sum_Orbits_4_Nodes(Network_Name, save_directory, len(commands))
    
    # Finally the splited data is deleted from the folder:
    
    for i in range(len(commands)):
        
        os.remove(str(save_directory) + "Orbits_G3_"+ "Split_" + str(i) + "_" + str(Network_Name))
        os.remove(str(save_directory) + "Orbits_G4_"+ "Split_" + str(i) + "_" + str(Network_Name))

        
# Funtion to Sum The Results from the Split Counting (G3):  
    
def Sum_Orbits_3_Nodes(Network_Name, save_directory, splits):
    
    d={}
    
    for x in range(splits):
            
            file = [str(save_directory) + "Orbits_G3_"+ "Split_" + str(x) + 
            "_"  + str(Network_Name)]
            
            d["string{0}".format(x)]= pd.read_csv(file[0], sep = " ", header = None,
              skiprows = [0])
    
    final = d["string{0}".format(0)].add(d["string{0}".format(1)], fill_value=0)
            
    for sums in range(2, splits):
        
        final = final.add(d["string{0}".format(sums)], fill_value = 0)
        
    final = final.drop( 0 , axis=1)
            
    final.to_csv(str(save_directory) + "_Orbits_G3_Sum_" + str(Network_Name), header = True, index = False)

    return("Sum Done")
    
    
# Funtion to Sum The Results from the Split Counting (G4): 
    
def Sum_Orbits_4_Nodes(Network_Name, save_directory, splits):
    
    d={}
    
    for x in range(splits):
            
            file = [str(save_directory) + "Orbits_G4_"+ "Split_" + str(x) + 
            "_"  + str(Network_Name)]
            
            d["string{0}".format(x)]= pd.read_csv(file[0], sep = " ", header = None,
              skiprows = [0])
    
    final = d["string{0}".format(0)].add(d["string{0}".format(1)], fill_value=0)
            
    for sums in range(2, splits):
        
        final = final.add(d["string{0}".format(sums)], fill_value = 0)
        
    final = final.drop( 0 , axis=1)
            
    final.to_csv(str(save_directory) + "_Orbits_G4_Sum_" + str(Network_Name), header = True, index = False)

    return("Sum Done")
    
# Function to prepare the results for distances:
        
def Prepare_Matrix_GDV(Network_Name, save_directory):
    
    G0 = pd.read_csv(str(save_directory) + "Orbits_G0_" + "Split_0_" + str(Network_Name),
                              sep = " ", header = None, skiprows=[0])
    
    G3 = pd.read_csv(str(save_directory) + "_Orbits_G3_Sum_" + str(Network_Name), 
                              sep = ",", header = None, skiprows=[0])
    
    G4 = pd.read_csv(str(save_directory) + "_Orbits_G4_Sum_" + str(Network_Name), 
                              sep = ",", header = None, skiprows=[0])
    
    GDV = pd.concat([G0, G3, G4], axis = 1)
    GDV.columns = ["NODE", "0", "1", "2" ,"3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9" ,"10" ,"11" ,"12" ,"13", "14"]
    
    GDV.to_csv(save_directory + "GDV_" + str(Network_Name), header = True, index = False)
    
    return(GDV)
    
def Prepare_To_Distance(GDV, Name_Network, save_directory):
    
    GDV = GDV.drop("NODE", axis = 1)
    np.savetxt(save_directory + "_Temporal_File" + str(Name_Network), np.array(GDV), 
               header =  f'{len(GDV)}' + ' ' +'15', fmt='%.7f')   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
    
    
    

