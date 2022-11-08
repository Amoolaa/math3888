from dm_algo import dm_algo, random_neighbors, best_neighbors, diamond
import utils
import pandas as pd
import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt
import os

# Tests performance of algorithms on recovering seeds. Plots Figure 2

selected_diseases = [os.fsdecode(file)[:-4] for file in os.listdir("data/diseases")]

# dataframe containing results for each algorithm
df = pd.DataFrame(columns = ["diamond", "dm_algo", "best", "random"], index = selected_diseases)

for disease_name in selected_diseases:
  disease_proteins = utils.get_disease_proteins(disease_name)
  G = nx.read_edgelist("data/ppin.txt",comments="#",nodetype=str)
  GC = G.subgraph(max(nx.connected_components(G), key=len)).copy()

  diamond_results = []
  dm_algo_results = []
  random_neighbors_results = []
  best_neighbors_results = []

  # Finding seeds present in giant component
  seeds = []
  for prot in disease_proteins:
    if prot in GC:
      seeds.append(prot)
  
  # seed degree dictionary
  ks = utils.get_seed_degree(GC, seeds)

  # Number of trials for each algorithm
  num_trials = 100
  for j in range(num_trials):

    # randomising
    arr = seeds.copy()
    random.shuffle(arr)

    deleted = []

    # Removing half
    n_removed = len(seeds) // 2
    for h in range(0, n_removed):
      deleted.append(arr.pop())

    # Running diamond
    new = diamond(GC, arr, n_removed)

    correct = 0
    false = 0
    for n in new:
      if n in deleted:
        correct += 1
      else:
        false += 1

    diamond_val= correct/len(new)
    diamond_results.append(diamond_val)


    # Running dm algorithm
    dm = dm_algo(GC, arr, ks, 3)

    # Finding set of new nodes
    new = [n for n in dm if n not in arr]

    # sorting set of new nodes and taking top
    l = [(n, ks[n]) for n in new]
    l.sort(reverse=True, key=lambda x:x[1])
    best_n = [n for n, v in l[:n_removed]]

    if best_n == []:
      dm_algo_val = correct/len(best_n)
    else:
      correct = 0
      false = 0
      for n in best_n:
        if n in deleted:
          correct += 1
        else:
          false += 1
      dm_algo_val = correct/len(best_n)
    
    dm_algo_results.append(dm_algo_val)

    # Running best neighbours algorithm
    new = best_neighbors(GC, arr, n_removed)

    correct = 0
    false = 0
    for n in new:
      if n in deleted:
        correct += 1
      else:
        false += 1

    best_neighbors_val = correct/len(new)
    best_neighbors_results.append(best_neighbors_val)

    # random neighbours
    new = random_neighbors(GC, arr, n_removed)
    correct = 0
    false = 0
    for node in new:
      if node in deleted:
        correct += 1
      else:
        false += 1

    random_neighbors_val = correct/(correct + false)
    random_neighbors_results.append(random_neighbors_val)
  
  diamond_avg = sum(diamond_results)/num_trials
  dm_algo_avg = sum(dm_algo_results)/num_trials
  best_avg = sum(best_neighbors_results)/num_trials
  random_avg = sum(random_neighbors_results)/num_trials

  print(f"diamond -> val: {diamond_avg}, disease: {disease_name}, num_seeds: {len(seeds)}")
  print(f"dm_algo -> val: {dm_algo_avg}, disease: {disease_name}, num_seeds: {len(seeds)}")
  print(f"best_neighbors -> val: {best_avg}, disease: {disease_name}, num_seeds: {len(seeds)}")
  print(f"random_neighbors -> val: {random_avg}, disease: {disease_name}, num_seeds: {len(seeds)}")
  print("==============")

  df.loc[disease_name] = [diamond_avg, dm_algo_avg, best_avg, random_avg]

df.to_csv("results.csv", sep=" ")

selected = ["blood coagulation disorders", "blood platelet disorders", "lipid metabolism disorders", "nutritional and metabolic diseases", "crohn disease"]

# Plotting
f = plt.figure(1, figsize=(50,20))
ax = f.add_subplot(1,1,1)
df.plot(kind='bar', ax = ax)
plt.ylabel("Average Prediction Accuracy")
plt.xlabel("Diseases")
ax.set_xticklabels(selected, rotation=0, fontsize=6)
plt.legend(["DIAMOnD", "PCSTP-DM", "Best Neighbours", "Random Neighbours"])
plt.title("Performance of Disease Module finding methods")
plt.show()