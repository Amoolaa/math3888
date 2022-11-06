from dm_algo import dm_algo, random_neighbours
import utils
import pandas as pd
import networkx as nx
import random
import numpy as np

df = pd.DataFrame(columns=["alpha", "dm_algo", "random_neighbours", "diamond"])

disease_proteins = utils.get_disease_proteins("data/disease_protein", 1)
G = nx.read_edgelist("data/ppin.txt",comments="#",nodetype=str)
GC = G.subgraph(max(nx.connected_components(G), key=len)).copy()

# Finding seeds present in giant component
seeds = []
for prot in disease_proteins:
  if prot in GC:
    seeds.append(prot)

num_trials = 5

for alpha in np.logspace(start = 0, stop = 4, num=5):
  dm_algo_results = []
  random_neighbour_results = []
  diamond_results = []
  for i in range(0, num_trials):

    # randomising
    arr = seeds.copy()
    random.shuffle(arr)

    deleted = []
    n_removed = len(seeds) // 4
    for i in range(0, n_removed):
      deleted.append(arr.pop())
    
    final = dm_algo(GC, arr, alpha, 0.95)
    

    correct = 0
    false = 0
    for node in final:
      if node in deleted:
        correct += 1
      elif node not in seeds:
        false += 1
    
    dm_algo_val = correct/(correct + false)
    dm_algo_results.append(dm_algo_val)
    print(f"dm_algo: {dm_algo_val}")

    # random neighbours here
    final = random_neighbours(GC, arr, correct + false)
    for node in final:
      if node in deleted:
        correct += 1
      elif node not in seeds:
        false += 1

    random_neighbour_val = correct/(correct + false)
    random_neighbour_results.append(random_neighbour_val)
    print(f"rn: {random_neighbour_val}")

  print("=============")
  df.loc[len(df.index)] = [alpha, sum(dm_algo_results)/num_trials, sum(random_neighbour_results)/num_trials, sum(diamond_results)/num_trials]


df.to_csv("results.csv", index = False)