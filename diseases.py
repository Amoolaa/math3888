from dm_algo import dm_algo, random_neighbors, best_neighbors
import utils
import pandas as pd
import networkx as nx
import random
import numpy as np

df = pd.DataFrame(columns = ["dm_algo", "best", "random"])

for i in [2, 7, 10, 15, 16]:
  disease_proteins = utils.get_disease_proteins("data/disease_protein", i)
  G = nx.read_edgelist("data/ppin.txt",comments="#",nodetype=str)
  GC = G.subgraph(max(nx.connected_components(G), key=len)).copy()

  dm_algo_results = []
  random_neighbors_results = []
  best_neighbors_results = []

  # Finding seeds present in giant component
  seeds = []
  for prot in disease_proteins:
    if prot in GC:
      seeds.append(prot)

  k = {}

  for node in GC.nodes():
    # Find number of links to seeds
    ks = 0
    for neighbour in GC.neighbors(node):
      if neighbour in seeds:
        ks += 1

    k[node] = ks/GC.degree(node)


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

    new = dm_algo(GC, arr, 0.95, k, n = n_removed)

    if new == []:
      dm_algo_results.append(0)
    else:
      correct = 0
      false = 0
      for n in new:
        if n in deleted:
          correct += 1
        else:
          false += 1

      dm_algo_val = correct/len(new)
      dm_algo_results.append(dm_algo_val)


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

    # random neighbours here
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
  
  dm_algo_avg = sum(dm_algo_results)/num_trials
  best_avg = sum(best_neighbors_results)/num_trials
  random_avg = sum(random_neighbors_results)/num_trials

  print(f"dm_algo -> val: {dm_algo_avg}, disease no: {i}, num_seeds: {len(seeds)}")
  print(f"best_neighbors -> val: {best_avg}, disease no: {i}, num_seeds: {len(seeds)}")
  print(f"random_neighbors -> val: {random_avg}, disease no: {i}, num_seeds: {len(seeds)}")

  df.loc[len(df.index)] = [dm_algo_avg, best_avg, random_avg]

  input("next disease\n")

print(df)

df.to_csv("diseases.csv", index=False)