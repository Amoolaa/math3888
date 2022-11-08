from dm_algo import dm_algo
import utils
import pandas as pd
import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt

# Plotting Figure 1 (a)

selected_diseases = ["lupus erythematosus", "blood platelet disorders", "lymphoma", "arthritis rheumatoid", "sarcoma"]
weights = np.linspace(0, 5, num = 100)
df = pd.DataFrame(index = weights)

for disease_name in selected_diseases:
  disease_proteins = utils.get_disease_proteins(disease_name)
  G = nx.read_edgelist("data/ppin.txt",comments="#",nodetype=str)
  GC = G.subgraph(max(nx.connected_components(G), key=len)).copy()

  # Finding seeds present in giant component
  seeds = []
  for protein in disease_proteins:
    if protein in GC:
      seeds.append(protein)

  # seed degree dictionary
  ks = utils.get_seed_degree(GC, seeds)

  # iterate over many values for seed weight
  num_seeds_kept = []
  for w_s in weights:
    dm = dm_algo(GC, seeds, ks, w_s)
    seeds_kept = set(dm).intersection(set(seeds))
    num_seeds_kept.append(len(seeds_kept))


  df[disease_name] = num_seeds_kept

df.plot.line()
plt.title("Effect of Changing Seed Weight")
plt.ylabel("Number of Seed Proteins Included")
plt.xlabel("Seed Weight")
plt.legend(selected_diseases)

plt.show()