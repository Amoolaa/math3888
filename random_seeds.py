from dm_algo import dm_algo
import utils
import pandas as pd
import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt

# Plotting figure 1(a)

selected_disease = "lymphoma"
weights = np.linspace(0, 5, num = 100)
df = pd.DataFrame(index = weights)

disease_proteins = utils.get_disease_proteins(selected_disease)
G = nx.read_edgelist("data/ppin.txt",comments="#",nodetype=str)
GC = G.subgraph(max(nx.connected_components(G), key=len)).copy()

# Finding seeds present in giant component
seeds = []
for protein in disease_proteins:
  if protein in GC:
    seeds.append(protein)

n = len(seeds)

# seed degree dictionary
ks = utils.get_seed_degree(GC, seeds)

# iterate over many values for seed weight
num_seeds_kept = []
for w_s in weights:
  dm = dm_algo(GC, seeds, ks, w_s)
  seeds_kept = set(dm).intersection(set(seeds))
  num_seeds_kept.append(len(seeds_kept))

df[selected_disease] = num_seeds_kept


# Random Selection
seeds = random.sample(GC.nodes(), n)

# seed degree dictionary
ks = utils.get_seed_degree(GC, seeds)

# iterate over many values for seed weight
num_seeds_kept = []
for w_s in weights:
  dm = dm_algo(GC, seeds, ks, w_s)
  seeds_kept = set(dm).intersection(set(seeds))
  num_seeds_kept.append(len(seeds_kept))


df["random"] = num_seeds_kept


# Bad selection
cc = nx.centrality.closeness_centrality(GC).items().sort(key=lambda x:x[1])
seeds = [n for (n, v) in cc[:n]]

# seed degree dictionary
ks = utils.get_seed_degree(GC, seeds)

# iterate over many values for seed weight
num_seeds_kept = []
for w_s in weights:
  dm = dm_algo(GC, seeds, ks, w_s)
  seeds_kept = set(dm).intersection(set(seeds))
  num_seeds_kept.append(len(seeds_kept))


df["bad"] = num_seeds_kept



df.plot.line(color=["green", "pink", "brown"])
plt.title("Comparison on Random Seed Selection")
plt.ylabel("Number of Seed Proteins Included")
plt.xlabel("Seed Weight")
plt.legend([selected_disease, "random", "bad"])

plt.show()