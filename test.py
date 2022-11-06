from dm_algo import dm_algo, random_neighbours
import utils
import pandas as pd
import networkx as nx
import numpy as np
import random

disease_proteins = utils.get_disease_proteins("data/disease_protein", 15)
G = nx.read_edgelist("data/ppin.txt",comments="#",nodetype=str)
GC = G.subgraph(max(nx.connected_components(G), key=len)).copy()

# Finding seeds present in giant component
seeds = []
for prot in disease_proteins:
  if prot in GC:
    seeds.append(prot)

# Randomly deleting a quarter of the seeds
arr = seeds.copy()
random.shuffle(arr)

deleted = []
for i in range(0, len(seeds) // 4):
  deleted.append(arr.pop())

print(f"deleted: {deleted}")

print(diamond.py())