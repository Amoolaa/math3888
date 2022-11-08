import networkx as nx
import pandas as pd
import utils
import pcst_fast
import random

# Contains algorithms discussed

# PCSTP-DM
def dm_algo(G, seeds, k, w_s):

  # pcst_fast requires integers as vertices
  node_id = {node: i for i, node in enumerate(G.nodes)}
  id_node = {v: k for k, v in node_id.items()}
  edge_id = [(node_id[u], node_id[v]) for u, v in G.edges]

  # Configuring prize and edge costs
  prizes = []
  for node in G.nodes():
    if node in seeds:
      prizes.append(max(k[node], w_s))
    else:
      prizes.append(k[node])
  edge_cost = [1 for i in range(G.number_of_edges())]

  root = -1

  vertices, edges = pcst_fast.pcst_fast(edge_id, prizes, edge_cost, root, 1, "gw", 0)

  dm = [id_node[v] for v in vertices]
  return dm


# Random Neighbours
def random_neighbors(GC, seeds, n):
  seed_neighbors = []
  for seed in seeds:
    seed_neighbors.extend(list(GC.neighbors(seed)))

  seed_neighbors = list(set(seed_neighbors))

  final = []

  for i in range(0, n):
    random.shuffle(seed_neighbors)
    new = seed_neighbors.pop()
    final.append(new)
    seed_neighbors.extend(list(GC.neighbors(new)))
    seed_neighbors = list(set(seed_neighbors))

  return final


# Best Neighbours
def best_neighbors(GC, seeds, n):

  temp = seeds.copy()

  seed_neighbors = {}
  for seed in temp:
    for neighbor in GC.neighbors(seed):
      if neighbor not in temp:
        # finding number of links to other seeds
        ks = 0
        for node in GC.neighbors(neighbor):
          if node in temp:
            ks += 1
      
        seed_neighbors[neighbor] = ks/GC.degree(neighbor)

  final = []

  for i in range(0, n):
  

    best, ks = sorted(seed_neighbors.items(), key=lambda x:x[1]).pop()

    seed_neighbors.pop(best)

    final.append(best)
    temp.append(best)

    # updating neighbors
    for node in GC.neighbors(best):
      # Find number of links to seeds
      if node not in temp:
        ks = 0
        for neighbor in GC.neighbors(node):
          if neighbor in temp:
            ks += 1

        seed_neighbors[node] = ks/GC.degree(node)

    final.append(best)

  return final

# Diamond
def diamond(GC, seeds, n):

  temp = seeds.copy()

  seed_neighbors = {}
  for seed in temp:
    for neighbor in GC.neighbors(seed):
      if neighbor not in temp:
        # finding number of links to other seeds
        ks = 0
        for node in GC.neighbors(neighbor):
          if node in temp:
            ks += 1
        seed_neighbors[neighbor] = utils.pvalue(ks, GC.degree(neighbor), GC.number_of_nodes(), len(temp))

  final = []

  for i in range(0, n):
  

    best, ks = sorted(seed_neighbors.items(), key=lambda x:x[1], reverse=True).pop()

    seed_neighbors.pop(best)

    final.append(best)
    temp.append(best)

    # updating neighbors
    for node in GC.neighbors(best):
      # Find number of links to seeds
      if node not in temp:
        ks = 0
        for neighbor in GC.neighbors(node):
          if neighbor in temp:
            ks += 1

        seed_neighbors[node] = utils.pvalue(ks, GC.degree(node), GC.number_of_nodes(), len(temp))

    final.append(best)

  return final