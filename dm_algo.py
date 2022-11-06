import networkx as nx
import pandas as pd
import utils
import pcst_fast
import random

# GC = giant component of network. If n != -1, return best n nodes.
def dm_algo(GC, seeds, t, k, n = -1):

  # pcst_fast requires integers as vertices
  node_id = {node: i for i, node in enumerate(GC.nodes)}
  id_node = {v: k for k, v in node_id.items()}
  edge_id = [(node_id[u], node_id[v]) for u, v in GC.edges]

  # Configuring prize and edge costs
  prizes = []
  for node in GC.nodes():
    if node in seeds:
      prizes.append(10)
    else:
      prizes.append(k[node])
  edge_cost = [1 for i in range(GC.number_of_edges())]

  # Finding PCSTs
  counts = {}
  for seed in seeds:
    root = node_id[seed]
    vertices, edges = pcst_fast.pcst_fast(edge_id, prizes, edge_cost, root, 1, "gw", 0)
    for v in vertices:
      if v in counts:
        counts[v] += 1
      else:
        counts[v] = 1

  final = [id_node[node] for node, count in counts.items() if count/len(seeds) >= t]
  new = [node for node in final if node not in seeds]

  # Return best n nodes
  if n != -1:
    l = [(node, k[node]) for node in new]
    l.sort(reverse=True, key=lambda x:x[1])
    return [node for (node, val) in l[:n]]
  else:
    return new

# find n random seed neighbours
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


# find n best neighbours
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