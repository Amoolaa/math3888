import pandas as pd
from math import comb

# Helper methods

def connectivity_significance(ks, k, N, s0):
  return comb(s0, ks) * comb(N-s0, k-ks) * 1/comb(N, k)

def pvalue(ks, k, N, s0):
    p = 0.0
    for ki in range(ks, k + 1):
        if ki > s0:
            break
        prob = connectivity_significance(ki, k, N, s0)
        # print prob
        p += prob

    if p > 1:
        return 1
    else:
        return p

def p_dict(GC, seeds):
  # Finding connectivity significance for every node
  N = GC.number_of_nodes()
  s0 = len(seeds)
  p = {}
  for node in GC.nodes():
    # Find number of links to seeds
    ks = 0
    for neighbour in GC.neighbors(node):
      if neighbour in seeds:
        ks += 1

    p[node] = pvalue(ks, GC.degree(node), N, s0)
  
  return p

# Helper to get disease proteins from data file with strange format quickly
def get_disease_proteins(disease_name):
  f = open("data/diseases/" + disease_name + ".txt", "r")
  proteins = f.read().split('\n')
  return proteins

def get_seed_degree(G, seeds):
  k = {}

  for node in G.nodes():
    # Find number of links to seeds
    ks = 0
    for neighbour in G.neighbors(node):
      if neighbour in seeds:
        ks += 1
    k[node] = ks/G.degree(node)
  
  return k
