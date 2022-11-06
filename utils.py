import pandas as pd
from math import comb

def connectivity_significance(ks, k, N, s0):
  return comb(s0, ks) * comb(N-s0, k-ks) * 1/comb(N, k)

def pvalue(ks, k, N, s0):
    """
    -------------------------------------------------------------------
    Computes the p-value for a node that has kb out of k links to
    seeds, given that there's a total of s seeds in a network of N nodes.
    p-val = \sum_{n=kb}^{k} HypergemetricPDF(n,k,N,s)
    -------------------------------------------------------------------
    """
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
def get_disease_proteins(f, disease_no):
  f = open("data/id_symbol.txt", "r")
  id_to_symbol = {}
  for line in f.readlines():
    id, symbol = line.split(" ")
    id_to_symbol[id] = symbol.strip()
  f.close()

  disease_protein_df = pd.read_csv("data/disease_protein.csv", sep=":")

  disease_proteins = []

  for id in disease_protein_df.iloc[disease_no]["proteins"].split("/"):
    disease_proteins.append(id_to_symbol[id])

  return disease_proteins