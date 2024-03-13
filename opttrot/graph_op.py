from typing import *

import networkx as nx
import pandas as pd
import numpy as np

from opttrot import Hamiltonian



def get_binary_graph(h:Hamiltonian):
    edge_df = h.graph_edge
    
    #edge_df["commute"] = (edge_df["commute"]-1).abs()
    edge_df= edge_df.loc[edge_df["commute"] == 1]
    
    G = nx.from_pandas_edgelist(
        edge_df,
        source="source",
        target="target",
        edge_attr="commute"
    )
    return G

def get_binary_H(h:Hamiltonian, mu0:float):
    # D-Wave routine
    edge_df = h.graph_edge
    N = h.qubit_num
    
    mu1 = N*mu0 +1
    edge_df["commute"] = (edge_df["commute"]-1).abs() # reverse
    edge_df= edge_df.loc[edge_df["commute"] == 1] # get anti commute pairs
    
    edges = edge_df.set_index(["source", "target"])["commute"].multiply(mu1).to_dict()
    nodes = h.local_decomposition.loc[h.local_decomposition["Pstring"]!=N*"I"]["Pstring"]
    H_dwave = {**edges}
    for n in nodes:
        H_dwave[n] = -mu0
    return H_dwave

def mu_cal(mu, N):
    mu2 = mu
    mu0 = 0.5*N*(N-1)*mu2*3 + 1
    mu1 = N*mu0 + 1
    return mu0, mu1, mu2

def get_basis_weight_graph(
    h:Hamiltonian, 
    mus:Tuple[float, float, float]):
    
    edge_df = h.graph_edge
    N = h.qubit_num
    mu0, mu1, mu2 = mus
    
    commute_column = (edge_df["commute"]-1).abs()
    basis_column = edge_df["basis_wieght"]
    
    edge_df["weight"] = mu2*basis_column + mu1*commute_column
    
    nodes = h.local_decomposition["Pstring"].loc[h.local_decomposition["Pstring"] !=N*"I"].unique()
    np_data = np.vstack([nodes, nodes.size*[-mu0]]).T
    nodes_graph = pd.DataFrame(np_data, columns=["node", "weight"])
    
    G = nx.from_pandas_edgelist(
        edge_df, 
        source="source", target="target",
        edge_attr="weight")
    G.add_nodes_from((n, dict(d)) for n, d in nodes_graph.iterrows())
    return G

def get_basis_weight_H(
    # D-Wave
    h:Hamiltonian, 
    mus:Tuple[float, float, float]):
    
    edge_df = h.graph_edge
    N = h.qubit_num
    mu0, mu1, mu2 = mus
    
    commute_column = (edge_df["commute"]-1).abs()
    basis_column = edge_df["basis_wieght"]
    
    edge_df["weight"] = mu2*basis_column + mu1*commute_column
    
    edges = edge_df.set_index(["source", "target"])["weight"].to_dict()
    nodes = h.local_decomposition.loc[h.local_decomposition["Pstring"]!=N*"I"]["Pstring"]
    H_dwave = {**edges}
    for n in nodes:
        H_dwave[n] = -mu0
    return H_dwave
