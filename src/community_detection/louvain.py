import networkx as nx
import community as community_louvain
import matplotlib.cm as cm
import numpy as np
import scipy as sp
from string import ascii_lowercase 

import csv
import matplotlib.pyplot as plt
font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 32}

plt.rc('font', **font)
G0 = nx.read_weighted_edgelist("4932.protein.links.v11.5.txt",comments="#",nodetype=str) #Delete the header first or this might not work!

print('number of nodes of G0:',G0.number_of_nodes())
#Remove the essential proteins from G0
with open("essential_yeast.csv") as f:
    essential = csv.reader(f)
    for row in essential:
        if "4932."+row[1] in G0.nodes:
            G0.remove_node("4932."+row[1])
print('number of nodes of G0:',G0.number_of_nodes())

# delete those edges with a combined score of <= threshold_score (small confidence)

#threshold_score = 950
#for edge in G0.edges: 
#    weight = list(G0.get_edge_data(edge[0],edge[1]).values())
#    if(weight[0] <= threshold_score):
#        G0.remove_edge(edge[0],edge[1])

partition = community_louvain.best_partition(G0)
number_of_communities = max(list(partition.values()))+1
print('# of partitions for Louvain modularity =',number_of_communities)
node_target = "4932.YGR241C"
print('The target protein YAP1802 belongs to community #',partition[node_target])
target_partition = []
for node in partition:
    if partition[node] == partition[node_target]:
        target_partition.append(node)
H = G0.subgraph(target_partition)
labels = {n:n for n in H}
pos = nx.spring_layout(H)
print("Position of YAP1802:", pos["4932.YGR241C"])
#Compute subgraph centrality of the community of YAP1802
H_SC = nx.subgraph_centrality(H)
#Sort the nodes in the community by subgraph centrality
H_SC = dict(sorted(H_SC.items(), key = lambda item: item[1], reverse = True))
H_SC_ls = []
for item in H_SC:
    H_SC_ls.append((item, H_SC[item]))
print(H_SC_ls[:5])
#Compute degree centrality of the community of YAP1802
H_DC = nx.degree_centrality(H)
#Sort the nodes in the community by degree centrality
H_DC = dict(sorted(H_DC.items(), key = lambda item: item[1], reverse = True))
H_DC_ls = []
for item in H_DC:
    H_DC_ls.append((item, H_DC[item]))
print(H_DC_ls[:5])
for i in range(0, len(H_SC_ls)):
    if H_SC_ls[i][0] == "4932.YGR241C":
        print("Subgraph centrality rank:", i)
    if H_DC_ls[i][0] == "4932.YGR241C":
        print("Degree centrality rank:", i)
#Find all cliques that YAP1802 is in
yap_cliques = []
yap_degree = H.degree("4932.YGR241C")
print("YAP1802 degree:", yap_degree)
print("Max degree:",H.degree(H_DC_ls[0]))
for clique in nx.find_cliques(H):
    if "4932.YGR241C" in clique:
        yap_cliques.append(clique)
        #print(clique)
#Length of cliques
clique_lengths = []
for clique in nx.find_cliques(H):
    clique_lengths.append(len(clique))
print(clique_lengths)
target_two = "4932.YLL001W"
nx.draw_networkx_nodes(H, pos, node_size = 30)
nx.draw_networkx_edges(H, pos, alpha=0.5)
nx.draw_networkx_labels(H,pos,labels,font_size=10)
plt.show()
